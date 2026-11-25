#!/usr/bin/env python3
"""
bvbrc2fasta.py — Convert BV-BRC raw JSON to FASTA with QC filtering.

Overview
--------
This script reads two BV-BRC "raw" JSON dumps:

    meta_raw.json        Metadata records (one per genome)
    sequences_raw.json   Sequence records  (one per genome)

and produces:

    <output>.fasta            FASTA of sequences that PASS all QC filters
    <output>_excluded.fasta   FASTA of sequences that FAIL QC (with reason)
    <output>.log              Full stderr transcript (via TeeLogger)

The script automatically detects Influenza A (Alpha) vs Influenza B (Beta)
from the metadata and applies the correct:

    * segment ID -> segment name mapping
    * TYPE_FIELD rule (subtype for A, lineage for B)
    * expected segment lengths for length-based QC

Pipeline
--------
    1. Load metadata and sequence JSON (plain or gzipped).
    2. Build lookups keyed by genome_id:
         metadata_lookup : genome_id -> metadata record
         sequence_lookup : genome_id -> sequence string
       Duplicate accessions / genome_ids are dropped (first record wins).
    3. For each metadata record:
         a. Detect virus type (alpha / beta / unknown).
         b. Extract and validate required fields.
         c. Apply sequence QC (length vs. expected, ambiguous bases).
         d. Write to passed FASTA or excluded FASTA.
    4. Print a summary and total runtime to stderr (and the log).

Required fields
---------------
    accession | genome_id | host | geographic_location |
    strain    | segment   | length | TYPE_FIELD

TYPE_FIELD
----------
    Influenza A -> subtype, matching regex H<num>N<num>(pdm<num>)?
    Influenza B -> lineage, one of {Victoria, Yamagata}

Header formats
--------------
    Passed:
      >accession host|segmentid|segmentName|TYPE_FIELD|
       geographic_location|strain|length

    Excluded:
      >accession=...|host=...|segment_id=...|segment_name=...|
       type=...|geographic_location=...|strain=...|length=...|
       reason=...

    Note: genome_id is used internally for lookups but is NOT written
    to either FASTA header. The accession is the sole record identifier.

Notes
-----
    * Segment lengths are based on the influenza reference described in:
      https://pmc.ncbi.nlm.nih.gov/articles/PMC3074182/
    * Length QC uses the metadata-supplied genome_length compared against
      a fraction (min_len_pct / max_len_pct) of the expected segment length.
    * Ambiguous bases counted: N R Y S W K M B D H V (IUPAC ambiguity codes).

Author: <add your name here>
License: <add license here>
"""

import argparse
import gzip
import json
import os
import re
import sys
import time


# =============================================================================
# Segment ID mappings
# =============================================================================
# Influenza A and B share segment *names* but assign them to different
# segment *numbers*. The maps below translate the BV-BRC numeric segment ID
# (a string like "1", "2", ...) to the canonical influenza segment name.

# Influenza A (Alphainfluenzavirus)
#   1=PB2, 2=PB1, 3=PA, 4=HA, 5=NP, 6=NA, 7=MP, 8=NS
SEGMENT_MAP_ALPHA = {
    "1": "PB2",
    "2": "PB1",
    "3": "PA",
    "4": "HA",
    "5": "NP",
    "6": "NA",
    "7": "MP",
    "8": "NS",
}

# Influenza B (Betainfluenzavirus)
#   1=PB1, 2=PB2, 3=PA, 4=HA, 5=NP, 6=NA, 7=MP, 8=NS
# Note: PB1 and PB2 are swapped relative to Influenza A.
SEGMENT_MAP_BETA = {
    "1": "PB1",
    "2": "PB2",
    "3": "PA",
    "4": "HA",
    "5": "NP",
    "6": "NA",
    "7": "MP",
    "8": "NS",
}


# Expected full-length segment sizes (in nucleotides) used by length QC.
# Reference:
#   https://pmc.ncbi.nlm.nih.gov/articles/PMC3074182/
SEGMENT_LENGTHS_MAP = {
    "1": 2341,
    "2": 2341,
    "3": 2233,
    "4": 1778,
    "5": 1565,
    "6": 1413,
    "7": 1027,
    "8": 890,
}


# =============================================================================
# Logging
# =============================================================================


class TeeLogger:
    """
    Mirror everything written to stderr into a log file as well.

    Used so that the full run transcript is preserved in `<output>.log`
    while still being visible on the terminal. Writes are buffered
    (1 MiB) for performance on large runs.
    """

    def __init__(self, log_file):
        self.log_file = log_file
        self.terminal = sys.stderr
        self.log_handle = open(
            log_file,
            "a",
            encoding="utf-8",
            buffering=1024 * 1024,
        )

    def write(self, message):
        """Forward a message to both stderr and the log file."""
        self.terminal.write(message)
        self.log_handle.write(message)
        self.terminal.flush()

    def flush(self):
        """Flush both underlying streams."""
        self.terminal.flush()
        self.log_handle.flush()

    def close(self):
        """Flush and close the log file handle."""
        self.log_handle.flush()
        self.log_handle.close()


# =============================================================================
# JSON loading
# =============================================================================


def load_json(file_path):
    """
    Load a JSON file, transparently handling `.gz` compression.

    Parameters
    ----------
    file_path : str
        Path to a `.json` or `.json.gz` file.

    Returns
    -------
    object
        Parsed JSON content (typically a list of records).
    """
    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rt") as handle:
            return json.load(handle)

    with open(file_path, "r") as handle:
        return json.load(handle)


# =============================================================================
# Metadata / sequence lookup
# =============================================================================


def build_metadata_lookup(metadata_records):
    """
    Build a `genome_id -> metadata record` dictionary.

    Duplicate handling
    ------------------
    * If two records share a `genbank_accessions` (or `accession`) value,
      only the first is kept. A warning is emitted for each duplicate.
    * If two records share a `genome_id`, only the first is kept.

    This guarantees a 1-to-1 relationship between genome_id and record
    downstream, which the rest of the pipeline relies on.

    Parameters
    ----------
    metadata_records : list of dict
        Raw metadata records loaded from `meta_raw.json`.

    Returns
    -------
    dict
        Mapping genome_id -> metadata record.
    """
    metadata_lookup = {}
    accession_seen = set()

    for record in metadata_records:
        genome_id = record.get("genome_id")

        if not genome_id:
            print(
                "WARNING: Metadata record missing genome_id; skipping record",
                file=sys.stderr,
            )
            continue

        # ---------------------------------------------------------------------
        # Normalise the accession field into a list of non-empty strings.
        # BV-BRC sometimes provides `genbank_accessions` (list or scalar)
        # and sometimes just `accession`.
        # ---------------------------------------------------------------------
        accession = record.get("genbank_accessions")

        if accession is None:
            accession = record.get("accession")

        if isinstance(accession, list):
            accessions = [
                str(value).strip()
                for value in accession
                if value is not None and str(value).strip()
            ]
        elif accession is not None:
            accession = str(accession).strip()
            accessions = [accession] if accession else []
        else:
            accessions = []

        # ---------------------------------------------------------------------
        # Skip duplicate accession. Keep the first record encountered.
        # ---------------------------------------------------------------------
        duplicate_accession = None

        for acc in accessions:
            if acc in accession_seen:
                duplicate_accession = acc
                break

        if duplicate_accession is not None:
            print(
                "WARNING: Duplicate accession in metadata; "
                f"keeping first record and skipping duplicate: "
                f"{duplicate_accession} "
                f"(genome_id: {genome_id})",
                file=sys.stderr,
            )
            continue

        # ---------------------------------------------------------------------
        # Skip duplicate genome_id. Keep the first record encountered.
        # ---------------------------------------------------------------------
        if genome_id in metadata_lookup:
            print(
                "WARNING: Duplicate genome_id in metadata; "
                f"keeping first record and skipping duplicate: {genome_id}",
                file=sys.stderr,
            )
            continue

        # ---------------------------------------------------------------------
        # Keep the record and register its accessions.
        # ---------------------------------------------------------------------
        metadata_lookup[genome_id] = record

        for acc in accessions:
            accession_seen.add(acc)

    return metadata_lookup


def build_sequence_lookup(sequence_records):
    """
    Build a `genome_id -> sequence` dictionary.

    Duplicate handling
    ------------------
    * Duplicate genome_ids: first record wins.
    * Records missing `genome_id` or `sequence` are skipped with a warning.

    Parameters
    ----------
    sequence_records : list of dict
        Raw sequence records loaded from `sequences_raw.json`.

    Returns
    -------
    dict
        Mapping genome_id -> sequence string.
    """
    sequence_lookup = {}

    for record in sequence_records:
        genome_id = record.get("genome_id")

        if not genome_id:
            print(
                "WARNING: Sequence record missing genome_id; skipping record",
                file=sys.stderr,
            )
            continue

        if genome_id in sequence_lookup:
            print(
                "WARNING: Duplicate genome_id in sequences; "
                f"keeping first record and skipping duplicate: {genome_id}",
                file=sys.stderr,
            )
            continue

        sequence = record.get("sequence")

        if sequence is None:
            print(
                "WARNING: Sequence record missing sequence for "
                f"genome_id {genome_id}; skipping record",
                file=sys.stderr,
            )
            continue

        sequence_lookup[genome_id] = sequence

    return sequence_lookup


# =============================================================================
# Influenza B lineage
# =============================================================================


def extract_lineage_from_other_clinical(metadata):
    """
    Extract a valid Influenza B lineage from the `other_clinical` field.

    BV-BRC stores Influenza B lineage information as free-text entries
    inside `other_clinical`, typically of the form:
        "Lineage: Victoria"
        "Lineage: Yamagata"

    Only the two currently circulating lineages are accepted. Any other
    value (e.g. ancestral B lineages) is ignored.

    Parameters
    ----------
    metadata : dict
        Metadata record.

    Returns
    -------
    str or None
        "Victoria" or "Yamagata" (title-cased), or None if not found.
    """
    other_clinical = metadata.get("other_clinical", [])

    if not isinstance(other_clinical, list):
        other_clinical = [other_clinical]

    for item in other_clinical:
        if not isinstance(item, str):
            continue

        item = item.strip()

        if not item.lower().startswith("lineage:"):
            continue

        lineage = item.split(":", 1)[1].strip()

        if lineage.lower() in {"victoria", "yamagata"}:
            return lineage.title()

    return None


# =============================================================================
# Virus detection
# =============================================================================


def detect_virus_type(metadata):
    """
    Detect whether a record is Influenza A (alpha) or Influenza B (beta).

    Detection order
    ---------------
    1. `genus` field: looks for "alphainfluenzavirus" / "betainfluenzavirus".
    2. `strain` field: falls back to the leading "A/" or "B/" prefix.

    Parameters
    ----------
    metadata : dict
        Metadata record.

    Returns
    -------
    str
        One of {"alpha", "beta", "unknown"}.
    """
    genus = metadata.get("genus", "")

    if genus:
        genus = str(genus).lower()

        if "betainfluenzavirus" in genus:
            return "beta"

        if "alphainfluenzavirus" in genus:
            return "alpha"

    strain = metadata.get("strain", "")

    if strain:
        strain = str(strain).lower()

        if strain.startswith("b/"):
            return "beta"

        if strain.startswith("a/"):
            return "alpha"

    return "unknown"


def get_segment_map(virus_type):
    """
    Return the appropriate segment ID -> name mapping for a virus type.

    Parameters
    ----------
    virus_type : str
        "alpha", "beta", or "unknown".

    Returns
    -------
    dict or None
        The relevant SEGMENT_MAP_* constant, or None if unknown.
    """
    if virus_type == "beta":
        return SEGMENT_MAP_BETA

    if virus_type == "alpha":
        return SEGMENT_MAP_ALPHA

    return None


# =============================================================================
# TYPE_FIELD (subtype for A, lineage for B)
# =============================================================================


def get_type_field(metadata, virus_type):
    """
    Return a valid TYPE_FIELD value for the record.

    Influenza A:  subtype matching `H<num>N<num>(pdm<num>)?`
    Influenza B:  lineage, one of {Victoria, Yamagata}

    Parameters
    ----------
    metadata : dict
        Metadata record.
    virus_type : str
        "alpha", "beta", or "unknown".

    Returns
    -------
    str or None
        Valid subtype/lineage, or None if the value is missing/invalid.
    """
    if not isinstance(metadata, dict):
        return None

    # -------------------------------------------------------------------------
    # Influenza A — validate against the HxNy (optionally pdmN) pattern.
    # -------------------------------------------------------------------------
    if virus_type == "alpha":
        subtype = metadata.get("subtype", "")

        if subtype is None:
            return None

        subtype = str(subtype).strip()

        if re.fullmatch(
            r"H\d+N\d+(?:pdm\d+)?",
            subtype,
            flags=re.IGNORECASE,
        ):
            return subtype

        return None

    # -------------------------------------------------------------------------
    # Influenza B — pull lineage from `other_clinical`.
    # -------------------------------------------------------------------------
    if virus_type == "beta":
        return extract_lineage_from_other_clinical(metadata)

    return None


# =============================================================================
# Metadata field extraction & normalisation
# =============================================================================


def get_accession(metadata):
    """
    Extract and normalise the accession value.

    Accepts either a scalar or a list under `genbank_accessions`
    (falling back to `accession`). Multiple accessions are joined with ",".

    Returns
    -------
    str or None
    """
    accession = metadata.get("genbank_accessions")

    if accession is None:
        accession = metadata.get("accession")

    if accession is None:
        return None

    if isinstance(accession, list):
        accession = [
            str(value).strip()
            for value in accession
            if value is not None and str(value).strip()
        ]

        if not accession:
            return None

        return ",".join(accession)

    accession = str(accession).strip()

    return accession if accession else None


def get_length(metadata):
    """
    Extract `genome_length` as an int.

    Returns
    -------
    int or None
        None if missing or not coercible to int.
    """
    length = metadata.get("genome_length")

    if length is None:
        return None

    try:
        return int(length)
    except (TypeError, ValueError):
        return None


def get_segment_id(metadata):
    """
    Extract a single segment ID from metadata.

    BV-BRC sometimes stores `segment` as a scalar and sometimes as a list.
    A one-element list is accepted; multi-element lists are rejected as
    ambiguous (flagged with the special error code "MULTIPLE_SEGMENTS").

    Returns
    -------
    tuple
        (segment_id, error) where:
            segment_id : str or None
            error      : "MULTIPLE_SEGMENTS" or None
    """
    segment = metadata.get("segment")

    if segment is None:
        return None, None

    if isinstance(segment, list):
        values = [
            str(value).strip()
            for value in segment
            if value is not None and str(value).strip()
        ]

        if not values:
            return None, None

        if len(values) > 1:
            return None, "MULTIPLE_SEGMENTS"

        return values[0], None

    segment = str(segment).strip()

    if not segment:
        return None, None

    return segment, None


def get_host(metadata):
    """
    Extract and standardise the first usable host value.

    Tries, in order: host_common_name, host_name, host_scientific_name.

    Discards obviously non-informative values such as:
        "null", "virus", "u", "unknow(n)", "not available",
        entries containing "unidentified", or age strings like "10 years".

    Normalisation
    -------------
    * Drop trailing "'s ..." possessives (e.g. "Human's" -> "Human").
    * Drop " sp." suffixes (e.g. "Anas sp. mallard" -> "Anas").
    * Replace internal whitespace with underscores and title-case.

    Returns
    -------
    str or None
    """
    for field in (
        "host_common_name",
        "host_name",
        "host_scientific_name",
    ):
        host = metadata.get(field)

        if not host:
            continue

        host = str(host).strip()

        if host.lower() in {
            "null",
            "virus",
            "u",
            "unknow",
            "unknown",
            "not available",
        }:
            continue

        if "unidentified" in host.lower():
            continue

        if re.fullmatch(r"\d+\s+years?", host, re.IGNORECASE):
            continue

        # Strip possessive suffix: "Human's" -> "Human"
        host = re.split(
            r"['’]s\b",
            host,
            maxsplit=1,
            flags=re.IGNORECASE,
        )[0]

        # Strip " sp." species-uncertainty suffix.
        host = re.sub(
            r"\s+sp\..*$",
            "",
            host,
            flags=re.IGNORECASE,
        )

        host = re.sub(r"\s+", "_", host)
        host = host.strip("_")

        if not host:
            continue

        host = host.title()

        return host

    return None


def get_geographic_location(metadata):
    """
    Extract and standardise the geographic location.

    * Truncates anything after the first colon (drops sub-locality detail).
    * Strips a leading "IDCDC-" prefix.
    * Underscores internal whitespace and title-cases the result.

    Returns
    -------
    str or None
    """
    location = metadata.get("geographic_location")

    if not location:
        return None

    location = str(location).strip()

    if ":" in location:
        location = location.split(":", 1)[0].strip()

    location = re.sub(
        r"^IDCDC-",
        "",
        location,
        flags=re.IGNORECASE,
    )

    location = re.sub(r"\s+", "_", location)
    location = location.strip("_")
    location = location.title()

    return location if location else None


def get_strain(metadata):
    """
    Extract and normalise the strain value (whitespace -> underscores).

    Returns
    -------
    str or None
    """
    strain = metadata.get("strain")

    if not strain:
        return None

    strain = str(strain).strip()

    if not strain:
        return None

    strain = re.sub(r"\s+", "_", strain)

    return strain


# =============================================================================
# Required-field validation
# =============================================================================


def validate_required_fields(
    genome_id,
    accession,
    host,
    geographic_location,
    strain,
    segment,
    length,
    type_field,
):
    """
    Return the list of required fields whose value is None.

    Parameters
    ----------
    All positional arguments correspond to the required-field values.
    `genome_id` is still validated here even though it is not written to
    the FASTA header — it is the key used for sequence lookup.

    Returns
    -------
    list of str
        Names of the missing fields (empty list if all present).
    """
    required_fields = {
        "genome_id": genome_id,
        "accession": accession,
        "host": host,
        "geographic_location": geographic_location,
        "strain": strain,
        "segment": segment,
        "length": length,
        "TYPE_FIELD": type_field,
    }

    return [field for field, value in required_fields.items() if value is None]


# =============================================================================
# FASTA headers
# =============================================================================


def create_fasta_header(
    accession,
    host,
    segment_id,
    segment_name,
    type_field,
    geographic_location,
    strain,
    length,
):
    """
    Build the FASTA header for a PASSED sequence.

    Format
    ------
        >accession host|segmentid|segmentName|TYPE_FIELD|
         geographic_location|strain|length

    Notes
    -----
    `genome_id` is intentionally NOT included in the header; the accession
    is the sole record identifier.
    """
    header_parts = [
        f"{accession} {host}".strip(),
        str(segment_id),
        str(segment_name),
        str(type_field),
        str(geographic_location),
        str(strain),
        str(length),
    ]

    return "|".join(header_parts)


def create_excluded_header(
    accession,
    host,
    segment_id,
    segment_name,
    type_field,
    geographic_location,
    strain,
    length,
    reason,
):
    """
    Build the FASTA header for an EXCLUDED sequence.

    Every field is given an explicit `key=value` label, and the exclusion
    reason is appended, so excluded records remain fully self-describing.

    Notes
    -----
    `genome_id` is intentionally NOT included; the accession is the sole
    record identifier.
    """
    fields = [
        f"accession={accession}",
        f"host={host}",
        f"segment_id={segment_id}",
        f"segment_name={segment_name}",
        f"type={type_field}",
        f"geographic_location={geographic_location}",
        f"strain={strain}",
        f"length={length}",
        f"reason={reason}",
    ]

    return "|".join(fields)


# =============================================================================
# Sequence QC
# =============================================================================


def filter_sequence(
    sequence,
    length,
    segment_id,
    min_len_pct,
    max_len_pct,
    max_ambigs,
):
    """
    Apply sequence-level QC filters to a single record.

    Checks performed (in order)
    ---------------------------
    1. EMPTY_SEQUENCE        — sequence is missing or blank.
    2. UNKNOWN_SEGMENT_LENGTH_<id>
                             — segment_id not in SEGMENT_LENGTHS_MAP.
    3. TOO_SHORT_SEGMENT_... — length < min_len_pct * expected length.
    4. TOO_LONG_SEGMENT_...  — length > max_len_pct * expected length
                               (only if max_len_pct is not None).
    5. TOO_MANY_AMBIGS_<n>   — ambiguous bases > max_ambigs
                               (only if max_ambigs >= 0).

    Parameters
    ----------
    sequence : str
        Sequence string (may be empty).
    length : int
        Sequence length as reported in metadata.
    segment_id : str
        Single influenza segment ID ("1".."8").
    min_len_pct : float
        Minimum allowed fraction of the expected segment length.
    max_len_pct : float or None
        Maximum allowed fraction of the expected segment length
        (None = no upper bound).
    max_ambigs : int
        Maximum allowed ambiguous bases (negative = skip the check).

    Returns
    -------
    tuple
        (passes, reason, ambig_count)
            passes      : bool
            reason      : str — "PASSED" or an error code
            ambig_count : int — ambiguous-base count (0 if not computed)
    """
    if sequence is None:
        sequence = ""

    sequence = str(sequence).strip()

    # -------------------------------------------------------------------------
    # 1. Empty sequence
    # -------------------------------------------------------------------------
    if not sequence:
        return False, "EMPTY_SEQUENCE", 0

    # -------------------------------------------------------------------------
    # 2/3/4. Segment-specific length filtering
    # -------------------------------------------------------------------------
    expected_length = SEGMENT_LENGTHS_MAP.get(str(segment_id))

    if expected_length is None:
        return (
            False,
            f"UNKNOWN_SEGMENT_LENGTH_{segment_id}",
            0,
        )

    min_allowed_length = int(expected_length * min_len_pct)

    if length < min_allowed_length:
        return (
            False,
            f"TOO_SHORT_SEGMENT_{segment_id}_{length}_MIN_{min_allowed_length}",
            0,
        )

    if max_len_pct is not None:
        max_allowed_length = int(expected_length * max_len_pct)

        if length > max_allowed_length:
            return (
                False,
                f"TOO_LONG_SEGMENT_{segment_id}_{length}_MAX_{max_allowed_length}",
                0,
            )

    # -------------------------------------------------------------------------
    # 5. Ambiguous bases (IUPAC ambiguity codes)
    # -------------------------------------------------------------------------
    if max_ambigs >= 0:
        ambig_count = sum(1 for base in sequence.upper() if base in "NRYSWKMBDHV")

        if ambig_count > max_ambigs:
            return (
                False,
                f"TOO_MANY_AMBIGS_{ambig_count}",
                ambig_count,
            )

    return True, "PASSED", 0


# =============================================================================
# FASTA formatting
# =============================================================================


def format_fasta_sequence(sequence):
    """
    Upper-case a sequence and wrap it at 60 characters per line,
    following the standard FASTA convention.
    """
    sequence = sequence.upper()

    return "\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60))


# =============================================================================
# Statistics
# =============================================================================


def initialize_stats(metadata_lookup, sequence_lookup):
    """
    Create the statistics dictionary used throughout the run.

    Counters include total sequences, passed/failed counts, per-reason
    failure tallies, per-field missing counts, TYPE_FIELD distribution,
    and a list of kept sequence lengths for summary statistics.
    """
    return {
        "total_sequences": len(metadata_lookup),
        "passed": 0,
        "failed": 0,
        "fail_reasons": {},
        "lengths": [],
        "genomes_processed": 0,
        "type_field_counts": {},
        "missing_field_counts": {},
        "metadata_records": len(metadata_lookup),
        "sequence_records": len(sequence_lookup),
        "missing_sequence_records": 0,
    }


# =============================================================================
# Main conversion
# =============================================================================


def bvbrc2fasta(
    metadata_lookup,
    sequence_lookup,
    output_file,
    excluded_file,
    stats,
    min_len_pct,
    max_len_pct,
    max_ambigs,
    verbose=True,
):
    """
    Convert BV-BRC JSON lookups into passed/excluded FASTA files.

    Iterates over the metadata lookup (metadata drives processing) and,
    for each genome_id:

        1. Looks up the sequence.
        2. Detects virus type and selects the right segment map.
        3. Extracts TYPE_FIELD and other required metadata.
        4. Validates required fields.
        5. Applies sequence QC.
        6. Writes to the passed or excluded FASTA accordingly.

    Records may fail for multiple reasons simultaneously; all reasons are
    concatenated with ";" into a single `reason` string.

    Notes
    -----
    `genome_id` is used internally as the lookup key, but it is NOT written
    to either FASTA header. The accession is the sole record identifier.

    Parameters
    ----------
    metadata_lookup : dict
        genome_id -> metadata record.
    sequence_lookup : dict
        genome_id -> sequence string.
    output_file : str
        Path to the passed FASTA.
    excluded_file : str
        Path to the excluded FASTA.
    stats : dict
        Statistics dict from `initialize_stats` (mutated in place).
    min_len_pct, max_len_pct, max_ambigs :
        QC thresholds (see `filter_sequence`).
    verbose : bool
        If True, print a per-genome progress line to stderr.

    Returns
    -------
    dict
        The updated `stats` dict.
    """
    with (
        open(output_file, "w") as f_out,
        open(excluded_file, "w") as f_excl,
    ):
        # Metadata drives processing.
        for genome_id, metadata in metadata_lookup.items():
            stats["genomes_processed"] += 1

            # -----------------------------------------------------------------
            # Sequence is looked up by genome_id; may be absent.
            # -----------------------------------------------------------------
            seq = sequence_lookup.get(genome_id)

            # -----------------------------------------------------------------
            # Determine virus type and corresponding segment map.
            # -----------------------------------------------------------------
            virus_type = detect_virus_type(metadata)
            segment_map = get_segment_map(virus_type)

            # -----------------------------------------------------------------
            # Extract and normalise all metadata fields once.
            # -----------------------------------------------------------------
            type_field = get_type_field(metadata, virus_type)

            if type_field:
                stats["type_field_counts"][type_field] = (
                    stats["type_field_counts"].get(type_field, 0) + 1
                )

            accession = get_accession(metadata)
            segment_id, segment_error = get_segment_id(metadata)
            length = get_length(metadata)
            host = get_host(metadata)
            geographic_location = get_geographic_location(metadata)
            strain = get_strain(metadata)

            segment_name = None

            if segment_id is not None and segment_map:
                segment_name = segment_map.get(segment_id)

            if verbose:
                type_label = virus_type.upper()
                print(
                    f"   Processing {genome_id} {accession} {type_label}",
                    file=sys.stderr,
                )

            # -----------------------------------------------------------------
            # Metadata validation.
            # -----------------------------------------------------------------
            missing_fields = validate_required_fields(
                genome_id,
                accession,
                host,
                geographic_location,
                strain,
                segment_id,
                length,
                type_field,
            )

            reasons = []

            # Multiple segment values are invalid, but not "missing".
            if segment_error == "MULTIPLE_SEGMENTS":
                if "segment" in missing_fields:
                    missing_fields.remove("segment")

                reasons.append("INVALID_SEGMENT_MULTIPLE_VALUES")

            if missing_fields:
                reason = "MISSING_REQUIRED_FIELDS_" + "_".join(missing_fields)
                reasons.append(reason)

                for field in missing_fields:
                    stats["missing_field_counts"][field] = (
                        stats["missing_field_counts"].get(field, 0) + 1
                    )

            # -----------------------------------------------------------------
            # Sequence availability and QC.
            # -----------------------------------------------------------------
            if genome_id not in sequence_lookup:
                stats["missing_sequence_records"] += 1
                reasons.append("MISSING_SEQUENCE_RECORD")
            else:
                passes_qc, qc_reason, ambig_count = filter_sequence(
                    sequence=seq,
                    length=length,
                    segment_id=segment_id,
                    min_len_pct=min_len_pct,
                    max_len_pct=max_len_pct,
                    max_ambigs=max_ambigs,
                )

                if not passes_qc:
                    reasons.append(qc_reason)

            # -----------------------------------------------------------------
            # Final pass/fail decision.
            # -----------------------------------------------------------------
            if reasons:
                passes = False
                reason = ";".join(reasons)
            else:
                passes = True
                reason = "PASSED"

            # -----------------------------------------------------------------
            # Write to the appropriate FASTA.
            # -----------------------------------------------------------------
            if passes:
                header = create_fasta_header(
                    accession,
                    host,
                    segment_id,
                    segment_name,
                    type_field,
                    geographic_location,
                    strain,
                    length,
                )

                f_out.write(f">{header}\n")
                f_out.write(format_fasta_sequence(seq) + "\n")

                stats["passed"] += 1
                stats["lengths"].append(length)
            else:
                header = create_excluded_header(
                    accession,
                    host,
                    segment_id,
                    segment_name,
                    type_field,
                    geographic_location,
                    strain,
                    length,
                    reason,
                )

                f_excl.write(f">{header}\n")

                # Only emit a sequence line if there was actually a sequence.
                if seq:
                    f_excl.write(format_fasta_sequence(seq) + "\n")

                stats["failed"] += 1

                stats["fail_reasons"][reason] = (
                    stats["fail_reasons"].get(reason, 0) + 1
                )

    return stats


# =============================================================================
# Summary
# =============================================================================


def print_summary(stats, output_file, excluded_file):
    """
    Print the end-of-run summary to stderr (and therefore to the log).

    Includes:
        * Missing sequence-record count
        * Missing required-field counts
        * Failure-reason breakdown
        * TYPE_FIELD distribution
        * Passed / excluded counts and percentages
        * Length statistics for kept sequences
        * Input record counts and output file paths
    """
    print("\n" + "=" * 80, file=sys.stderr)
    print("📊 SUMMARY", file=sys.stderr)
    print("=" * 80, file=sys.stderr)

    total = stats["passed"] + stats["failed"]

    # -------------------------------------------------------------------------
    # Missing sequence records
    # -------------------------------------------------------------------------
    if stats["missing_sequence_records"] > 0:
        print(
            "\n   Missing sequence records:",
            file=sys.stderr,
        )
        print(
            f"     {stats['missing_sequence_records']:,}",
            file=sys.stderr,
        )

    # -------------------------------------------------------------------------
    # Missing required fields
    # -------------------------------------------------------------------------
    if stats["missing_field_counts"]:
        print(
            "\n   Missing required fields:",
            file=sys.stderr,
        )

        for field, count in sorted(
            stats["missing_field_counts"].items(),
            key=lambda item: item[1],
            reverse=True,
        ):
            print(
                f"     {field}: {count:,}",
                file=sys.stderr,
            )

    # -------------------------------------------------------------------------
    # Failure reasons
    # -------------------------------------------------------------------------
    if stats["failed"] > 0:
        print(
            "\n   Failure reasons:",
            file=sys.stderr,
        )

        for reason, count in sorted(
            stats["fail_reasons"].items(),
            key=lambda item: item[1],
            reverse=True,
        ):
            print(
                f"     {reason}: {count:,}",
                file=sys.stderr,
            )

    # -------------------------------------------------------------------------
    # TYPE_FIELD distribution
    # -------------------------------------------------------------------------
    if stats["type_field_counts"]:
        print(
            "\n   Type field distribution:",
            file=sys.stderr,
        )

        for field, count in sorted(
            stats["type_field_counts"].items(),
            key=lambda item: item[1],
            reverse=True,
        ):
            print(
                f"     {field}: {count:,}",
                file=sys.stderr,
            )

    # -------------------------------------------------------------------------
    # Processing summary
    # -------------------------------------------------------------------------
    print(
        f"\n   Sequences processed: {stats['genomes_processed']:,}",
        file=sys.stderr,
    )
    print(
        f"   Sequences kept:      {stats['passed']:,}",
        file=sys.stderr,
    )
    print(
        f"   Sequences excluded:  {stats['failed']:,}",
        file=sys.stderr,
    )

    if total > 0:
        passed_pct = stats["passed"] / total * 100
        failed_pct = stats["failed"] / total * 100
    else:
        passed_pct = 0
        failed_pct = 0

    print(
        f"   Kept:                {passed_pct:.1f}%",
        file=sys.stderr,
    )
    print(
        f"   Excluded:            {failed_pct:.1f}%",
        file=sys.stderr,
    )

    # -------------------------------------------------------------------------
    # Length statistics for kept sequences
    # -------------------------------------------------------------------------
    if stats["lengths"]:
        lengths = stats["lengths"]

        print(
            "\n   Length statistics (kept):",
            file=sys.stderr,
        )
        print(
            f"     Min:   {min(lengths):,} bp",
            file=sys.stderr,
        )
        print(
            f"     Max:   {max(lengths):,} bp",
            file=sys.stderr,
        )
        print(
            f"     Avg:   {sum(lengths) / len(lengths):.1f} bp",
            file=sys.stderr,
        )
        print(
            f"     Total: {sum(lengths):,} bp",
            file=sys.stderr,
        )

    # -------------------------------------------------------------------------
    # Record counts
    # -------------------------------------------------------------------------
    print(
        f"\n   Metadata records: {stats['metadata_records']:,}",
        file=sys.stderr,
    )
    print(
        f"   Sequence records: {stats['sequence_records']:,}",
        file=sys.stderr,
    )

    # -------------------------------------------------------------------------
    # Output files
    # -------------------------------------------------------------------------
    print(
        f"   Output FASTA: {output_file}",
        file=sys.stderr,
    )

    if stats["failed"] > 0:
        print(
            f"   Excluded FASTA: {excluded_file}",
            file=sys.stderr,
        )


# =============================================================================
# Command-line interface
# =============================================================================


def main():
    """
    Parse arguments, load data, build lookups, run conversion, print summary.

    The full run is logged to `<output>.log` (or `--log`) via TeeLogger,
    which mirrors stderr to both the terminal and the log file.

    Only a single total-runtime measurement is taken for the whole script.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Convert BV-BRC raw JSON to FASTA with QC. "
            "Handles both Alpha and Beta Influenza."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:

  # Basic conversion
  bvbrc2fasta.py -m meta_raw.json -s sequences_raw.json -o output.fasta

  # QC filtering
  bvbrc2fasta.py -m meta_raw.json -s sequences_raw.json -o output.fasta \\
      --min-len-pct 0.9 --max-len-pct 1.1 --max-ambigs 5

HEADER FORMAT:

  Passed:
  >accession host|segmentid|segmentName|TYPE_FIELD|geographic_location|strain|length

  Excluded:
  >accession=...|host=...|segment_id=...|segment_name=...|type=...|geographic_location=...|strain=...|length=...|reason=...

  Note: genome_id is used internally for lookups but is NOT written to FASTA headers.

REQUIRED FIELDS:

  accession | genome_id | host | geographic_location |
  strain | segment | length | TYPE_FIELD

TYPE_FIELD:

  Influenza A = SUBTYPE
  Influenza B = LINEAGE

SEGMENTS:

  Influenza A:
    1=PB2  2=PB1  3=PA  4=HA
    5=NP   6=NA   7=MP  8=NS

  Influenza B:
    1=PB1  2=PB2  3=PA  4=HA
    5=NP   6=NA   7=MP  8=NS
        """,
    )

    # =========================================================================
    # Input / output
    # =========================================================================
    parser.add_argument(
        "-m",
        "--metadata",
        required=True,
        help="Input metadata JSON file (meta_raw.json)",
    )
    parser.add_argument(
        "-s",
        "--sequences",
        required=True,
        help="Input sequences JSON file (sequences_raw.json)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output FASTA file (passed sequences)",
    )
    parser.add_argument(
        "--excluded",
        default=None,
        help="Output FASTA file for excluded sequences "
             "(default: <output>_excluded.fasta)",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Log file for stderr output (default: <output>.log)",
    )

    # =========================================================================
    # QC
    # =========================================================================
    parser.add_argument(
        "--min-len-pct",
        type=float,
        default=0.9,
        help=(
            "Minimum length range multiplier for expected segment length "
            "(default: 0.9)"
        ),
    )
    parser.add_argument(
        "--max-len-pct",
        type=float,
        default=None,
        help=(
            "Maximum length range multiplier for expected segment length "
            "(default: no limit)"
        ),
    )
    parser.add_argument(
        "--max-ambigs",
        type=int,
        default=0,
        help="Maximum allowed ambiguous bases (default: 0)",
    )

    # =========================================================================
    # Miscellaneous
    # =========================================================================
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress per-genome verbose output",
    )

    args = parser.parse_args()

    # =========================================================================
    # Validate QC arguments
    # =========================================================================
    if args.min_len_pct <= 0:
        parser.error("--min-len-pct must be > 0")

    if args.max_len_pct is not None and args.max_len_pct <= 0:
        parser.error("--max-len-pct must be greater than 0")

    if args.max_len_pct is not None and args.min_len_pct > args.max_len_pct:
        parser.error("--min-len-pct cannot be greater than --max-len-pct")

    # =========================================================================
    # Validate input files
    # =========================================================================
    for file_path in (args.metadata, args.sequences):
        if not os.path.exists(file_path):
            print(
                f"❌ Error: File not found: {file_path}",
                file=sys.stderr,
            )
            sys.exit(1)

    # =========================================================================
    # Determine excluded output path
    # =========================================================================
    if args.excluded is None:
        stem, suffix = os.path.splitext(args.output)
        excluded_file = f"{stem}_excluded{suffix}"
    else:
        excluded_file = args.excluded

    # =========================================================================
    # Determine log output path
    # =========================================================================
    if args.log is None:
        stem, _ = os.path.splitext(args.output)
        log_file = f"{stem}.log"
    else:
        log_file = args.log

    # =========================================================================
    # Start total runtime timer (whole script)
    # =========================================================================
    script_start = time.perf_counter()

    original_stderr = sys.stderr
    logger = None

    try:
        # =====================================================================
        # Start logger (mirrors stderr to the log file)
        # =====================================================================
        logger = TeeLogger(log_file)
        sys.stderr = logger

        print("=" * 80, file=sys.stderr)
        print("🧬 BV-BRC TO FASTA CONVERTER", file=sys.stderr)
        print("=" * 80, file=sys.stderr)

        print("\n📁 Input files:", file=sys.stderr)
        print(f"   Metadata: {args.metadata}", file=sys.stderr)
        print(f"   Sequences: {args.sequences}", file=sys.stderr)
        print(f"   Output: {args.output}", file=sys.stderr)
        print(f"   Excluded: {excluded_file}", file=sys.stderr)
        print(f"   Log: {log_file}", file=sys.stderr)

        # =====================================================================
        # Load JSON
        # =====================================================================
        print("\n🔍 LOADING DATA", file=sys.stderr)
        print("-" * 60, file=sys.stderr)

        metadata_records = load_json(args.metadata)
        sequence_records = load_json(args.sequences)

        print(
            f"   Metadata records: {len(metadata_records):,}",
            file=sys.stderr,
        )
        print(
            f"   Sequence records: {len(sequence_records):,}",
            file=sys.stderr,
        )

        # =====================================================================
        # Build lookups
        # =====================================================================
        print("\n🔗 BUILDING LOOKUPS", file=sys.stderr)
        print("-" * 60, file=sys.stderr)

        metadata_lookup = build_metadata_lookup(metadata_records)
        sequence_lookup = build_sequence_lookup(sequence_records)

        print(
            f"   Metadata lookup: {len(metadata_lookup):,} genomes",
            file=sys.stderr,
        )
        print(
            f"   Sequence lookup: {len(sequence_lookup):,} genomes",
            file=sys.stderr,
        )

        # =====================================================================
        # Initialize statistics
        # =====================================================================
        stats = initialize_stats(metadata_lookup, sequence_lookup)

        # =====================================================================
        # Report QC settings
        # =====================================================================
        print("\n🔬 QC FILTERS", file=sys.stderr)
        print("-" * 60, file=sys.stderr)

        print(
            f"   Minimum length range: {args.min_len_pct}",
            file=sys.stderr,
        )
        print(
            f"   Maximum length range: {args.max_len_pct}",
            file=sys.stderr,
        )
        print(
            "   Expected segment lengths: "
            + ", ".join(
                f"{segment}={expected:,} bp"
                for segment, expected in SEGMENT_LENGTHS_MAP.items()
            ),
            file=sys.stderr,
        )
        print(
            f"   Max ambiguous bases: {args.max_ambigs}",
            file=sys.stderr,
        )
        print(
            "   Required fields: "
            "accession, genome_id, host, "
            "geographic_location, strain, "
            "segment, length, TYPE_FIELD",
            file=sys.stderr,
        )
        print("-" * 60, file=sys.stderr)

        # =====================================================================
        # Conversion
        # =====================================================================
        stats = bvbrc2fasta(
            metadata_lookup,
            sequence_lookup,
            output_file=args.output,
            excluded_file=excluded_file,
            stats=stats,
            min_len_pct=args.min_len_pct,
            max_len_pct=args.max_len_pct,
            max_ambigs=args.max_ambigs,
            verbose=not args.quiet,
        )

        # =====================================================================
        # Summary
        # =====================================================================
        print_summary(stats, args.output, excluded_file)

        # =====================================================================
        # Runtime (total only)
        # =====================================================================
        total_time = time.perf_counter() - script_start

        print("\n" + "=" * 80, file=sys.stderr)
        print("⏱️ RUNTIME", file=sys.stderr)
        print("=" * 80, file=sys.stderr)

        print(
            f"   Total script runtime:    {total_time:.2f} seconds",
            file=sys.stderr,
        )

        if total_time >= 60:
            print(
                f"   Total script runtime:    {total_time / 60:.2f} minutes",
                file=sys.stderr,
            )

        if total_time >= 3600:
            print(
                f"   Total script runtime:    {total_time / 3600:.2f} hours",
                file=sys.stderr,
            )

        if total_time > 0:
            print(
                "   Processing rate:         "
                f"{stats['genomes_processed'] / total_time:,.2f} "
                "genomes/sec",
                file=sys.stderr,
            )

        print("\n" + "=" * 80, file=sys.stderr)
        print("✅ CONVERSION COMPLETE", file=sys.stderr)
        print("=" * 80, file=sys.stderr)

    except Exception:
        print(
            "\n❌ ERROR: Conversion failed.",
            file=sys.stderr,
        )
        raise

    finally:
        # Calculate final runtime before closing the logger.
        total_time = time.perf_counter() - script_start

        if logger is not None:
            print(
                f"\nTotal script runtime: {total_time:.2f} seconds",
                file=sys.stderr,
            )
            logger.flush()

        # Restore normal stderr.
        sys.stderr = original_stderr

        # Close log file.
        if logger is not None:
            logger.close()


# =============================================================================
# Entry point
# =============================================================================


if __name__ == "__main__":
    main()
