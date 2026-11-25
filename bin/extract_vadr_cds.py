#!/usr/bin/env python3

# betainfluenza virus structure
# https://viralzone.expasy.org/80

# alphainfluenza virus structure
# https://viralzone.expasy.org/6


"""
Extract full influenza A segment/gene sequences from VADR output.

The script uses:

    1. FASTA sequence name
    2. VADR .sqc file, matched by sequence name
    3. VADR .ftr file, matched by sequence name

Reference extraction:

    - Only PASS gene features are used.
    - All gene features for a sequence are considered.
    - The reference region is the enclosing genomic interval:
          minimum gene start .. maximum gene end
    - If all genes are on the '+' strand:
          extract the interval as-is.
    - If all genes are on the '-' strand:
          extract the interval and reverse-complement it.
    - If genes occur on both '+' and '-' strands:
          report an ERROR and skip the sequence.

Sequence name:

    The sequence name is the first whitespace-delimited field
    in the FASTA header.

Example:

    >PV717827|11320.1516956 Wild_Waterbird|2|PB1|H2N3|Australia|...

Sequence name:

    PV717827|11320.1516956

The same sequence name is expected in the SQC and FTR files.

Original segment:

    The original segment is obtained from the original FASTA
    description, before any VADR debug information is appended.

Example:

    >OL371557|11320.304589 Wild_Waterbird|2|PB1|H2N3|Australia|...

    description fields:

        Wild_Waterbird | 2 | PB1 | H2N3 | ...

    Therefore original segment = 2.

The output header contains additional VADR information for debugging.
"""

import argparse
import logging
import re
import sys
from pathlib import Path
from logging.handlers import RotatingFileHandler
from typing import Optional, Union

# https://www.nature.com/articles/nature04239/tables/1
# https://pmc.ncbi.nlm.nih.gov/articles/PMC3074182/
# https://www.nature.com/articles/nature04239/tables/1
SEGMENT_LENGTHS_MAP = {
    "1": 2274,
    "2": 2280,
    "3": 2151,
    "4": 1701,
    "5": 1497,
    "6": 1410,
    "7": 982,
    "8": 838,
}


# Create a module-level logger for main.py
# This satisfies your linter and captures the file context correctly.
logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------
# Create a module-level logger for main.py
# This satisfies your linter and captures the file context correctly.
logger = logging.getLogger(__name__)


def setup_logging(log_file=None, error_log_file=None):
    """Configure logging to console, main log, and error log.

    Console:      INFO and above
    Main log:     INFO and above
    Error log:    ERROR and above
    """
    # 1. Access root logger and set global threshold to lowest required level
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # 2. Clear existing handlers cleanly
    if logger.hasHandlers():
        logger.handlers.clear()

    # 3. Use an easy-to-read timestamp format
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # 4. Console handler (Using stdout instead of stderr for standard info)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Helper function to prevent repeating file handler boilerplate
    def _add_file_handler(file_path: Union[str, Path], level: int) -> None:
        path = Path(file_path)
        path.parent.mkdir(parents=True, exist_ok=True)

        # 5. Use RotatingFileHandler to protect disk space (max 10MB per file, keeps 5 backups)
        file_handler = RotatingFileHandler(
            path, maxBytes=10 * 1024 * 1024, backupCount=5, encoding="utf-8"
        )
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # 6. Main log
    if log_file:
        _add_file_handler(log_file, logging.INFO)

    # 7. Error log
    if error_log_file:
        _add_file_handler(error_log_file, logging.ERROR)


# ----------------------------------------------------------------------
# FASTA
# ----------------------------------------------------------------------


def read_fasta(fasta_file):
    """
    Read FASTA file.

    Returns:

        {
            seq_name: {
                "header": original_header,
                "sequence": sequence,
            }
        }

    Sequence name is the first whitespace-delimited field.
    """

    sequences = {}

    header = None
    sequence_parts = []

    with open(fasta_file) as fh:
        for line_number, line in enumerate(fh, start=1):
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                # Save previous record.
                if header is not None:
                    _store_fasta_record(
                        sequences,
                        header,
                        sequence_parts,
                        fasta_file,
                    )

                header = line[1:].strip()
                sequence_parts = []

            else:
                if header is None:
                    logger.error(
                        "FASTA sequence encountered before first header "
                        "in %s at line %d",
                        fasta_file,
                        line_number,
                    )
                    continue

                sequence_parts.append(line.strip())

    # Save final record.
    if header is not None:
        _store_fasta_record(
            sequences,
            header,
            sequence_parts,
            fasta_file,
        )

    logger.info(
        "Read %d FASTA sequences from %s",
        len(sequences),
        fasta_file,
    )

    return sequences


def _store_fasta_record(
    sequences,
    header,
    sequence_parts,
    fasta_file,
):
    """
    Store one FASTA record and detect duplicate sequence names.
    """

    sequence_name = header.split(None, 1)[0]
    sequence = "".join(sequence_parts).upper()

    if sequence_name in sequences:
        logger.error(
            "Duplicate FASTA sequence name: %s (file: %s). Record will be skipped.",
            sequence_name,
            fasta_file,
        )
        return

    sequences[sequence_name] = {
        "header": header,
        "sequence": sequence,
    }


# ----------------------------------------------------------------------
# SQC
# ----------------------------------------------------------------------


def parse_sqc(sqc_file):
    """
    Parse VADR .sqc file.

    Returns:

        {
            seq_name: {
                "status": ...,
                "model": ...,
                "segment": ...,
            }
        }

    Sequence name is fields[1].
    """

    annotations = {}

    with open(sqc_file) as fh:
        for line_number, line in enumerate(fh, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 7:
                logger.error(
                    "Malformed SQC line %d in %s: %s",
                    line_number,
                    sqc_file,
                    line,
                )
                continue

            try:
                seq_name = fields[1]
                status = fields[3]
                model = fields[5]
                segment = fields[6]

            except IndexError:
                logger.error(
                    "Unable to parse SQC line %d in %s: %s",
                    line_number,
                    sqc_file,
                    line,
                )
                continue

            if seq_name in annotations:
                logger.error(
                    "Duplicate SQC sequence name: %s (line %d in %s).",
                    seq_name,
                    line_number,
                    sqc_file,
                )
                continue

            annotations[seq_name] = {
                "status": status,
                "model": model,
                "segment": segment,
            }

    logger.info(
        "Read %d SQC annotations from %s",
        len(annotations),
        sqc_file,
    )

    return annotations


# ----------------------------------------------------------------------
# FTR coordinate parsing
# ----------------------------------------------------------------------


def parse_coords(coords):
    # Pre-compile the regex pattern for rapid lookups and execution speed
    COORDS_RE = re.compile(r"(\d+)\.\.(\d+):([+-])")

    """Parse VADR sequence coordinates. e.g., '26..51:+,740..1007:+'"""
    if not coords or coords in (".", "-"):
        return []

    # Find all valid match segments in the string
    matches = COORDS_RE.findall(coords)

    # If the number of reconstructed parts doesn't match the original structure, it's invalid
    if not matches or len(coords.split(",")) != len(matches):
        raise ValueError(f"Invalid VADR coordinates: {coords}")

    parsed = []
    for start_str, end_str, strand in matches:
        start, end = int(start_str), int(end_str)

        # Enforce start <= end orientation
        if start > end:
            start, end = end, start

        parsed.append((start, end, strand))

    return parsed


# ----------------------------------------------------------------------
# FTR
# ----------------------------------------------------------------------


def parse_ftr(ftr_file):
    """Parse VADR .ftr file. Only PASS gene features are retained."""
    annotations = {}

    with open(ftr_file) as fh:
        for line_number, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 24:
                logger.error(
                    "Malformed FTR line %d in %s: expected at least 24 fields, got %d: %s",
                    line_number,
                    ftr_file,
                    len(fields),
                    line,
                )
                continue

            # Descriptive tuple unpacking handles index assignment cleanly
            ftr_index, seq_name = fields[0], fields[1]
            status, feature_type, feature_name = fields[3], fields[5], fields[6]
            seq_coords = fields[23]

            # Early filtration guard clauses
            if status != "PASS" or feature_type.lower() != "gene":
                continue

            try:
                parsed_coords = parse_coords(seq_coords)
                if not parsed_coords:
                    raise ValueError("Coordinates empty")
            except ValueError as exc:
                logger.error(
                    "Invalid/missing coordinates for gene %s (sequence %s, FTR %s, line %d): %s",
                    feature_name,
                    seq_name,
                    ftr_index,
                    line_number,
                    exc,
                )
                continue

            # Ensure uniform strandedness (+ or -)
            strands = {strand for _, _, strand in parsed_coords}
            if len(strands) != 1:
                logger.error(
                    "Gene feature has mixed strands within coordinates: sequence=%s gene=%s coords=%s FTR=%s",
                    seq_name,
                    feature_name,
                    seq_coords,
                    ftr_index,
                )
                continue
            strand = strands.pop()

            # Pythonic math tracking
            start = min(coord[0] for coord in parsed_coords)
            end = max(coord[1] for coord in parsed_coords)
            coord_length = sum(e - s + 1 for s, e, _ in parsed_coords)

            # Dictionary population using simple dictionary layout
            annotations.setdefault(seq_name, []).append(
                {
                    "feature": feature_name,
                    "coordinates": seq_coords,
                    "strand": strand,
                    "start": start,
                    "end": end,
                    "length": coord_length,
                    "status": status,
                    "ftr_index": ftr_index,
                }
            )

    logger.info(
        "Read PASS gene annotations for %d sequences from %s",
        len(annotations),
        ftr_file,
    )
    return annotations


def get_original_segment(header):
    """
    Get the original segment from the original FASTA header.

    Expected original header format:

        >SEQNAME HOST|SEGMENT|GENE|SUBTYPE|LOCATION|STRAIN|LENGTH

    Example:

        >OL371557|11320.304589
        Wild_Waterbird|2|PB1|H2N3|Australia|...|2274

    Therefore segment = 2.

    This function intentionally does NOT search for:
        fluA-seg1
        vadr|...

    because those values may have been added later by this script.
    """

    parts = header.split(None, 1)

    if len(parts) < 2:
        return None

    description = parts[1]

    fields = description.split("|")

    if len(fields) < 2:
        return None

    segment_value = fields[1].strip()

    if re.fullmatch(r"[1-8]", segment_value):
        return int(segment_value)

    return None


# ----------------------------------------------------------------------
# Gene/reference region
# ----------------------------------------------------------------------


def determine_reference_region(seq_name, genes):
    """
    Determine the full genomic reference region from all PASS genes.

    Rules:

        all '+':
            min(start) .. max(end), '+'

        all '-':
            min(start) .. max(end), '-'

        mixed '+/-':
            ERROR and return None
    """

    if not genes:
        logger.error(
            "No PASS gene annotations for sequence %s.",
            seq_name,
        )
        return None

    strands = {gene["strand"] for gene in genes}

    if len(strands) > 1:
        logger.error(
            "Mixed gene strands for sequence %s. Reference extraction skipped.",
            seq_name,
        )

        for gene in genes:
            logger.error(
                "  gene=%s coords=%s strand=%s FTR=%s",
                gene["feature"],
                gene["coordinates"],
                gene["strand"],
                gene["ftr_index"],
            )

        return None

    strand = next(iter(strands))

    start = min(gene["start"] for gene in genes)

    end = max(gene["end"] for gene in genes)

    gene_names = [gene["feature"] for gene in genes]

    return {
        "start": start,
        "end": end,
        "strand": strand,
        "coordinates": f"{start}..{end}:{strand}",
        "gene_names": gene_names,
    }


# ----------------------------------------------------------------------
# Sequence extraction
# ----------------------------------------------------------------------


def reverse_complement(sequence):
    """
    Return reverse complement of a nucleotide sequence.
    """

    complement = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan",
    )

    return sequence.translate(complement)[::-1]


def extract_reference(
    sequence,
    start,
    end,
    strand,
):
    """
    Extract one continuous genomic interval.

    VADR coordinates are 1-based inclusive.

    Therefore:

        start=1, end=10

    becomes:

        sequence[0:10]

    For '-' strand, the entire extracted interval is
    reverse-complemented.

    IMPORTANT:
        We do NOT concatenate split coordinate pieces here.

        For example:

            26..51:+,740..1007:+

        produces:

            26..1007

        rather than:

            26..51 + 740..1007
    """

    if start < 1:
        raise ValueError(f"Invalid start coordinate: {start}")

    if end > len(sequence):
        raise ValueError(
            f"End coordinate {end} exceeds sequence length {len(sequence)}"
        )

    if start > end:
        raise ValueError(
            f"Start coordinate {start} is greater than end coordinate {end}"
        )

    subsequence = sequence[start - 1 : end]

    if strand == "-":
        subsequence = reverse_complement(subsequence)

    return subsequence


# ----------------------------------------------------------------------
# QC
# ----------------------------------------------------------------------


def count_ambiguous(sequence):
    """
    Count bases other than A/T/G/C.
    """

    return sum(base not in "ATGC" for base in sequence.upper())


# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------


def write_fasta_record(
    output_handle,
    header,
    sequence,
):
    """
    Write one FASTA record.
    """

    output_handle.write(f">{header}|len={len(sequence)}\n")

    for i in range(0, len(sequence), 80):
        output_handle.write(sequence[i : i + 80] + "\n")


def extract_sequences(
    fasta_file,
    sqc_file,
    ftr_file,
    output_file,
    min_len_pct=0.90,
    max_len_pct=1.10,
    maxambigs=None,
):
    """
    Extract reference sequences.

    FASTA, SQC, and FTR are joined exclusively by sequence name.
    """
    fasta_sequences = read_fasta(fasta_file)
    sqc_annotations = parse_sqc(sqc_file)
    ftr_annotations = parse_ftr(ftr_file)

    total = len(fasta_sequences)
    extracted = 0
    excluded = 0

    logger.info("Starting extraction of %d FASTA sequences.", total)

    with open(output_file, "w") as output_handle:
        for seq_name, fasta_record in fasta_sequences.items():
            header = fasta_record["header"]
            sequence = fasta_record["sequence"]

            logger.info("Processing sequence: %s", seq_name)

            # ----------------------------------------------------------
            # SQC lookup
            # ----------------------------------------------------------
            sqc = sqc_annotations.get(seq_name)
            if sqc is None:
                logger.error(
                    "No SQC annotation found for sequence %s. Header: %s",
                    seq_name,
                    header,
                )
                excluded += 1
                continue

            if sqc["status"] != "PASS":
                logger.error(
                    "SQC status is not PASS for sequence %s: status=%s model=%s segment=%s. Header: %s",
                    seq_name,
                    sqc["status"],
                    sqc["model"],
                    sqc["segment"],
                    header,
                )
                excluded += 1
                continue

            # ----------------------------------------------------------
            # FTR lookup
            # ----------------------------------------------------------
            genes = ftr_annotations.get(seq_name)
            if not genes:
                logger.error(
                    "No PASS gene annotations found for sequence %s. Header: %s",
                    seq_name,
                    header,
                )
                excluded += 1
                continue

            # ----------------------------------------------------------
            # Determine enclosing reference region
            # ----------------------------------------------------------
            region = determine_reference_region(seq_name, genes)
            if region is None:
                excluded += 1
                continue

            start = region["start"]
            end = region["end"]
            strand = region["strand"]

            # ----------------------------------------------------------
            # Extract sequence
            # ----------------------------------------------------------
            try:
                extracted_sequence = extract_reference(
                    sequence=sequence,
                    start=start,
                    end=end,
                    strand=strand,
                )
            except ValueError as exc:
                logger.error(
                    "Failed to extract reference region for sequence %s: %s. Header: %s",
                    seq_name,
                    exc,
                    header,
                )
                excluded += 1
                continue

            extracted_length = len(extracted_sequence)

            # ----------------------------------------------------------
            # Segment-Specific Length QC
            # ----------------------------------------------------------
            seg_key = str(get_original_segment(header))

            if seg_key not in SEGMENT_LENGTHS_MAP:
                logger.error(
                    "Segment %s found for sequence %s is missing from standard lengths configuration map. Header: %s",
                    seg_key,
                    seq_name,
                    header,
                )
                excluded += 1
                continue

            expected_len = SEGMENT_LENGTHS_MAP[seg_key]
            min_allowed = expected_len * min_len_pct
            #max_allowed = expected_len * max_len_pct

            if extracted_length < min_allowed:
                logger.error(
                    "Sequence %s excluded: extracted length %d is below minimum threshold of %.1f (expected: %d). Region=%s Header=%s",
                    seq_name,
                    extracted_length,
                    min_allowed,
                    expected_len,
                    region["coordinates"],
                    header,
                )
                excluded += 1
                continue

            if max_len_pct is not None:
                max_allowed = int(expected_len * max_len_pct)

                if extracted_length > max_allowed:
                    logger.error(
                        "Sequence %s excluded: extracted length %d is "
                        "above maximum threshold of %d "
                        "(expected=%d).",
                        seq_name,
                        extracted_length,
                        max_allowed,
                        expected_len,
                    )
                    excluded += 1
                    continue

            # ----------------------------------------------------------
            # Ambiguous-base QC
            # ----------------------------------------------------------
            ambiguous_count = count_ambiguous(extracted_sequence)
            if maxambigs is not None and ambiguous_count > maxambigs:
                logger.error(
                    "Sequence %s excluded: %d ambiguous bases (maximum %d). Region=%s Header=%s",
                    seq_name,
                    ambiguous_count,
                    maxambigs,
                    region["coordinates"],
                    header,
                )
                excluded += 1
                continue

            # ----------------------------------------------------------
            # Output header
            # ----------------------------------------------------------
            gene_names = ",".join(region["gene_names"])

            output_header = (
                f"{header}"
                f" vadr|genes={gene_names}"
                f"|coords={region['coordinates']}"
                f"|strand={strand}"
            )

            write_fasta_record(
                output_handle,
                output_header,
                extracted_sequence,
            )

            extracted += 1

            logger.info(
                "Extracted %s: segment=%s genes=%s region=%s strand=%s length=%d ambiguous=%d",
                seq_name,
                sqc["segment"],
                gene_names,
                region["coordinates"],
                strand,
                extracted_length,
                ambiguous_count,
            )

    return {"total": total, "extracted": extracted, "excluded": excluded}


# ----------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Extract full influenza A segment/gene sequences from VADR output."
    )

    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file.")
    parser.add_argument("--ftr", required=True, help="VADR .ftr file.")
    parser.add_argument("--sqc", required=True, help="VADR .sqc file.")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file.")

    parser.add_argument(
        "--min-len-pct",
        type=float,
        default=0.95,
        help="Minimum len coverage (default: 0.95).",
    )
    parser.add_argument(
        "--max-len-pct",
        type=float,
        default=None,
        help="Maximum sequence length coverage (default: no limit).",
    )
    parser.add_argument(
        "--maxambigs",
        type=int,
        default=0,
        help="Maximum number of ambiguous bases (default: 0).",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Main log file. INFO and ERROR messages are written here.",
    )
    parser.add_argument(
        "--error-log",
        default=None,
        help="Error log file. ERROR messages are written here.",
    )


    args = parser.parse_args()

    # --------------------------------------------------------------
    # Default log names
    # --------------------------------------------------------------
    output_path = Path(args.output)

    # Use pure paths instead of forcing str conversion here
    if args.log is None:
        args.log = output_path.with_suffix(output_path.suffix + ".log")
    else:
        args.log = Path(args.log)

    if args.error_log is None:
        args.error_log = output_path.with_suffix(output_path.suffix + ".errors.log")
    else:
        args.error_log = Path(args.error_log)

    # Initialize the custom logging configuration
    setup_logging(log_file=args.log, error_log_file=args.error_log)

    # --------------------------------------------------------------
    # Log Arguments (Swapped 'logger.info' for 'logger.info')
    # --------------------------------------------------------------
    logger.info("extract_vadr_cds.py")
    logger.info("FASTA: %s", args.fasta)
    logger.info("SQC: %s", args.sqc)
    logger.info("FTR: %s", args.ftr)
    logger.info("Output: %s", args.output)
    logger.info("Minimum length coverage: %s", args.min_len_pct)
    logger.info("Maximum length coverage: %s", args.max_len_pct)
    logger.info("Maximum ambiguous bases: %s", args.maxambigs)
    logger.info("Main log: %s", args.log)
    logger.info("Error log: %s", args.error_log)

    # --------------------------------------------------------------
    # Run extraction
    # --------------------------------------------------------------
    try:
        summary = extract_sequences(
            fasta_file=args.fasta,
            sqc_file=args.sqc,
            ftr_file=args.ftr,
            output_file=args.output,
            min_len_pct=args.min_len_pct,
            max_len_pct=args.max_len_pct,
            maxambigs=args.maxambigs,
        )
    except Exception:
        # Replaced with logger.exception to trace back cleanly through my_app
        logger.exception("Fatal error during extraction.")
        return 1

    logger.info(
        "Final summary: total=%d extracted=%d excluded=%d sequences.",
        summary["total"],
        summary["extracted"],
        summary["excluded"],
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
