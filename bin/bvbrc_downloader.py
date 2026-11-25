#!/usr/bin/env python3
"""
bvbrc_downloader.py — Download raw metadata and sequences from BV-BRC.

Purpose
-------
Download BV-BRC genome **metadata** and **sequences** in JSON form, driven
either by a taxonomy search (taxon ID + optional keywords), by an explicit
list of genome IDs, or by reusing/retrying a previous run.

The output is intentionally "raw": no filtering, no reformatting. The
companion `bvbrc2fasta.py` script consumes these files to produce QC'd
FASTA.

Search modes (mutually exclusive)
---------------------------------
    --taxon  <id>[,<id>...]     Search by NCBI taxon ID (+ optional keyword)
    --genome-id <id>[,<id>...]  Fetch specific genome IDs directly
    --metadata <file.json>      Reuse existing metadata; skip search
    --retry-failed <file.json>  Retry only genomes listed in a failures file

Output files
------------
For a taxon search, the output prefix is:

    <prefix>_<taxid>_<keyword>_<YYYYMMDD>

and the script writes:

    <prefix>_meta_raw.json                 Raw metadata records
    <prefix>_sequences_raw.json            Raw sequence records
    <prefix>_sequences_raw.json.complete   Marker: sequence download finished
    <prefix>_failed_genomes.json           Genomes that failed during the run
    <prefix>_failed_after_retry.json       Genomes still failing after retries
    <prefix>_download_<timestamp>_pid<PID>_<runid>.log

The log file name includes a timestamp, PID, and random suffix so that
concurrent runs with identical search criteria never overwrite each
other's logs.

The three throttling knobs (READ THIS)
--------------------------------------
These are the most commonly confused options:

    --batch N
        Client-side batch size. Two effects:
          * Metadata pagination: how many metadata records to request
            per page.
          * Sequence fetching: how many genome IDs are packed into one
            sequence request.

    --seq-limit N
        Server-side cap on rows returned per sequence request. Passed
        as `limit=N` to the BV-BRC sequence search API.

        IMPORTANT: keep `--seq-limit >= --batch`. If `--seq-limit` is
        smaller than `--batch`, the server silently truncates each
        response and you lose sequences without any error being raised.

    --parallel-batches N
        Concurrency level. Number of sequence requests in flight at the
        same time (ThreadPoolExecutor workers). This is the main
        throughput lever.

Retry / cache semantics
-----------------------
    * Metadata cache: keyed by (taxon list, keywords, keyword logic,
      max records). Cache hits are used unless --no-cache is set.
    * Deterministic output reuse: if `<prefix>_meta_raw.json` already
      exists, the script reuses it instead of re-downloading.
    * Sequence completion marker: if
      `<prefix>_sequences_raw.json.complete` exists, the sequence file
      is treated as complete and skipped.
    * Auto-retry: after the initial run, failed metadata and sequence
      genome IDs are retried up to `--retry-attempts` times, with
      exponential backoff and jitter on 429 / 5xx responses.
    * Manual retry: `--retry-failed <file>` re-attempts only the
      genome IDs listed in that file, in a fresh output prefix.

Author: <add your name here>
License: <add license here>
"""

import argparse
import hashlib
import json
import os
import random
import re
import sys
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import bvbrc as bv


# -----------------------------------------------------------------------------
# BV-BRC clients
# -----------------------------------------------------------------------------
# Two clients are needed:
#   * genome_client    → metadata search / fetch (GenomeClient)
#   * sequence_client  → sequence search / fetch (GenomeSequenceClient)
# Both are instantiated once at import time and reused throughout the run.

genome_client = bv.GenomeClient()
sequence_client = bv.GenomeSequenceClient()


# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------


class Logger:
    """
    Mirror stdout to both the terminal and a unique log file.

    Unlike a stderr tee, this logger intentionally writes to stdout so
    that the log file captures exactly what the user sees on screen.
    Each write opens the log file in append mode with line buffering,
    which is slow but keeps the log readable if the process is killed.
    """

    def __init__(self, log_file):
        self.log_file = log_file
        # sys.__stdout__ is the *original* stdout, captured before any
        # redirection, so we never recurse into ourselves.
        self.console = sys.__stdout__

    def write(self, message):
        """Forward a message to both the console and the log file."""
        self.console.write(message)
        self.console.flush()

        with open(self.log_file, "a", buffering=1) as f:
            f.write(message)

    def flush(self):
        """Flush the console stream (log file is opened per-write)."""
        self.console.flush()


def make_unique_log_file(output_prefix):
    """
    Build a log filename safe for concurrent identical runs.

    Format:
        <prefix>_download_<timestamp>_pid<PID>_<runid>.log

    The timestamp + PID + random suffix guarantees that two processes
    started with the same arguments will never write to the same log.
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pid = os.getpid()
    run_id = uuid.uuid4().hex[:8]

    return (
        f"{output_prefix}_download_"
        f"{timestamp}_pid{pid}_{run_id}.log"
    )


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------


def sanitize_filename_part(value, max_length=60):
    """
    Make a value safe for use in a filename.

    Collapses whitespace to underscores, strips characters that are not
    [A-Za-z0-9_.-], collapses repeated underscores, and truncates to
    `max_length`. Returns "unknown" if nothing usable remains.
    """
    value = str(value).strip()
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    value = value.strip("_.-")

    if not value:
        return "unknown"

    return value[:max_length]


def parse_taxa(taxon_string):
    """
    Parse a comma-separated list of taxon IDs.

    Returns a de-duplicated list preserving the user's original order.
    Raises ValueError on invalid or non-positive IDs, or if the list
    ends up empty.
    """
    taxa = []

    for item in taxon_string.split(","):
        item = item.strip()

        if not item:
            continue

        try:
            taxon_id = int(item)
        except ValueError:
            raise ValueError(
                f"Invalid taxon ID: {item!r}"
            )

        if taxon_id <= 0:
            raise ValueError(
                f"Taxon ID must be > 0: {taxon_id}"
            )

        taxa.append(taxon_id)

    if not taxa:
        raise ValueError(
            "At least one taxon ID is required"
        )

    # Preserve user order but remove duplicates.
    return list(dict.fromkeys(taxa))


def parse_keywords(keyword_str):
    """
    Parse keyword input into a list of tokens.

    If the string contains a boolean operator (AND / OR / NOT), it is
    preserved verbatim as a single expression, because the user has
    already told us how the terms should combine. Otherwise the string
    is split on whitespace and commas.
    """
    if not keyword_str:
        return []

    keyword_str = keyword_str.strip()

    if re.search(
        r"\b(?:OR|NOT|AND)\b",
        keyword_str,
        flags=re.IGNORECASE,
    ):
        return [keyword_str]

    return [
        kw.strip()
        for kw in re.split(r"[\s,]+", keyword_str)
        if kw.strip()
    ]


def build_keyword_query(keywords, combine_with_or=False):
    """
    Build a BV-BRC keyword query from a list of keywords.

    A single keyword becomes a single bv.keyword(...) node. Multiple
    keywords are either joined with OR into one expression, or returned
    as a list of nodes to be AND-ed by the caller (see
    fetch_genomes_paginated).
    """
    if not keywords:
        return None

    if len(keywords) == 1:
        return bv.keyword(keywords[0])

    if combine_with_or:
        return bv.keyword(
            " OR ".join(keywords)
        )

    return [
        bv.keyword(kw)
        for kw in keywords
    ]


def make_output_prefix(
    base_prefix,
    taxon_list,
    keyword,
    date_string=None,
):
    """
    Build the deterministic output prefix:

        <base_prefix>_<taxid>[-<taxid>...]_<keyword>_<YYYYMMDD>

    The keyword component is omitted if no keyword was supplied. The
    date component defaults to today's date.
    """
    date_string = (
        date_string
        or datetime.now().strftime("%Y%m%d")
    )

    taxon_part = "-".join(
        str(t)
        for t in taxon_list
    )

    taxon_part = sanitize_filename_part(
        taxon_part,
        max_length=80,
    )

    parts = [
        sanitize_filename_part(base_prefix)
    ]

    parts.append(taxon_part)

    if keyword:
        keyword_part = sanitize_filename_part(
            keyword,
            max_length=80,
        )
        parts.append(keyword_part)

    parts.append(date_string)

    return "_".join(parts)


# -----------------------------------------------------------------------------
# Cache
# -----------------------------------------------------------------------------


def get_cache_key(
    taxon_list=None,
    keywords=None,
    keyword_logic="and",
    max_records=None,
):
    """
    Generate a stable cache key from the metadata search parameters.

    Any change to the taxon list, keywords, keyword logic, or max
    records produces a different key, so a cached file is never reused
    for a different search.
    """
    parts = [
        "taxon="
        + ",".join(
            str(t)
            for t in sorted(taxon_list or [])
        ),
        "kw=" + (keywords or ""),
        "logic=" + keyword_logic,
        "max="
        + (
            str(max_records)
            if max_records is not None
            else "all"
        ),
    ]

    return hashlib.md5(
        "|".join(parts).encode()
    ).hexdigest()


def load_cached_metadata(cache_file):
    """
    Load cached metadata if the cache file exists.

    Returns the parsed records on success, or None on any failure
    (missing file, unreadable JSON, etc.). Never raises.
    """
    if not os.path.exists(cache_file):
        return None

    try:
        with open(cache_file, "r") as f:
            data = json.load(f)

        print(
            f"\n📂 Loaded cached metadata: {cache_file}"
        )

        print(
            f"   Records: {len(data):,}"
        )

        return data

    except Exception as e:

        print(
            f"   ⚠️ Could not load cache "
            f"{cache_file}: {e}"
        )

        return None


def save_cached_metadata(records, cache_file):
    """
    Save metadata to a cache file atomically.

    Writes to a temporary file (unique per PID + random suffix) and
    then os.replace()s it into place, so a crash mid-write cannot
    leave a partially written cache file.
    """
    try:

        tmp_file = (
            f"{cache_file}.tmp."
            f"{os.getpid()}."
            f"{uuid.uuid4().hex[:8]}"
        )

        with open(tmp_file, "w") as f:
            json.dump(
                records,
                f,
                indent=2,
                default=str,
            )

        os.replace(
            tmp_file,
            cache_file,
        )

        print(
            f"\n💾 Cached metadata saved: "
            f"{cache_file}"
        )

        return True

    except Exception as e:

        print(
            f"   ⚠️ Error saving cache: {e}"
        )

        return False


# -----------------------------------------------------------------------------
# File helpers
# -----------------------------------------------------------------------------


def atomic_json_dump(data, output_file):
    """
    Write JSON to `output_file` atomically.

    A crash mid-write cannot leave a truncated output file: the data
    goes to a temp file first, then os.replace() renames it into place.
    """
    tmp_file = (
        f"{output_file}.tmp."
        f"{os.getpid()}."
        f"{uuid.uuid4().hex[:8]}"
    )

    with open(tmp_file, "w") as f:
        json.dump(
            data,
            f,
            indent=2,
            default=str,
        )

    os.replace(
        tmp_file,
        output_file,
    )


def save_raw_metadata(records, output_file):
    """
    Save metadata records to the given output file and print a summary.
    """
    atomic_json_dump(
        records,
        output_file,
    )

    file_size = (
        os.path.getsize(output_file)
        / (1024 * 1024)
    )

    print(
        f"\n💾 Metadata saved: {output_file}"
    )

    print(
        f"   Records: {len(records):,}"
    )

    print(
        f"   Size: {file_size:.2f} MB"
    )

    return output_file


def save_failed_genomes(
    failed_genomes,
    output_file,
    attempt=0,
):
    """
    Save the list of failed genome IDs to a JSON file.

    De-duplicates the input, sorts the result for stability, and
    records a timestamp + attempt number for traceability.
    """
    if not failed_genomes:
        return

    failed_genomes = list(
        dict.fromkeys(failed_genomes)
    )

    atomic_json_dump(
        {
            "timestamp": datetime.now().isoformat(),
            "failed_genomes": sorted(
                failed_genomes
            ),
            "count": len(failed_genomes),
            "attempt": attempt,
        },
        output_file,
    )

    print(
        f"\n💾 Failed genomes saved: "
        f"{output_file}"
    )

    print(
        f"   Failed genomes: "
        f"{len(failed_genomes):,}"
    )


def load_failed_genomes(input_file):
    """
    Load a failed-genomes JSON file.

    Returns the list of genome IDs, or an empty list on any error.
    """
    try:

        with open(input_file, "r") as f:
            data = json.load(f)

        return data.get(
            "failed_genomes",
            [],
        )

    except Exception as e:

        print(
            f"   ⚠️ Error loading failed genomes: "
            f"{e}"
        )

        return []


def load_metadata_from_file(metadata_file):
    """
    Load a metadata JSON file.

    Returns the parsed list on success, or None on any error.
    """
    try:

        with open(metadata_file, "r") as f:
            records = json.load(f)

        print(
            f"\n📂 Loaded metadata: "
            f"{metadata_file}"
        )

        print(
            f"   Records: {len(records):,}"
        )

        return records

    except Exception as e:

        print(
            f"   ❌ Error loading metadata: "
            f"{e}"
        )

        return None


def load_sequences_from_file(sequence_file):
    """
    Load a sequence JSON file.

    Returns the parsed list on success, or None on any error.
    """
    try:

        with open(sequence_file, "r") as f:
            sequences = json.load(f)

        print(
            f"\n📂 Loaded sequences: "
            f"{sequence_file}"
        )

        print(
            f"   Records: {len(sequences):,}"
        )

        return sequences

    except Exception as e:

        print(
            f"   ⚠️ Could not load sequences: "
            f"{e}"
        )

        return None


# -----------------------------------------------------------------------------
# Genome fetching
# -----------------------------------------------------------------------------


def fetch_genome_with_retry(
    genome_id,
    max_retries=3,
):
    """
    Fetch a single genome record with retry handling.

    Retry policy
    ------------
    * 200 → success, return (True, record, None).
    * 404 → hard failure ("Not found"), no retry.
    * 429 → rate limited; back off linearly (10 * attempt + jitter).
    * 5xx → server error; back off exponentially (2^attempt + jitter).
    * other 4xx → hard failure, no retry.
    * Network exception → exponential backoff, retried up to max_retries.

    Returns
    -------
    tuple
        (success, record_or_None, error_or_None)
    """
    for attempt in range(max_retries):

        try:

            response = genome_client.get(
                genome_id
            )

            if response.status_code == 200:
                return (
                    True,
                    response.json(),
                    None,
                )

            if response.status_code == 404:
                return (
                    False,
                    None,
                    "Not found",
                )

            if response.status_code == 429:

                wait_time = (
                    10 * (attempt + 1)
                    + random.randint(0, 3)
                )

            elif response.status_code >= 500:

                wait_time = (
                    2 ** attempt
                    + random.randint(0, 3)
                )

            else:

                return (
                    False,
                    None,
                    f"HTTP {response.status_code}",
                )

            if attempt < max_retries - 1:
                time.sleep(wait_time)

        except Exception as e:

            if attempt == max_retries - 1:
                return (
                    False,
                    None,
                    str(e),
                )

            time.sleep(
                2 ** attempt
            )

    return (
        False,
        None,
        "Max retries exceeded",
    )


def fetch_genomes_batch_retry(
    failed_genomes,
    batch_size=200,
    max_retries=3,
):
    """
    Retry failed genome metadata in batches.

    Strategy
    --------
    For each batch of up to `batch_size` genome IDs:
        1. Build an OR query: genome_id == gid1 | genome_id == gid2 | ...
        2. Search BV-BRC for that batch with retry / backoff.
        3. If any genome IDs from the batch are missing from the
           response, fall back to fetching them individually.

    If the batch search fails entirely, every genome ID in the batch
    is fetched individually, and any still failing are returned.

    Returns
    -------
    tuple
        (recovered_records, still_failed_ids)
    """
    if not failed_genomes:
        return [], []

    failed_genomes = list(
        dict.fromkeys(failed_genomes)
    )

    all_records = []
    still_failed = []

    total_batches = (
        len(failed_genomes)
        + batch_size
        - 1
    ) // batch_size

    print(
        f"\n📦 Retry metadata batches: "
        f"{len(failed_genomes):,} genomes, "
        f"batch size {batch_size:,}, "
        f"{total_batches} batches"
    )

    for batch_start in range(
        0,
        len(failed_genomes),
        batch_size,
    ):

        batch = failed_genomes[
            batch_start:
            batch_start + batch_size
        ]

        batch_num = (
            batch_start // batch_size
            + 1
        )

        batch_records = []
        success = False

        for attempt in range(max_retries):

            try:

                query = None

                for gid in batch:

                    condition = (
                        genome_client.genome_id
                        == gid
                    )

                    query = (
                        condition
                        if query is None
                        else query | condition
                    )

                response = genome_client.search(
                    query,
                    limit=len(batch),
                    start=0,
                )

                if response.status_code == 200:

                    batch_records = (
                        response.json()
                    )

                    success = True
                    break

                if response.status_code == 429:

                    wait_time = (
                        10 * (attempt + 1)
                        + random.randint(0, 3)
                    )

                elif response.status_code >= 500:

                    wait_time = (
                        2 ** attempt
                        + random.randint(0, 3)
                    )

                else:

                    print(
                        f"      ⚠️ Batch {batch_num}: "
                        f"HTTP {response.status_code}"
                    )

                    break

                if attempt < max_retries - 1:

                    time.sleep(
                        wait_time
                    )

            except Exception as e:

                if attempt == max_retries - 1:

                    print(
                        f"      ⚠️ Batch {batch_num} "
                        f"failed: {e}"
                    )

                else:

                    time.sleep(
                        2 ** attempt
                    )

        # ---------------------------------------------------------------------
        # Reconcile batch result: fetch any genome IDs the batch search
        # didn't return individually, and record the ones that still fail.
        # ---------------------------------------------------------------------
        if success:

            fetched_ids = {
                r.get("genome_id")
                for r in batch_records
                if r.get("genome_id")
            }

            missing_ids = [
                gid
                for gid in batch
                if gid not in fetched_ids
            ]

            for gid in missing_ids:

                ok, record, error = (
                    fetch_genome_with_retry(
                        gid,
                        max_retries=2,
                    )
                )

                if ok and record:

                    batch_records.append(
                        record
                    )

                else:

                    still_failed.append(
                        gid
                    )

            all_records.extend(
                batch_records
            )

        else:

            # Batch search failed entirely → try each genome ID alone.
            for gid in batch:

                ok, record, error = (
                    fetch_genome_with_retry(
                        gid,
                        max_retries=2,
                    )
                )

                if ok and record:

                    all_records.append(
                        record
                    )

                else:

                    still_failed.append(
                        gid
                    )

        print(
            f"   Retry batch "
            f"{batch_num}/{total_batches}: "
            f"recovered "
            f"{len(batch_records):,}, "
            f"failed "
            f"{len(still_failed):,} cumulative"
        )

        # Small politeness delay between batches.
        if (
            batch_start + batch_size
            < len(failed_genomes)
        ):

            time.sleep(0.5)

    return (
        all_records,
        list(dict.fromkeys(still_failed)),
    )


def fetch_genomes_by_genome_ids(
    genome_ids,
    use_cache=True,
    cache_dir="cache",
):
    """
    Fetch genome metadata by explicit genome IDs.

    Uses a cache keyed by the sorted genome ID list. If every requested
    genome ID is present in the cache, the cached records are returned
    without any API calls.

    Returns
    -------
    tuple
        (records, failed_ids)
    """
    os.makedirs(
        cache_dir,
        exist_ok=True,
    )

    genome_ids = list(
        dict.fromkeys(genome_ids)
    )

    cache_key = hashlib.md5(
        ",".join(
            sorted(genome_ids)
        ).encode()
    ).hexdigest()

    cache_file = os.path.join(
        cache_dir,
        f"genome_ids_{cache_key}.json",
    )

    if use_cache:

        cached = load_cached_metadata(
            cache_file
        )

        if cached is not None:

            cached_ids = {
                r.get("genome_id")
                for r in cached
                if r.get("genome_id")
            }

            # Only trust the cache if it covers every requested ID.
            if all(
                gid in cached_ids
                for gid in genome_ids
            ):
                return cached, []

    all_records = []
    failed = []

    print(
        "\n📥 FETCHING GENOMES BY GENOME_ID"
    )

    print(
        f"   Total: {len(genome_ids):,}"
    )

    for i, genome_id in enumerate(
        genome_ids,
        1,
    ):

        print(
            f"   {i}/{len(genome_ids)}: "
            f"{genome_id}"
        )

        ok, record, error = (
            fetch_genome_with_retry(
                genome_id
            )
        )

        if ok:

            all_records.append(
                record
            )

        else:

            failed.append(
                genome_id
            )

            print(
                f"      ❌ Failed: {error}"
            )

    print(
        f"   ✅ Success: "
        f"{len(all_records):,}"
    )

    print(
        f"   ❌ Failed: "
        f"{len(failed):,}"
    )

    if all_records and use_cache:

        save_cached_metadata(
            all_records,
            cache_file,
        )

    return (
        all_records,
        failed,
    )


# -----------------------------------------------------------------------------
# Metadata pagination
# -----------------------------------------------------------------------------


def fetch_genomes_paginated(
    taxon_list,
    keywords=None,
    keyword_logic="and",
    batch_size=1000,
    max_records=None,
    use_cache=True,
    cache_dir="cache",
):
    """
    Fetch raw genome metadata by taxon ID using pagination.

    For each taxon ID:
        * Build a query: taxon_id == <id> [AND <keyword expression>]
        * Page through results with `limit=batch_size` and increasing
          `start=`, until:
            - an empty page is returned,
            - the server returns fewer than `batch_size` rows (last page),
            - or `max_records` has been reached for this taxon.

    The result is cached only if no page failed; a partial result is
    never cached, so a subsequent run will refetch the missing data.

    Returns
    -------
    tuple
        (all_records, []) — the second element is always an empty list
        for API compatibility with the other fetchers.
    """
    os.makedirs(
        cache_dir,
        exist_ok=True,
    )

    cache_key = get_cache_key(
        taxon_list,
        keywords,
        keyword_logic,
        max_records,
    )

    cache_file = os.path.join(
        cache_dir,
        f"metadata_{cache_key}.json",
    )

    if use_cache:

        cached = load_cached_metadata(
            cache_file
        )

        if cached is not None:
            return cached, []

    keyword_list = (
        parse_keywords(keywords)
        if keywords
        else []
    )

    keyword_query = build_keyword_query(
        keyword_list,
        keyword_logic == "or",
    )

    all_records = []
    failed_pages = []

    max_retries = 3

    print(
        "\n📥 FETCHING RAW GENOME DATA"
    )

    print(
        "   Search by: taxon_id"
    )

    print(
        "   Taxa: "
        + ", ".join(
            str(t)
            for t in taxon_list
        )
    )

    if keyword_list:

        print(
            "   Keywords: "
            + ", ".join(keyword_list)
        )

        print(
            f"   Keyword logic: "
            f"{keyword_logic.upper()}"
        )

    print(
        f"   Metadata/sequence batch size: "
        f"{batch_size:,}"
    )

    print(
        f"   Max records per taxon: "
        f"{max_records if max_records else 'ALL'}"
    )

    print("-" * 60)

    for taxon_id in taxon_list:

        print(
            f"\n🔄 Fetching taxon_id={taxon_id}"
        )

        start = 0
        term_fetched = 0
        page_num = 0
        done = False

        while not done:

            if (
                max_records is not None
                and term_fetched >= max_records
            ):
                break

            # Never request more than the user's per-taxon cap.
            current_limit = batch_size

            if max_records is not None:

                current_limit = min(
                    batch_size,
                    max_records
                    - term_fetched,
                )

            page_num += 1
            success = False

            for attempt in range(
                max_retries
            ):

                try:

                    query = (
                        genome_client.taxon_id
                        == taxon_id
                    )

                    if keyword_query:

                        if isinstance(
                            keyword_query,
                            list,
                        ):

                            for kw in keyword_query:
                                query = (
                                    query & kw
                                )

                        else:

                            query = (
                                query
                                & keyword_query
                            )

                    print(
                        f"   Page {page_num}: "
                        f"start={start}, "
                        f"limit={current_limit}"
                    )

                    response = (
                        genome_client.search(
                            query,
                            limit=current_limit,
                            start=start,
                        )
                    )

                    if response.status_code == 200:

                        success = True

                        batch_data = (
                            response.json()
                        )

                        batch_count = len(
                            batch_data
                        )

                        if batch_count == 0:

                            done = True
                            break

                        all_records.extend(
                            batch_data
                        )

                        term_fetched += (
                            batch_count
                        )

                        print(
                            f"   ✅ Fetched "
                            f"{batch_count:,} "
                            f"(taxon total: "
                            f"{term_fetched:,}, "
                            f"overall: "
                            f"{len(all_records):,})"
                        )

                        start += batch_count

                        # A short page means we've hit the end.
                        if (
                            batch_count
                            < current_limit
                        ):

                            done = True

                        elif (
                            max_records
                            is not None
                            and term_fetched
                            >= max_records
                        ):

                            done = True

                        break

                    if response.status_code == 429:

                        wait_time = (
                            10 * (attempt + 1)
                            + random.randint(0, 3)
                        )

                    elif response.status_code >= 500:

                        wait_time = (
                            2 ** attempt
                            + random.randint(0, 5)
                        )

                    else:

                        print(
                            f"   ❌ HTTP "
                            f"{response.status_code}; "
                            f"stopping this taxon"
                        )

                        done = True
                        break

                    if attempt < max_retries - 1:

                        print(
                            f"   ⚠️ HTTP "
                            f"{response.status_code}; "
                            f"retry "
                            f"{attempt + 1}/"
                            f"{max_retries} "
                            f"in {wait_time}s"
                        )

                        time.sleep(
                            wait_time
                        )

                except Exception as e:

                    if attempt < max_retries - 1:

                        wait_time = (
                            2 ** attempt
                            + random.randint(0, 5)
                        )

                        print(
                            f"   ⚠️ Exception: "
                            f"{e}; retry "
                            f"{attempt + 1}/"
                            f"{max_retries} "
                            f"in {wait_time}s"
                        )

                        time.sleep(
                            wait_time
                        )

                    else:

                        print(
                            f"   ❌ Exception after "
                            f"{max_retries} "
                            f"attempts: {e}"
                        )

            if (
                not success
                and not done
            ):

                failed_pages.append(
                    (
                        taxon_id,
                        page_num,
                    )
                )

                done = True

            if not done:
                time.sleep(0.5)

    print("-" * 60)

    print(
        f"✅ Total metadata records: "
        f"{len(all_records):,}"
    )

    if failed_pages:

        print(
            f"⚠️ Failed pages: "
            f"{failed_pages}"
        )

    # Never cache an incomplete pagination result.
    if (
        all_records
        and use_cache
        and not failed_pages
    ):

        save_cached_metadata(
            all_records,
            cache_file,
        )

    return (
        all_records,
        [],
    )


# -----------------------------------------------------------------------------
# Sequence fetching
# -----------------------------------------------------------------------------


def fetch_sequences_parallel_batch(
    records,
    output_file,
    batch_size=100,
    parallel_batches=2,
    seq_limit=5000,
    max_retries=5,
    genome_ids_filter=None,
):
    """
    Fetch genome sequences from BV-BRC, using batched concurrent requests.

    Throttling semantics (see module docstring for the full story)
    --------------------------------------------------------------
    * `batch_size`       → number of genome IDs packed into one request.
    * `parallel_batches` → number of requests in flight concurrently.
    * `seq_limit`        → server-side row cap per request.
                           MUST be >= batch_size, otherwise the server
                           silently truncates each response.

    If `genome_ids_filter` is provided, only those genome IDs are
    fetched (used for retrying only the failed subset).

    Returns
    -------
    tuple
        (sequence_count, failed_genome_ids)
    """
    all_genome_ids = list(
        dict.fromkeys(
            r.get("genome_id")
            for r in records
            if r.get("genome_id")
        )
    )

    if genome_ids_filter is not None:

        wanted = set(
            genome_ids_filter
        )

        genome_ids = [
            gid
            for gid in all_genome_ids
            if gid in wanted
        ]

    else:

        genome_ids = all_genome_ids

    if not genome_ids:

        # Ensure the output file exists even when there's nothing to do,
        # so downstream consumers can rely on its presence.
        if not os.path.exists(output_file):

            atomic_json_dump(
                [],
                output_file,
            )

        print(
            "\n⚠️ No genome IDs to fetch."
        )

        return 0, []

    batches = [
        genome_ids[i:i + batch_size]
        for i in range(
            0,
            len(genome_ids),
            batch_size,
        )
    ]

    print(
        "\n🧬 FETCHING SEQUENCES"
    )

    print("=" * 60)

    print(
        f"   Genomes: {len(genome_ids):,}"
    )

    print(
        f"   Sequence batch size: "
        f"{batch_size:,}"
    )

    print(
        f"   Parallel requests: "
        f"{parallel_batches}"
    )

    print(
        f"   Sequence limit/request: "
        f"{seq_limit:,}"
    )

    print(
        f"   Max retries: "
        f"{max_retries}"
    )

    print("=" * 60)

    start_time = time.time()

    total_batches = len(batches)

    all_sequences = []
    failed_batches = []

    def fetch_batch(batch_data):
        """
        Worker: fetch one batch of genome IDs.

        Retry policy matches the genome fetcher: 429 and 5xx are
        retried with exponential backoff + jitter, 4xx (other than
        429) are hard failures. Returns a tuple:
            (batch_idx, sequences, error_or_None)
        """
        batch_idx, batch = batch_data

        query = None

        for gid in batch:

            condition = (
                sequence_client.genome_id
                == gid
            )

            query = (
                condition
                if query is None
                else query | condition
            )

        for attempt in range(
            max_retries
        ):

            try:

                response = (
                    sequence_client.search(
                        query,
                        limit=seq_limit,
                        start=0,
                    )
                )

                status = response.status_code

                if status == 200:

                    data = response.json()

                    return (
                        batch_idx,
                        data,
                        None,
                    )

                if status == 429:

                    if (
                        attempt
                        >= max_retries - 1
                    ):

                        return (
                            batch_idx,
                            [],
                            f"HTTP 429 after "
                            f"{max_retries} attempts",
                        )

                    wait_time = min(
                        60,
                        (2 ** attempt) * 5
                        + random.uniform(0, 3),
                    )

                    time.sleep(
                        wait_time
                    )

                    continue

                if 500 <= status < 600:

                    if (
                        attempt
                        >= max_retries - 1
                    ):

                        return (
                            batch_idx,
                            [],
                            f"HTTP {status} after "
                            f"{max_retries} attempts",
                        )

                    wait_time = min(
                        60,
                        (2 ** attempt) * 3
                        + random.uniform(0, 3),
                    )

                    time.sleep(
                        wait_time
                    )

                    continue

                if 400 <= status < 500:

                    return (
                        batch_idx,
                        [],
                        f"HTTP {status} "
                        f"(non-retryable)",
                    )

                return (
                    batch_idx,
                    [],
                    f"HTTP {status}",
                )

            except Exception as e:

                if (
                    attempt
                    >= max_retries - 1
                ):

                    return (
                        batch_idx,
                        [],
                        f"{type(e).__name__}: {e}",
                    )

                wait_time = min(
                    60,
                    (2 ** attempt) * 3
                    + random.uniform(0, 3),
                )

                time.sleep(
                    wait_time
                )

        return (
            batch_idx,
            [],
            f"Failed after "
            f"{max_retries} attempts",
        )

    # -------------------------------------------------------------------------
    # Controlled concurrency: submit all batches, drain as they finish.
    # -------------------------------------------------------------------------

    with ThreadPoolExecutor(
        max_workers=parallel_batches
    ) as executor:

        futures = {
            executor.submit(
                fetch_batch,
                (idx, batch),
            ): idx
            for idx, batch
            in enumerate(batches)
        }

        completed = 0

        for future in as_completed(
            futures
        ):

            batch_idx = futures[future]

            completed += 1

            try:

                (
                    _,
                    sequences,
                    error,
                ) = future.result()

                if error:

                    failed_batches.append(
                        (
                            batch_idx,
                            error,
                        )
                    )

                else:

                    all_sequences.extend(
                        sequences
                    )

            except Exception as e:

                failed_batches.append(
                    (
                        batch_idx,
                        f"{type(e).__name__}: {e}",
                    )
                )

            # Progress line every 10 batches (and on the last one).
            if (
                completed % 10 == 0
                or completed == total_batches
            ):

                elapsed = (
                    time.time()
                    - start_time
                )

                rate = (
                    completed / elapsed
                    if elapsed
                    else 0
                )

                print(
                    f"   📊 "
                    f"{completed}/"
                    f"{total_batches} batches "
                    f"({len(all_sequences):,} "
                    f"sequences, "
                    f"{rate:.2f} "
                    f"batches/sec)"
                )

    # -------------------------------------------------------------------------
    # Collect failed genome IDs from failed batches.
    # -------------------------------------------------------------------------

    failed_genomes = []

    for batch_idx, error in failed_batches:

        failed_genomes.extend(
            batches[batch_idx]
        )

        print(
            f"   ⚠️ Batch "
            f"{batch_idx + 1}: "
            f"{error}"
        )

    failed_genomes = list(
        dict.fromkeys(
            failed_genomes
        )
    )

    # -------------------------------------------------------------------------
    # Save sequences and report.
    # -------------------------------------------------------------------------

    atomic_json_dump(
        all_sequences,
        output_file,
    )

    elapsed = (
        time.time()
        - start_time
    )

    file_size = (
        os.path.getsize(output_file)
        / (1024 * 1024)
    )

    print(
        f"\n✅ Sequence download finished "
        f"in {elapsed:.1f}s"
    )

    print(
        f"💾 {output_file}"
    )

    print(
        f"   Records: "
        f"{len(all_sequences):,}"
    )

    print(
        f"   Size: "
        f"{file_size:.2f} MB"
    )

    print(
        f"   Failed genomes: "
        f"{len(failed_genomes):,}"
    )

    return (
        len(all_sequences),
        failed_genomes,
    )


# -----------------------------------------------------------------------------
# Sequence merge / retry
# -----------------------------------------------------------------------------


def merge_sequence_files(
    existing_file,
    new_sequences,
):
    """
    Merge sequence records by genome_id.

    New records replace existing records with the same genome_id, so
    retries can safely overwrite stale or partial entries.
    """
    existing = []

    if os.path.exists(existing_file):

        loaded = load_sequences_from_file(
            existing_file
        )

        if loaded is not None:
            existing = loaded

    merged = {}

    for record in existing:

        gid = record.get("genome_id")

        if gid:
            merged[gid] = record

    for record in new_sequences:

        gid = record.get("genome_id")

        if gid:
            merged[gid] = record

    return list(
        merged.values()
    )


def fetch_and_merge_sequence_retry(
    genome_records,
    failed_genomes,
    output_file,
    batch_size,
    parallel_batches,
    seq_limit,
    max_retries,
):
    """
    Retry only failed sequence genome IDs and merge recovered records
    into the existing sequence file.

    Uses a temporary retry file so a crash mid-retry cannot corrupt
    the main sequence file; the merge happens only after the retry
    file has been fully written.

    Returns
    -------
    tuple
        (total_sequences_in_output, still_failed_ids)
    """
    if not failed_genomes:

        existing = load_sequences_from_file(
            output_file
        )

        return (
            len(existing)
            if existing is not None
            else 0,
            [],
        )

    print(
        "\n🔄 RETRYING FAILED SEQUENCES"
    )

    print(
        f"   Genome IDs to retry: "
        f"{len(failed_genomes):,}"
    )

    retry_file = (
        f"{output_file}.retry."
        f"{os.getpid()}."
        f"{uuid.uuid4().hex[:8]}.json"
    )

    try:

        _, still_failed = (
            fetch_sequences_parallel_batch(
                genome_records,
                retry_file,
                batch_size=batch_size,
                parallel_batches=parallel_batches,
                seq_limit=seq_limit,
                max_retries=max_retries,
                genome_ids_filter=failed_genomes,
            )
        )

        retry_sequences = (
            load_sequences_from_file(
                retry_file
            )
        )

        if retry_sequences is None:
            retry_sequences = []

        merged = merge_sequence_files(
            output_file,
            retry_sequences,
        )

        atomic_json_dump(
            merged,
            output_file,
        )

        print(
            f"   Recovered sequences: "
            f"{len(retry_sequences):,}"
        )

        print(
            f"   Still failed: "
            f"{len(still_failed):,}"
        )

        return (
            len(merged),
            still_failed,
        )

    finally:

        # Always clean up the temp retry file, even on exception.
        if os.path.exists(retry_file):

            try:
                os.remove(
                    retry_file
                )
            except OSError:
                pass


# -----------------------------------------------------------------------------
# Parser
# -----------------------------------------------------------------------------


def build_parser():
    """
    Build the argparse parser for the downloader.

    Search modes (--taxon / --genome-id / --metadata / --retry-failed)
    are mutually exclusive and one is required.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Download raw metadata and "
            "sequences from BV-BRC by "
            "taxon ID."
        ),
        formatter_class=(
            argparse.RawDescriptionHelpFormatter
        ),
        epilog="""
EXAMPLES:

  # Taxon + keyword
  bvbrc_downloader.py \\
      --taxon 11320 \\
      --keyword "Influenza A" \\
      --batch 1000

  # Multiple taxa
  bvbrc_downloader.py \\
      --taxon 11320,11520 \\
      --keyword influenza \\
      --batch 1000

  # Reuse existing metadata
  bvbrc_downloader.py \\
      --metadata bvbrc_11320_influenza_20260828_meta_raw.json

  # Retry failed genomes
  bvbrc_downloader.py \\
      --retry-failed bvbrc_11320_influenza_20260828_failed_genomes.json

OUTPUT:

  <prefix>_meta_raw.json
  <prefix>_sequences_raw.json
  <prefix>_failed_genomes.json
  <prefix>_failed_after_retry.json
  <prefix>_download_<timestamp>_pid<PID>_<runid>.log

The output prefix for a taxon search is:

  <prefix>_<taxid>_<keyword>_<YYYYMMDD>

IMPORTANT — THROTTLING:

  --batch controls BOTH metadata page size
  AND sequence batch size (genome IDs per request).

  --seq-limit is the SERVER-SIDE cap on rows returned
  per sequence request. Keep --seq-limit >= --batch,
  otherwise the server silently truncates responses
  and you lose sequences without any error.

  --parallel-batches controls how many sequence
  requests are sent at the same time.
""",
    )

    # -------------------------------------------------------------------------
    # Search/input modes (mutually exclusive, one required)
    # -------------------------------------------------------------------------

    search_group = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    search_group.add_argument(
        "--taxon",
        "-t",
        type=str,
        help="Taxon ID(s), comma-separated",
    )

    search_group.add_argument(
        "--genome-id",
        type=str,
        help="Genome ID(s), comma-separated",
    )

    search_group.add_argument(
        "--metadata",
        "-m",
        type=str,
        help=(
            "Use an existing metadata JSON "
            "file; skip metadata search"
        ),
    )

    search_group.add_argument(
        "--retry-failed",
        type=str,
        help=(
            "Retry failed genomes from a "
            "JSON file"
        ),
    )

    # -------------------------------------------------------------------------
    # Search parameters
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--keyword",
        "-k",
        type=str,
        default=None,
        help="Keyword filter",
    )

    parser.add_argument(
        "--keyword-logic",
        "-l",
        choices=["and", "or"],
        default="and",
        help=(
            "How multiple keywords are "
            "combined (default: and)"
        ),
    )

    parser.add_argument(
        "--max",
        type=int,
        default=None,
        help="Maximum records per taxon",
    )

    parser.add_argument(
        "--batch",
        "-b",
        type=int,
        default=200,
        help=(
            "Metadata page size AND sequence "
            "batch size (genome IDs per request). "
            "Default: 200"
        ),
    )

    parser.add_argument(
        "--parallel-batches",
        type=int,
        default=2,
        help=(
            "Number of parallel sequence "
            "requests (default: 2)"
        ),
    )

    parser.add_argument(
        "--seq-limit",
        type=int,
        default=200,
        help=(
            "Server-side cap on sequences returned "
            "per request. Keep >= --batch. "
            "Default: 200"
        ),
    )

    parser.add_argument(
        "--retry-attempts",
        type=int,
        default=3,
        help=(
            "Automatic retry attempts "
            "(default: 3)"
        ),
    )

    parser.add_argument(
        "--prefix",
        "-p",
        type=str,
        default="bvbrc",
        help=(
            "Base output prefix "
            "(default: bvbrc)"
        ),
    )

    parser.add_argument(
        "--inspect",
        "-i",
        action="store_true",
        help=(
            "Print first metadata record "
            "and exit"
        ),
    )

    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Skip metadata cache",
    )

    parser.add_argument(
        "--cache-dir",
        type=str,
        default="cache",
        help=(
            "Metadata cache directory "
            "(default: cache)"
        ),
    )

    return parser


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main():
    """
    Entry point.

    Responsibilities:
        1. Parse and validate CLI arguments.
        2. Derive the output prefix and unique log filename.
        3. Redirect stdout to the Logger (console + log file).
        4. Run one of:
             * manual retry mode (--retry-failed), or
             * normal download mode (taxon / genome-id / metadata).
        5. Perform automatic retries for any failures.
        6. Write a final failed-genomes file if any remain.
        7. Print a summary and restore stdout.
    """
    parser = build_parser()
    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # Validate arguments
    # -------------------------------------------------------------------------

    if args.batch < 1:
        parser.error(
            "--batch must be >= 1"
        )

    if args.parallel_batches < 1:
        parser.error(
            "--parallel-batches must be >= 1"
        )

    if args.seq_limit < 1:
        parser.error(
            "--seq-limit must be >= 1"
        )

    if args.retry_attempts < 1:
        parser.error(
            "--retry-attempts must be >= 1"
        )

    if (
        args.max is not None
        and args.max < 1
    ):
        parser.error(
            "--max must be >= 1"
        )

    # -------------------------------------------------------------------------
    # Determine search information and output prefix
    # -------------------------------------------------------------------------

    if args.taxon:

        try:

            taxon_list = parse_taxa(
                args.taxon
            )

        except ValueError as e:

            parser.error(
                str(e)
            )

        output_prefix = make_output_prefix(
            args.prefix,
            taxon_list,
            args.keyword,
        )

    elif args.metadata:

        taxon_list = []

        base = os.path.splitext(
            os.path.basename(
                args.metadata
            )
        )[0]

        output_prefix = (
            f"{args.prefix}_"
            f"{sanitize_filename_part(base, 100)}"
        )

    elif args.genome_id:

        taxon_list = []

        output_prefix = (
            f"{args.prefix}_genomeid_"
            f"{sanitize_filename_part(args.genome_id, 80)}"
        )

    else:

        taxon_list = []

        base = os.path.splitext(
            os.path.basename(
                args.retry_failed
            )
        )[0]

        output_prefix = (
            f"{args.prefix}_"
            f"{sanitize_filename_part(base, 100)}"
        )

    # -------------------------------------------------------------------------
    # Unique log per run (timestamp + PID + random suffix)
    # -------------------------------------------------------------------------

    log_file = make_unique_log_file(
        output_prefix
    )

    logger = Logger(
        log_file
    )

    original_stdout = sys.stdout

    sys.stdout = logger

    try:

        print("=" * 80)
        print("🦠 BV-BRC DATA DOWNLOADER")
        print("=" * 80)

        print(
            f"   Start time: "
            f"{datetime.now().isoformat(timespec='seconds')}"
        )

        print(
            f"   PID: {os.getpid()}"
        )

        print(
            f"   Log: {log_file}"
        )

        # =====================================================================
        # MANUAL RETRY MODE
        # =====================================================================
        # Re-fetch only the genome IDs listed in a previous failures file.
        # Writes to a fresh "_retry_*" prefix so it can't clobber the
        # original run's outputs.

        if args.retry_failed:

            print(
                "\n🔄 MANUAL RETRY"
            )

            print(
                f"   Failed file: "
                f"{args.retry_failed}"
            )

            failed_genomes = (
                load_failed_genomes(
                    args.retry_failed
                )
            )

            if not failed_genomes:

                print(
                    "   No failed genomes "
                    "to retry."
                )

                return

            # -----------------------------------------------------------------
            # Retry metadata first
            # -----------------------------------------------------------------

            retry_records, still_failed = (
                fetch_genomes_batch_retry(
                    failed_genomes,
                    batch_size=args.batch,
                    max_retries=3,
                )
            )

            if not retry_records:

                print(
                    "   ❌ No genomes recovered."
                )

                return

            metadata_file = (
                f"{output_prefix}_retry_meta_raw.json"
            )

            sequences_file = (
                f"{output_prefix}_retry_sequences_raw.json"
            )

            save_raw_metadata(
                retry_records,
                metadata_file,
            )

            # -----------------------------------------------------------------
            # Fetch sequences for recovered metadata
            # -----------------------------------------------------------------

            (
                sequence_count,
                seq_failures,
            ) = fetch_sequences_parallel_batch(
                retry_records,
                sequences_file,
                batch_size=args.batch,
                parallel_batches=args.parallel_batches,
                seq_limit=args.seq_limit,
            )

            all_failures = list(
                dict.fromkeys(
                    still_failed
                    + seq_failures
                )
            )

            if all_failures:

                save_failed_genomes(
                    all_failures,
                    (
                        f"{output_prefix}_"
                        f"retry_failed_genomes.json"
                    ),
                )

            print(
                "\n📊 MANUAL RETRY SUMMARY"
            )

            print(
                f"   Metadata recovered: "
                f"{len(retry_records):,}"
            )

            print(
                f"   Sequences recovered: "
                f"{sequence_count:,}"
            )

            print(
                f"   Still failed: "
                f"{len(all_failures):,}"
            )

            return

        # =====================================================================
        # NORMAL MODE
        # =====================================================================

        print(
            "\n⚙️ CONFIGURATION"
        )

        if args.taxon:

            print(
                f"   Taxon(s): "
                f"{', '.join(str(t) for t in taxon_list)}"
            )

            print(
                f"   Keyword: "
                f"{args.keyword or 'NONE'}"
            )

            print(
                f"   Keyword logic: "
                f"{args.keyword_logic.upper()}"
            )

        elif args.metadata:

            print(
                f"   Existing metadata: "
                f"{args.metadata}"
            )

        elif args.genome_id:

            print(
                f"   Genome ID(s): "
                f"{args.genome_id}"
            )

        print(
            f"   Batch size: "
            f"{args.batch:,}"
        )

        print(
            f"   Parallel sequence requests: "
            f"{args.parallel_batches}"
        )

        print(
            f"   Sequence limit/request: "
            f"{args.seq_limit:,}"
        )

        print(
            f"   Max records/taxon: "
            f"{args.max or 'ALL'}"
        )

        print(
            f"   Retry attempts: "
            f"{args.retry_attempts}"
        )

        print(
            f"   Cache enabled: "
            f"{not args.no_cache}"
        )

        print(
            f"   Output prefix: "
            f"{output_prefix}"
        )

        print(
            f"   Log file: "
            f"{log_file}"
        )

        failed_genomes = []

        # ---------------------------------------------------------------------
        # Metadata acquisition
        # ---------------------------------------------------------------------

        metadata_file = (
            f"{output_prefix}_meta_raw.json"
        )

        genome_records = None

        # Reuse the deterministic metadata file if it already exists,
        # unless --no-cache was explicitly requested.
        if (
            os.path.exists(metadata_file)
            and not args.no_cache
        ):

            print(
                f"\n📂 Existing output metadata "
                f"found: {metadata_file}"
            )

            print(
                "   Skipping metadata download "
                "and continuing to sequence fetching."
            )

            genome_records = (
                load_metadata_from_file(
                    metadata_file
                )
            )

            if genome_records is None:

                print(
                    "   ⚠️ Existing metadata "
                    "could not be loaded; "
                    "downloading again."
                )

        # ---------------------------------------------------------------------
        # Existing metadata supplied by user
        # ---------------------------------------------------------------------

        if (
            genome_records is None
            and args.metadata
        ):

            genome_records = (
                load_metadata_from_file(
                    args.metadata
                )
            )

            if genome_records is None:

                raise RuntimeError(
                    "Failed to load existing "
                    "metadata file"
                )

            save_raw_metadata(
                genome_records,
                metadata_file,
            )

        # ---------------------------------------------------------------------
        # Genome IDs supplied directly
        # ---------------------------------------------------------------------

        elif (
            genome_records is None
            and args.genome_id
        ):

            genome_ids = [
                g.strip()
                for g in args.genome_id.split(",")
                if g.strip()
            ]

            (
                genome_records,
                failed_genomes,
            ) = fetch_genomes_by_genome_ids(
                genome_ids,
                use_cache=not args.no_cache,
                cache_dir=args.cache_dir,
            )

            if genome_records:

                save_raw_metadata(
                    genome_records,
                    metadata_file,
                )

        # ---------------------------------------------------------------------
        # Taxon search
        # ---------------------------------------------------------------------

        elif (
            genome_records is None
            and args.taxon
        ):

            (
                genome_records,
                failed_genomes,
            ) = fetch_genomes_paginated(
                taxon_list=taxon_list,
                keywords=args.keyword,
                keyword_logic=args.keyword_logic,
                batch_size=args.batch,
                max_records=args.max,
                use_cache=not args.no_cache,
                cache_dir=args.cache_dir,
            )

            if genome_records:

                save_raw_metadata(
                    genome_records,
                    metadata_file,
                )

        # ---------------------------------------------------------------------
        # No metadata obtained → bail out
        # ---------------------------------------------------------------------

        if not genome_records:

            print(
                "\n❌ No records retrieved."
            )

            if failed_genomes:

                save_failed_genomes(
                    failed_genomes,
                    (
                        f"{output_prefix}_"
                        f"failed_genomes.json"
                    ),
                )

            raise RuntimeError(
                "No metadata records available"
            )

        # ---------------------------------------------------------------------
        # Inspect mode: print the first record and exit
        # ---------------------------------------------------------------------

        if args.inspect:

            print(
                "\n🔍 FIRST RECORD"
            )

            print(
                json.dumps(
                    genome_records[0],
                    indent=2,
                    default=str,
                )[:5000]
            )

            return

        # ---------------------------------------------------------------------
        # Save any initial metadata failures for traceability
        # ---------------------------------------------------------------------

        if failed_genomes:

            save_failed_genomes(
                failed_genomes,
                (
                    f"{output_prefix}_"
                    f"failed_genomes.json"
                ),
            )

        # ---------------------------------------------------------------------
        # Sequence acquisition
        # ---------------------------------------------------------------------

        sequences_file = (
            f"{output_prefix}_sequences_raw.json"
        )

        # A separate ".complete" marker distinguishes "download finished
        # successfully" from "file exists but may be partial".
        sequence_complete_file = (
            f"{sequences_file}.complete"
        )

        sequence_count = 0
        seq_failures = []

        # Reuse completed sequence file if the marker is present.
        if (
            os.path.exists(sequences_file)
            and os.path.exists(
                sequence_complete_file
            )
            and not args.no_cache
        ):

            print(
                f"\n📂 Existing completed "
                f"sequence file found: "
                f"{sequences_file}"
            )

            existing_sequences = (
                load_sequences_from_file(
                    sequences_file
                )
            )

            if (
                existing_sequences
                is not None
            ):

                sequence_count = len(
                    existing_sequences
                )

                print(
                    f"   Records: "
                    f"{sequence_count:,}"
                )

                print(
                    "   Skipping sequence download."
                )

            else:

                print(
                    "   ⚠️ Existing sequence "
                    "file could not be loaded."
                )

                sequence_count, seq_failures = (
                    fetch_sequences_parallel_batch(
                        genome_records,
                        sequences_file,
                        batch_size=args.batch,
                        parallel_batches=args.parallel_batches,
                        seq_limit=args.seq_limit,
                    )
                )

        else:

            sequence_count, seq_failures = (
                fetch_sequences_parallel_batch(
                    genome_records,
                    sequences_file,
                    batch_size=args.batch,
                    parallel_batches=args.parallel_batches,
                    seq_limit=args.seq_limit,
                )
            )

        # ---------------------------------------------------------------------
        # Write or remove the completion marker
        # ---------------------------------------------------------------------

        if not seq_failures:

            with open(
                sequence_complete_file,
                "w",
            ) as f:

                f.write(
                    "Completed: "
                    f"{datetime.now().isoformat()}\n"
                )

            print(
                f"   ✅ Sequence completion "
                f"marker: "
                f"{sequence_complete_file}"
            )

        else:

            # Remove stale completion marker if one exists, since the
            # sequence file is now known to be incomplete.
            if os.path.exists(
                sequence_complete_file
            ):

                try:
                    os.remove(
                        sequence_complete_file
                    )
                except OSError:
                    pass

            print(
                f"   ⚠️ Sequence download "
                f"incomplete; "
                f"{len(seq_failures):,} "
                f"genomes failed."
            )

        # =====================================================================
        # AUTOMATIC RETRY
        # =====================================================================
        # Both metadata and sequence failures are retried up to
        # (retry_attempts - 1) additional times, with a wait between
        # attempts.

        current_metadata_failures = list(
            dict.fromkeys(
                failed_genomes
            )
        )

        current_sequence_failures = list(
            dict.fromkeys(
                seq_failures
            )
        )

        for attempt in range(
            1,
            args.retry_attempts,
        ):

            if (
                not current_metadata_failures
                and not current_sequence_failures
            ):
                break

            print(
                "\n" + "=" * 80
            )

            print(
                f"🔄 AUTO-RETRY ATTEMPT "
                f"{attempt}/"
                f"{args.retry_attempts - 1}"
            )

            print(
                "=" * 80
            )

            # -----------------------------------------------------------------
            # Retry metadata failures
            # -----------------------------------------------------------------

            if current_metadata_failures:

                print(
                    f"\n📋 Metadata failures: "
                    f"{len(current_metadata_failures):,}"
                )

                (
                    recovered_records,
                    still_metadata_failed,
                ) = fetch_genomes_batch_retry(
                    current_metadata_failures,
                    batch_size=args.batch,
                    max_retries=3,
                )

                existing_ids = {
                    r.get("genome_id")
                    for r in genome_records
                    if r.get("genome_id")
                }

                for record in recovered_records:

                    gid = record.get(
                        "genome_id"
                    )

                    if (
                        gid
                        and gid not in existing_ids
                    ):

                        genome_records.append(
                            record
                        )

                        existing_ids.add(
                            gid
                        )

                current_metadata_failures = (
                    list(
                        dict.fromkeys(
                            still_metadata_failed
                        )
                    )
                )

                print(
                    f"   Recovered metadata: "
                    f"{len(recovered_records):,}"
                )

                print(
                    f"   Still failed metadata: "
                    f"{len(current_metadata_failures):,}"
                )

                # Refresh the deterministic metadata file so a
                # subsequent run picks up the recovered records.
                save_raw_metadata(
                    genome_records,
                    metadata_file,
                )

            # -----------------------------------------------------------------
            # Retry sequence failures
            # -----------------------------------------------------------------

            if current_sequence_failures:

                print(
                    f"\n🧬 Sequence failures: "
                    f"{len(current_sequence_failures):,}"
                )

                (
                    sequence_count,
                    current_sequence_failures,
                ) = fetch_and_merge_sequence_retry(
                    genome_records,
                    current_sequence_failures,
                    sequences_file,
                    batch_size=args.batch,
                    parallel_batches=args.parallel_batches,
                    seq_limit=args.seq_limit,
                    max_retries=5,
                )

                if not current_sequence_failures:

                    with open(
                        sequence_complete_file,
                        "w",
                    ) as f:

                        f.write(
                            "Completed: "
                            f"{datetime.now().isoformat()}\n"
                        )

                    print(
                        "   ✅ All sequence "
                        "failures recovered."
                    )

                else:

                    print(
                        f"   ⚠️ Still failed "
                        f"sequences: "
                        f"{len(current_sequence_failures):,}"
                    )

            # -----------------------------------------------------------------
            # Wait before the next attempt (linear backoff)
            # -----------------------------------------------------------------

            if (
                current_metadata_failures
                or current_sequence_failures
            ):

                if attempt < (
                    args.retry_attempts - 1
                ):

                    wait_time = 10 * attempt

                    print(
                        f"\n   ⏳ Waiting "
                        f"{wait_time}s before "
                        f"next retry..."
                    )

                    time.sleep(
                        wait_time
                    )

        # =====================================================================
        # Final failures
        # =====================================================================

        all_failures = list(
            dict.fromkeys(
                current_metadata_failures
                + current_sequence_failures
            )
        )

        if all_failures:

            final_failed_file = (
                f"{output_prefix}_"
                f"failed_after_retry.json"
            )

            save_failed_genomes(
                all_failures,
                final_failed_file,
                attempt=args.retry_attempts,
            )

        # =====================================================================
        # Final summary
        # =====================================================================

        print(
            "\n" + "=" * 80
        )

        print(
            "📁 OUTPUT FILES"
        )

        print(
            "=" * 80
        )

        print(
            f"  {metadata_file}"
        )

        print(
            f"  {sequences_file}"
        )

        print(
            f"  {log_file}"
        )

        if all_failures:

            print(
                f"  {final_failed_file}"
            )

        print(
            "\n📊 FINAL SUMMARY"
        )

        print(
            f"   Total genomes: "
            f"{len(genome_records):,}"
        )

        print(
            f"   Total sequences: "
            f"{sequence_count:,}"
        )

        print(
            f"   Batch size: "
            f"{args.batch:,}"
        )

        print(
            f"   Parallel sequence requests: "
            f"{args.parallel_batches}"
        )

        print(
            f"   Retry attempts: "
            f"{args.retry_attempts}"
        )

        if os.path.exists(
            metadata_file
        ):

            print(
                f"   Metadata size: "
                f"{os.path.getsize(metadata_file) / (1024 * 1024):.2f} MB"
            )

        if os.path.exists(
            sequences_file
        ):

            print(
                f"   Sequence size: "
                f"{os.path.getsize(sequences_file) / (1024 * 1024):.2f} MB"
            )

        if all_failures:

            print(
                f"\n⚠️ Still failed: "
                f"{len(all_failures):,}"
            )

        else:

            print(
                "\n✅ All genomes/sequences "
                "downloaded successfully!"
            )

        print(
            "\n" + "=" * 80
        )

        print(
            "✅ DOWNLOAD COMPLETE"
        )

        print(
            "=" * 80
        )

    finally:

        # Always restore stdout, even if an exception propagated.
        sys.stdout = original_stdout

        print(
            f"\n📄 Log saved to: "
            f"{log_file}"
        )


if __name__ == "__main__":
    main()
