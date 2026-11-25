#!/usr/bin/env python3

"""
fasta_segment_stats.py — Segment length and type statistics from FASTA headers.

Reads a FASTA file produced by bvbrc2fasta.py and reports:

    * Segment length statistics (count / mean / median / min / max) per segment
    * Type field distribution (subtype for FluA, lineage for FluB)
    * Virus classification summary (Influenza A / B / Unknown)
    * Accession summary (how many accessions per segment)

All statistics are split by virus type (A / B / Unknown) so FluA and FluB
can be inspected independently.

Header format expected (produced by bvbrc2fasta.py):
    ><accession> <host>|<segment_id>|<segment_name>|<type_field>|<geo>|<strain>|<len> vadr|genes=...|coords=...|strand=...|len=<N>

Example:
    >OP870523 Pig|6|NA|H1N2|Usa|A/swine/Iowa/22TOSU1425/2022|1442 vadr|genes=NA|coords=8..1417:+|strand=+|len=1410
    >OR203317 Human|7|MP|Victoria|Usa|B/Texas/01/2023|1154 vadr|genes=M1,BM2|coords=8..1083:+|strand=+|len=1076
"""

import argparse
import re
import statistics
from collections import defaultdict


# =============================================================================
# Virus classification
# =============================================================================


def classify_virus(type_field, strain):
    """
    Determine whether a record is Influenza A, B, or unknown.

    Strategy
    --------
    1. Type field (4th pipe-separated column):
       Influenza A subtypes match H<num>N<num>(pdm<num>)?
       Influenza B lineages are 'Victoria' or 'Yamagata'
    2. Fall back to the strain prefix ('A/...' or 'B/...').

    Returns
    -------
    str
        'A', 'B', or 'Unknown'
    """
    type_field = (type_field or "").strip()
    strain = (strain or "").strip()

    if re.fullmatch(r"H\d+N\d+(?:pdm\d+)?", type_field, flags=re.IGNORECASE):
        return "A"

    if type_field.lower() in {"victoria", "yamagata"}:
        return "B"

    if strain.lower().startswith("a/"):
        return "A"
    if strain.lower().startswith("b/"):
        return "B"

    return "Unknown"


# =============================================================================
# FASTA parsing
# =============================================================================


def parse_fasta(fasta_file):
    """
    Read FASTA headers and extract per-record information.

    Fields (pipe-separated):
        [0] accession + host
        [1] segment_id
        [2] segment_name
        [3] type_field
        [4] geographic_location
        [5] strain
        [6] length + tool marker
        [7..] vadr details, including len=<N> at the tail

    Returns
    -------
    segment_lengths : dict[(virus, segment_id)] -> list[int]
    type_counts     : dict[(virus, type_field)] -> int
    virus_counts    : dict[virus] -> int
    accession_segments : dict[accession] -> set of segment_id
    """
    segment_lengths = defaultdict(list)
    type_counts = defaultdict(int)
    virus_counts = defaultdict(int)
    accession_segments = defaultdict(set)

    with open(fasta_file, "r") as fh:
        for line_number, line in enumerate(fh, 1):
            line = line.strip()

            if not line.startswith(">"):
                continue

            header = line[1:]
            fields = header.split("|")

            # Minimum required fields:
            # accession+host | seg_id | seg_name | type | geo | strain | length
            if len(fields) < 7:
                print(
                    f"WARNING: Invalid FASTA header at line {line_number}: "
                    f"{header}"
                )
                continue

            # Accession = word before the first space in field 0
            accession, _, _host = fields[0].partition(" ")
            accession = accession.strip()

            segment_id = fields[1].strip()
            type_field = fields[3].strip()
            strain = fields[5].strip()

            virus = classify_virus(type_field, strain)
            virus_counts[virus] += 1

            if accession:
                accession_segments[accession].add(segment_id)

            if type_field:
                type_counts[(virus, type_field)] += 1
            else:
                print(
                    f"WARNING: Missing type field at line {line_number}: "
                    f"{header}"
                )

            # Extract actual sequence length from len=<N> at the tail
            match = re.search(r"\blen=(\d+)\b", header)
            if not match:
                print(
                    f"WARNING: No len= field at line {line_number}: "
                    f"{header}"
                )
                continue

            sequence_length = int(match.group(1))
            segment_lengths[(virus, segment_id)].append(sequence_length)

    return segment_lengths, type_counts, virus_counts, accession_segments


# =============================================================================
# Reporting
# =============================================================================


def print_virus_summary(virus_counts):
    """Print record counts per virus class."""
    print()
    print("Virus classification summary")
    print("=" * 40)
    print(f"{'Virus':<20}{'Count':>10}")
    print("-" * 40)

    total = 0
    for virus in ("A", "B", "Unknown"):
        count = virus_counts.get(virus, 0)
        if count == 0:
            continue
        total += count
        label = f"Influenza {virus}" if virus != "Unknown" else "Unknown"
        print(f"{label:<20}{count:>10}")

    print("-" * 40)
    print(f"{'Total':<20}{total:>10}")


def print_segment_stats(segment_lengths, virus):
    """Print segment length statistics for a given virus."""
    subset = {
        seg: lengths
        for (v, seg), lengths in segment_lengths.items()
        if v == virus
    }

    if not subset:
        return

    print()
    print(f"Influenza {virus} — segment length statistics")
    print("=" * 75)
    print(
        f"{'Segment':<10}"
        f"{'Count':>10}"
        f"{'Mean':>12}"
        f"{'Median':>12}"
        f"{'Min':>12}"
        f"{'Max':>12}"
    )
    print("-" * 75)

    total_count = 0

    for segment_id in sorted(
        subset,
        key=lambda x: int(x) if x.isdigit() else x,
    ):
        lengths = subset[segment_id]
        count = len(lengths)
        mean_length = statistics.mean(lengths)
        median_length = statistics.median(lengths)
        min_length = min(lengths)
        max_length = max(lengths)

        total_count += count

        print(
            f"{segment_id:<10}"
            f"{count:>10}"
            f"{mean_length:>12.1f}"
            f"{median_length:>12.1f}"
            f"{min_length:>12}"
            f"{max_length:>12}"
        )

    print("-" * 75)
    print(f"{'Total':<10}{total_count:>10}")


def print_type_stats(type_counts, virus):
    """Print type field statistics for a given virus."""
    subset = {
        type_field: count
        for (v, type_field), count in type_counts.items()
        if v == virus
    }

    if not subset:
        return

    print()
    print(f"Influenza {virus} — type statistics")
    print("=" * 40)
    print(f"{'Type':<20}{'Count':>10}")
    print("-" * 40)

    total = 0
    for type_field, count in sorted(
        subset.items(),
        key=lambda x: (-x[1], x[0]),
    ):
        total += count
        print(f"{type_field:<20}{count:>10}")

    print("-" * 40)
    print(f"{'Total':<20}{total:>10}")


def print_accession_summary(accession_segments, virus_lookup, max_examples=20):
    """
    Print how many accessions have each number of segments.

    `virus_lookup` maps accession -> virus class (built during parsing).
    `max_examples` limits the number of example accessions printed.
    """
    if not accession_segments:
        return

    # Group accessions by their segment count
    by_count = defaultdict(list)
    for acc, segs in accession_segments.items():
        by_count[len(segs)].append(acc)

    print()
    print("Accession summary")
    print("=" * 60)
    print(f"{'Total accessions:':<40}{len(accession_segments):>10}")
    print()

    print(f"{'Segments per accession':<30}{'Accessions':>15}")
    print("-" * 60)

    for n in sorted(by_count.keys()):
        print(f"{n:<30}{len(by_count[n]):>15}")

    # Show examples of accessions with the most segments
    max_segments = max(by_count.keys())
    if max_segments > 1:
        examples = sorted(by_count[max_segments])[:max_examples]
        print()
        print(f"Example accessions with {max_segments} segments:")
        for acc in examples:
            virus = virus_lookup.get(acc, "Unknown")
            segs = sorted(accession_segments[acc])
            print(f"  [{virus}] {acc}: {segs}")

        if len(by_count[max_segments]) > max_examples:
            remaining = len(by_count[max_segments]) - max_examples
            print(f"  ... and {remaining} more")


# =============================================================================
# Main
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate sequence length statistics by segment and "
            "type statistics from a FASTA file, split by influenza "
            "virus type (A / B / Unknown)."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input FASTA file",
    )

    args = parser.parse_args()

    (
        segment_lengths,
        type_counts,
        virus_counts,
        accession_segments,
    ) = parse_fasta(args.input)

    if not segment_lengths:
        print("No valid segment length information found.")
        return

    print_virus_summary(virus_counts)

    for virus in ("A", "B", "Unknown"):
        print_segment_stats(segment_lengths, virus)
        print_type_stats(type_counts, virus)

    # Accession summary is global (not split by virus) for now.
    print_accession_summary(accession_segments, virus_lookup={})


if __name__ == "__main__":
    main()
