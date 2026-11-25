#!/usr/bin/env python

import argparse
import re
import sys
import logging
import traceback
from collections import defaultdict


def parse_cd_hit_clstr(clstr_file):
    """
    Parse a CD-HIT .clstr file.

    Returns a list of dictionaries:

        [
            {
                "cluster_id": ">Cluster 0",
                "sequences": [...]
            },
            ...
        ]
    """
    clusters = []
    current_cluster = None

    try:
        with open(clstr_file, "r") as file:
            for line in file:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">Cluster"):
                    if current_cluster is not None:
                        clusters.append(current_cluster)

                    current_cluster = {
                        "cluster_id": line,
                        "sequences": []
                    }

                else:
                    if current_cluster is not None:
                        current_cluster["sequences"].append(line)

        if current_cluster is not None:
            clusters.append(current_cluster)

    except FileNotFoundError:
        logging.error(
            f"The file {clstr_file} does not exist."
        )
        sys.exit(1)

    except Exception as e:
        logging.error(
            f"Error reading CLSTR file {clstr_file}: {e}"
        )
        logging.error(traceback.format_exc())
        sys.exit(1)

    return clusters


def parse_fasta_header(header):
    """
    Parse a FASTA header.

    Expected format:

        >SEQ_ID DESCRIPTION

    Example:

        >PZ610282|11320.1516667 Gull|2|PB1|H5N1|Germany|A/gull/Germany-SN/VD202605383/2026|2306

    The split occurs only at the first whitespace.

    Returns:

        seq_id:
            PZ610282|11320.1516667

        description:
            Gull|2|PB1|H5N1|Germany|A/gull/Germany-SN/VD202605383/2026|2306
    """
    header = header.lstrip(">")

    parts = header.split(None, 1)

    seq_id = parts[0] if parts else ""
    description = parts[1] if len(parts) > 1 else ""

    return seq_id, description


def parse_fasta_description(description):
    """
    Split the description into pipe-separated fields.

    Fields:

        0 = host
        1 = segment
        2 = gene
        3 = subtype
        4 = location
        5 = strain
        6 = length
    """
    if not description:
        return []

    return description.split("|")


def extract_subtype_from_fasta(fasta_file):
    """
    Build a mapping:

        sequence_id -> subtype
    """
    id2subtype = {}

    try:
        with open(fasta_file, "r") as file:
            for line_number, line in enumerate(file, start=1):

                line = line.rstrip("\r\n")

                if not line.startswith(">"):
                    continue

                seq_id, description = parse_fasta_header(line)
                fields = parse_fasta_description(description)

                # Subtype is the fourth field:
                #
                # Gull | 2 | PB1 | H5N1 | Germany | strain | length
                #   0     1     2      3

                if len(fields) <= 3:
                    logging.warning(
                        f"FASTA header at line {line_number} has fewer "
                        f"than 4 description fields: {line}"
                    )
                    continue

                subtype = fields[3].strip()

                if not subtype:
                    logging.warning(
                        f"Empty subtype for sequence {seq_id} "
                        f"at line {line_number}"
                    )
                    continue

                id2subtype[seq_id] = subtype

    except FileNotFoundError:
        logging.error(
            f"The FASTA file {fasta_file} does not exist."
        )
        sys.exit(1)

    except Exception as e:
        logging.error(
            f"Error reading FASTA file {fasta_file}: {e}"
        )
        logging.error(traceback.format_exc())
        sys.exit(1)

    logging.info(
        f"Loaded subtype information for "
        f"{len(id2subtype)} sequences"
    )

    return id2subtype


def extract_fasta_headers(fasta_file):
    """
    Build a mapping:

        sequence_id -> complete original FASTA header

    The complete header is retained exactly as it appears in the
    input FASTA, including the leading '>'.

    This is used to report complete FASTA information for excluded
    singleton clusters.
    """
    id2header = {}

    try:
        with open(fasta_file, "r") as file:
            for line_number, line in enumerate(file, start=1):

                line = line.rstrip("\r\n")

                if not line.startswith(">"):
                    continue

                seq_id, _ = parse_fasta_header(line)

                if not seq_id:
                    logging.warning(
                        f"Empty sequence ID in FASTA header "
                        f"at line {line_number}: {line}"
                    )
                    continue

                id2header[seq_id] = line

    except FileNotFoundError:
        logging.error(
            f"The FASTA file {fasta_file} does not exist."
        )
        sys.exit(1)

    except Exception as e:
        logging.error(
            f"Error reading FASTA headers from "
            f"{fasta_file}: {e}"
        )
        logging.error(traceback.format_exc())
        sys.exit(1)

    logging.info(
        f"Loaded complete FASTA headers for "
        f"{len(id2header)} sequences"
    )

    return id2header


def extract_sequence_details(seq):
    """
    Extract sequence length and sequence ID from a CD-HIT cluster line.

    Typical CD-HIT lines:

        0       2233nt, >NC_007371.1... *
        1       2233nt, >NC_007376.1... at 1:2233:1:2233/+/93.01%

    The sequence ID may contain '|', therefore the ID is captured
    until the CD-HIT '...' marker.

    Returns:

        (sequence_length, seqid)

    or:

        None
    """

    pattern = r"(\d+)nt,\s*>(.+?)\.\.\."

    match = re.search(pattern, seq)

    if not match:
        logging.warning(
            f"Could not parse CD-HIT sequence line: {seq}"
        )
        return None

    sequence_length = int(match.group(1))
    seqid = match.group(2)

    return sequence_length, seqid


def get_representative_sequence(
    cluster,
    id2subtype,
    id2header,
    exclusion_counts
):
    """
    Select the representative sequence for a CD-HIT cluster.

    Selection rules:

    1. Singleton clusters are excluded.
    2. Determine subtype composition of the cluster.
    3. Identify the dominant subtype.
    4. The dominant subtype must represent at least 70%.
    5. If the dominant subtype is <70%, reject the entire cluster.
    6. If the dominant subtype is >=70%, select the longest
       sequence belonging to that dominant subtype.
    7. There is NO fallback to the CD-HIT default representative.

    The input FASTA has already been QC-filtered to remove sequences
    containing ambiguous bases, so ambiguous-base count is not used
    as a representative-selection criterion.
    """

    cluster_id = cluster["cluster_id"]
    cluster_size = len(cluster["sequences"])

    # ---------------------------------------------------------------
    # Exclude singleton clusters
    # ---------------------------------------------------------------
    if cluster_size == 1:

        exclusion_counts["singleton"] += 1

        seq_details = extract_sequence_details(
            cluster["sequences"][0]
        )

        logging.info(
            f"EXCLUDED {cluster_id}: "
            f"reason=singleton"
        )

        if seq_details is not None:

            sequence_length, seqid = seq_details

            # Get the complete original FASTA header.
            full_header = id2header.get(
                seqid,
                f">{seqid}"
            )

            # Get subtype if available.
            subtype = id2subtype.get(
                seqid,
                "UNKNOWN"
            )

            logging.info(
                f"  FASTA header: {full_header}"
            )

            logging.info(
                f"  subtype={subtype}, "
                f"length={sequence_length}"
            )

        else:

            logging.info(
                f"  Could not parse sequence information "
                f"for singleton cluster."
            )

        return None

    # ---------------------------------------------------------------
    # Process sequences in the cluster
    # ---------------------------------------------------------------
    seq_with_subtype = []
    subtype_count = defaultdict(int)

    for seq in cluster["sequences"]:

        seq_details = extract_sequence_details(seq)

        if seq_details is None:
            logging.warning(
                f"{cluster_id}: Could not parse sequence line. "
                f"Sequence skipped."
            )
            continue

        sequence_length, seqid = seq_details

        subtype = id2subtype.get(seqid)

        if subtype is None:
            logging.warning(
                f"{cluster_id}: subtype not found for "
                f"sequence ID {seqid}. "
                f"Sequence skipped for subtype selection."
            )
            continue

        subtype_count[subtype] += 1

        seq_with_subtype.append(
            (
                sequence_length,
                seq,
                subtype
            )
        )

    # ---------------------------------------------------------------
    # No usable sequences
    # ---------------------------------------------------------------
    if not subtype_count:

        exclusion_counts["no_subtype"] += 1

        logging.info(
            f"EXCLUDED {cluster_id}: "
            f"reason=no_subtype_information, "
            f"cluster_size={cluster_size}"
        )

        return None

    total_count = sum(subtype_count.values())

    # ---------------------------------------------------------------
    # Determine dominant subtype
    # ---------------------------------------------------------------
    sorted_subtype_count = sorted(
        subtype_count.items(),
        key=lambda x: x[1],
        reverse=True
    )

    abundant_subtype, abundant_count = sorted_subtype_count[0]

    abundant_proportion = abundant_count / total_count

    # ---------------------------------------------------------------
    # Reject mixed clusters
    # ---------------------------------------------------------------
    if abundant_proportion < 0.70:

        exclusion_counts["mixed_lt_70"] += 1

        logging.info(
            f"EXCLUDED {cluster_id}: "
            f"reason=dominant_subtype_below_70_percent, "
            f"dominant_subtype={abundant_subtype}, "
            f"proportion={abundant_proportion:.2f}, "
            f"count={abundant_count}/{total_count}, "
            f"cluster_size={cluster_size}"
        )

        logging.info(
            f"  Subtype composition for {cluster_id}:"
        )

        for subtype, count in sorted_subtype_count:

            proportion = count / total_count

            logging.info(
                f"    subtype={subtype}, "
                f"count={count}/{total_count}, "
                f"proportion={proportion:.2f}"
            )

        return None

    # ---------------------------------------------------------------
    # Select sequences belonging to dominant subtype
    # ---------------------------------------------------------------
    dominant_subtype_sequences = [
        (
            sequence_length,
            seq,
            subtype
        )
        for sequence_length, seq, subtype in seq_with_subtype
        if subtype == abundant_subtype
    ]

    # ---------------------------------------------------------------
    # Safety check
    # ---------------------------------------------------------------
    if not dominant_subtype_sequences:

        exclusion_counts["no_dominant_sequence"] += 1

        logging.info(
            f"EXCLUDED {cluster_id}: "
            f"reason=no_sequence_for_dominant_subtype, "
            f"dominant_subtype={abundant_subtype}"
        )

        return None

    # ---------------------------------------------------------------
    # Select longest sequence
    # ---------------------------------------------------------------
    dominant_subtype_sequences.sort(
        key=lambda x: x[0],
        reverse=True
    )

    (
        representative_length,
        representative_seq,
        representative_subtype
    ) = dominant_subtype_sequences[0]

    representative_details = extract_sequence_details(
        representative_seq
    )

    representative_seqid = (
        representative_details[1]
        if representative_details is not None
        else "UNKNOWN"
    )

    logging.info(
        f"SELECTED {cluster_id}: "
        f"dominant_subtype={abundant_subtype}, "
        f"proportion={abundant_proportion:.2f}, "
        f"count={abundant_count}/{total_count}, "
        f"sequence={representative_seqid}, "
        f"length={representative_length}"
    )

    return representative_seq


def extract_fasta_subset(
    input_file,
    seqids,
    output_file=None
):
    """
    Extract sequences from a FASTA file based on sequence IDs.

    The sequence ID is everything before the first whitespace.

    The original description is retained unchanged.
    """

    extracted_sequences = []

    seqids = set(seqids)

    current_seqid = None
    current_sequence_desc = ""
    current_sequence = []

    try:
        with open(input_file, "r") as file:

            for line in file:
                line = line.rstrip("\r\n")

                if line.startswith(">"):

                    # Save previous sequence
                    if (
                        current_seqid is not None
                        and current_seqid in seqids
                    ):
                        extracted_sequences.append(
                            f">{current_seqid}"
                            f" {current_sequence_desc}\n"
                            f"{''.join(current_sequence)}"
                        )

                    # Parse new header
                    (
                        current_seqid,
                        current_sequence_desc
                    ) = parse_fasta_header(line)

                    current_sequence = []

                else:
                    current_sequence.append(
                        line.strip()
                    )

            # Handle final sequence
            if (
                current_seqid is not None
                and current_seqid in seqids
            ):
                extracted_sequences.append(
                    f">{current_seqid}"
                    f" {current_sequence_desc}\n"
                    f"{''.join(current_sequence)}"
                )

    except FileNotFoundError:
        logging.error(
            f"The FASTA file {input_file} does not exist."
        )
        sys.exit(1)

    except Exception as e:
        logging.error(
            f"Error reading input FASTA file: {e}"
        )
        logging.error(traceback.format_exc())
        sys.exit(1)

    if output_file:

        try:
            with open(output_file, "w") as out_file:

                if extracted_sequences:
                    out_file.write(
                        "\n".join(extracted_sequences)
                    )
                    out_file.write("\n")

        except Exception as e:
            logging.error(
                f"Error writing output FASTA "
                f"{output_file}: {e}"
            )
            logging.error(traceback.format_exc())
            sys.exit(1)

    else:
        return "\n".join(extracted_sequences)


def parse_clusters_and_select_representatives(
    clusters,
    id2subtype,
    id2header,
    fasta_file,
    output_file
):
    """
    Select representatives from all CD-HIT clusters and extract
    them from the original FASTA file.

    Also keeps track of why clusters were excluded.
    """

    representative_sequences = []
    representative_seqids = []

    exclusion_counts = defaultdict(int)

    for cluster in clusters:

        rep_seq = get_representative_sequence(
            cluster,
            id2subtype,
            id2header,
            exclusion_counts
        )

        if rep_seq is None:
            continue

        seq_details = extract_sequence_details(
            rep_seq
        )

        if seq_details is None:

            exclusion_counts[
                "representative_parse_error"
            ] += 1

            logging.info(
                f"EXCLUDED {cluster['cluster_id']}: "
                f"reason=could_not_parse_selected_representative"
            )

            continue

        _, seqid = seq_details

        representative_sequences.append(
            {
                "cluster_id": cluster["cluster_id"],
                "seq": rep_seq
            }
        )

        representative_seqids.append(seqid)

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    selected_count = len(
        representative_seqids
    )

    excluded_count = sum(
        exclusion_counts.values()
    )

    logging.info("")
    logging.info("=" * 70)
    logging.info(
        "REPRESENTATIVE SELECTION SUMMARY"
    )
    logging.info("=" * 70)

    logging.info(
        f"Total CD-HIT clusters: {len(clusters)}"
    )

    logging.info(
        f"Selected representatives: {selected_count}"
    )

    logging.info(
        f"Excluded clusters: {excluded_count}"
    )

    if exclusion_counts:

        logging.info("")
        logging.info(
            "Exclusion reasons:"
        )

        for reason, count in sorted(
            exclusion_counts.items()
        ):
            logging.info(
                f"  {reason}: {count}"
            )

    logging.info("=" * 70)

    # ---------------------------------------------------------------
    # Extract selected sequences from FASTA
    # ---------------------------------------------------------------
    extract_fasta_subset(
        fasta_file,
        representative_seqids,
        output_file
    )

    return representative_sequences


def write_output_file(
    output_file,
    representative_sequences
):
    """
    Write selected CD-HIT representative information.
    """

    try:
        with open(output_file, "w") as out_file:

            for rep in representative_sequences:

                out_file.write(
                    f"Cluster ID: {rep['cluster_id']}, "
                    f"Sequence: {rep['seq']}\n"
                )

        logging.info(
            f"Representative sequences saved to "
            f"{output_file}"
        )

    except Exception as e:
        logging.error(
            f"Error writing to file {output_file}: {e}"
        )
        logging.error(traceback.format_exc())
        sys.exit(1)


def print_representative_sequences(
    representative_sequences
):
    """
    Print representative sequences to stdout.
    """

    for rep in representative_sequences:

        print(
            f"Cluster ID: {rep['cluster_id']}, "
            f"Sequence: {rep['seq']}"
        )


def main():

    description = """
    Select representative sequences from CD-HIT clusters.

    Selection logic:

      1. Exclude singleton CD-HIT clusters.
      2. Determine subtype composition for each remaining cluster.
      3. Identify the dominant subtype.
      4. Require the dominant subtype to represent at least 70%
         of sequences with available subtype information.
      5. If the dominant subtype is <70%, reject the entire cluster.
      6. If the dominant subtype is >=70%, select the longest
         sequence belonging to that subtype.
      7. Do not fall back to the CD-HIT default representative.

    The input FASTA is assumed to have already passed QC for
    ambiguous bases.

    FASTA header format:

      >SEQ_ID DESCRIPTION

    Example:

      >PZ610282|11320.1516667 Gull|2|PB1|H5N1|Germany|A/gull/Germany-SN/VD202605383/2026|2306

    The sequence ID is:

      PZ610282|11320.1516667

    The description is:

      Gull|2|PB1|H5N1|Germany|A/gull/Germany-SN/VD202605383/2026|2306

    Description fields are separated by '|'.
    Subtype is field 4 (zero-based index 3).
    """

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-c",
        "--clstr",
        required=True,
        help="Input CLSTR file from cd-hit-est output."
    )

    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        help="Input FASTA file used to generate the cd-hit-est results."
    )

    parser.add_argument(
        "--out_clstr",
        required=True,
        help="Output file containing selected representatives."
    )

    parser.add_argument(
        "--out_fasta",
        required=True,
        help="Output FASTA containing selected representative sequences."
    )

    parser.add_argument(
        "--print",
        action="store_true",
        help="Print representative sequences to the console."
    )

    args = parser.parse_args()

    # ---------------------------------------------------------------
    # Logging
    # ---------------------------------------------------------------
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    logging.info(
        "Starting representative sequence selection"
    )

    logging.info(
        f"FASTA: {args.fasta}"
    )

    logging.info(
        f"CLSTR: {args.clstr}"
    )

    logging.info(
        "Selection rules: "
        "singleton clusters excluded; "
        "dominant subtype must be >=70%; "
        "longest sequence from dominant subtype selected."
    )

    # ---------------------------------------------------------------
    # Extract subtype information
    # ---------------------------------------------------------------
    id2subtype = extract_subtype_from_fasta(
        args.fasta
    )

    # ---------------------------------------------------------------
    # Extract complete FASTA headers
    # ---------------------------------------------------------------
    id2header = extract_fasta_headers(
        args.fasta
    )

    # ---------------------------------------------------------------
    # Parse CD-HIT clusters
    # ---------------------------------------------------------------
    clusters = parse_cd_hit_clstr(
        args.clstr
    )

    logging.info(
        f"Loaded {len(clusters)} CD-HIT clusters"
    )

    # ---------------------------------------------------------------
    # Select representatives
    # ---------------------------------------------------------------
    representative_sequences = (
        parse_clusters_and_select_representatives(
            clusters,
            id2subtype,
            id2header,
            args.fasta,
            args.out_fasta
        )
    )

    # ---------------------------------------------------------------
    # Write representative CLSTR information
    # ---------------------------------------------------------------
    write_output_file(
        args.out_clstr,
        representative_sequences
    )

    # ---------------------------------------------------------------
    # Optional console output
    # ---------------------------------------------------------------
    if args.print:

        print_representative_sequences(
            representative_sequences
        )

    logging.info(
        f"Finished. Selected "
        f"{len(representative_sequences)} representatives."
    )


if __name__ == "__main__":
    main()
