#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

VENV_DIR="${PIPELINE_DIR}/.venv"
REQUIREMENTS="${PIPELINE_DIR}/requirements.txt"

PIPELINE_NF="${PIPELINE_DIR}/main.nf"
DOWNLOADER="${PIPELINE_DIR}/bin/bvbrc_downloader.py"

# -----------------------------------------------------------------------------
# Defaults for CLI parameters
# -----------------------------------------------------------------------------

BATCH=200
SEQ_LIMIT=200
MAX_RECORDS=""

# -----------------------------------------------------------------------------
# Help function
# -----------------------------------------------------------------------------

show_help() {
    cat << EOF
Usage: bash run_pipeline.sh [OPTIONS] [NEXTFLOW_OPTIONS]

Options:
  -h, --help              Show this help message and exit.
  --outdir <path>         Output directory (default: current directory).
  --batch <n>             Metadata page size / sequence batch size (default: 200).
  --seq-limit <n>         Server-side sequence cap per request (default: 200).
  --max <n>               Max records per taxon (default: no limit).

Environment Variables:
  SKIP_DOWNLOAD=1         Skip the BV-BRC data download step.
  DB_VERSION=<string>     Override the database version tag (default: YYYYMMDD).

Nextflow Options:
  Any unrecognized arguments are forwarded directly to `nextflow run`.
  Examples:
    -profile <name>       Use specific Nextflow profiles (e.g., -profile slurm)
    -c <file>             Use a custom Nextflow config file
    -resume               Resume a previous run
    -w <path>             Set a custom work directory

Examples:
  bash run_pipeline.sh --outdir results
  bash run_pipeline.sh --outdir results --max 500 -profile singularity,slurm
  SKIP_DOWNLOAD=1 bash run_pipeline.sh --outdir results -resume
EOF
}

# -----------------------------------------------------------------------------
# Parse CLI arguments
# -----------------------------------------------------------------------------

OUTDIR=""
POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            show_help
            exit 0
            ;;
        --outdir)
            [[ -z "${2:-}" ]] && { echo "ERROR: --outdir requires a value" >&2; exit 1; }
            OUTDIR="$2"
            shift 2
            ;;
        --outdir=*)
            OUTDIR="${1#*=}"
            shift
            ;;
        --batch)
            [[ -z "${2:-}" ]] && { echo "ERROR: --batch requires a value" >&2; exit 1; }
            BATCH="$2"
            shift 2
            ;;
        --batch=*)
            BATCH="${1#*=}"
            shift
            ;;
        --seq-limit)
            [[ -z "${2:-}" ]] && { echo "ERROR: --seq-limit requires a value" >&2; exit 1; }
            SEQ_LIMIT="$2"
            shift 2
            ;;
        --seq-limit=*)
            SEQ_LIMIT="${1#*=}"
            shift
            ;;
        --max)
            [[ -z "${2:-}" ]] && { echo "ERROR: --max requires a value" >&2; exit 1; }
            MAX_RECORDS="$2"
            shift 2
            ;;
        --max=*)
            MAX_RECORDS="${1#*=}"
            shift
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

# -----------------------------------------------------------------------------
# Validate numeric parameters
# -----------------------------------------------------------------------------

for var_name in BATCH SEQ_LIMIT; do
    value="${!var_name}"
    if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: --${var_name,,} must be a positive integer (got: ${value})" >&2
        exit 1
    fi
    if [[ "${value}" -lt 1 ]]; then
        echo "ERROR: --${var_name,,} must be >= 1 (got: ${value})" >&2
        exit 1
    fi
done

if [[ -n "${MAX_RECORDS}" ]]; then
    if ! [[ "${MAX_RECORDS}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: --max must be a positive integer (got: ${MAX_RECORDS})" >&2
        exit 1
    fi
    if [[ "${MAX_RECORDS}" -lt 1 ]]; then
        echo "ERROR: --max must be >= 1 (got: ${MAX_RECORDS})" >&2
        exit 1
    fi
fi

# -----------------------------------------------------------------------------
# Resolve OUTDIR
# -----------------------------------------------------------------------------

if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="${OUTDIR_ENV:-$(pwd)}"
fi

OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
DATA_DIR="${OUTDIR}/data"

# -----------------------------------------------------------------------------
# Other configuration
# -----------------------------------------------------------------------------

TAXON="11320,11520"
PREFIX="bvbrc-influenza_ab"
SKIP_DOWNLOAD=${SKIP_DOWNLOAD:-0}

export PIP_CACHE_DIR="${PIPELINE_DIR}/.pip-cache"

# -----------------------------------------------------------------------------
# Step 0 — Ensure virtualenv exists
# -----------------------------------------------------------------------------

ensure_venv() {
    if [ -x "${VENV_DIR}/bin/python" ]; then
        echo "=== Using existing virtualenv: ${VENV_DIR} ==="
        return
    fi

    echo "=== Creating virtualenv: ${VENV_DIR} ==="
    python3 -m venv "${VENV_DIR}"
    "${VENV_DIR}/bin/pip" install --upgrade pip

    if [ -f "${REQUIREMENTS}" ]; then
        echo "=== Installing dependencies from ${REQUIREMENTS} ==="
        "${VENV_DIR}/bin/pip" install -r "${REQUIREMENTS}"
    else
        echo "=== No requirements.txt; installing bvbrc directly ==="
        "${VENV_DIR}/bin/pip" install bvbrc
    fi
}

# -----------------------------------------------------------------------------
# Step 1 — Download metadata + sequences
# -----------------------------------------------------------------------------

download_data() {
    echo "=== Downloading from BV-BRC ==="
    echo "    Target:    ${DATA_DIR}"
    echo "    Batch:     ${BATCH}"
    echo "    Seq limit: ${SEQ_LIMIT}"
    if [[ -n "${MAX_RECORDS}" ]]; then
        echo "    Max:       ${MAX_RECORDS} per taxon"
    else
        echo "    Max:       (no limit)"
    fi

    if [ ! -f "${DOWNLOADER}" ]; then
        echo "ERROR: downloader not found: ${DOWNLOADER}" >&2
        exit 1
    fi

    mkdir -p "${DATA_DIR}"

    local original_dir
    original_dir="$(pwd)"
    cd "${DATA_DIR}"

    # shellcheck disable=SC1091
    source "${VENV_DIR}/bin/activate"

    local max_arg=()
    if [[ -n "${MAX_RECORDS}" ]]; then
        max_arg=(--max "${MAX_RECORDS}")
    fi

    python "${DOWNLOADER}" \
        --taxon "${TAXON}" \
        --batch "${BATCH}" \
        --seq-limit "${SEQ_LIMIT}" \
        "${max_arg[@]}" \
        --prefix "${PREFIX}"

    deactivate

    META_SRC=$(ls -t "${PREFIX}"_*_meta_raw.json 2>/dev/null | head -1 || true)
    SEQ_SRC=$(ls -t "${PREFIX}"_*_sequences_raw.json 2>/dev/null | head -1 || true)

    if [ -z "${META_SRC}" ] || [ -z "${SEQ_SRC}" ]; then
        echo "ERROR: download did not produce expected JSON files" >&2
        echo "       Looked in: $(pwd)" >&2
        cd "${original_dir}"
        exit 1
    fi

    mv "${META_SRC}" "meta_raw.json"
    mv "${SEQ_SRC}" "sequences_raw.json"

    cd "${original_dir}"

    echo "Downloaded:"
    echo "  ${DATA_DIR}/meta_raw.json"
    echo "  ${DATA_DIR}/sequences_raw.json"
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

echo "=== Configuration ==="
echo "    Pipeline root : ${PIPELINE_DIR}"
echo "    Output dir    : ${OUTDIR}"
echo "    Data dir      : ${DATA_DIR}"
echo "    Batch         : ${BATCH}"
echo "    Seq limit     : ${SEQ_LIMIT}"
if [[ -n "${MAX_RECORDS}" ]]; then
    echo "    Max records   : ${MAX_RECORDS}"
else
    echo "    Max records   : (no limit)"
fi
echo ""

ensure_venv

if [ "${SKIP_DOWNLOAD}" != "1" ]; then
    download_data
else
    echo "=== Step 1: SKIPPED (SKIP_DOWNLOAD=1) ==="
fi

echo ""
echo "=== Running Nextflow pipeline ==="

TODAY="$(date +%Y%m%d)"
DB_VERSION="${DB_VERSION:-${TODAY}}"
echo "    db_version: ${DB_VERSION}"

eval "$(conda shell.bash hook)"

set +u
conda activate nf-core-4
set -u

nextflow run "${PIPELINE_NF}" -profile singularity -resume \
    --db_version "${DB_VERSION}" \
    --outdir "${OUTDIR}" \
    --meta_json_file "${DATA_DIR}/meta_raw.json" \
    --seq_json_file "${DATA_DIR}/sequences_raw.json" \
    "${POSITIONAL_ARGS[@]}"

conda deactivate
