#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=99999:00:00
#PBS -l mem=100G
#PBS -j oe

set -euo pipefail

# --- Variables from qsub ---
TISSUE_NAME="${TISSUE_NAME:-}"
REP1="${REP1:-}"
REP2="${REP2:-}"

# --- Configuration ---
UTIL_SCRIPT="scripts/utils/combine_isoform_v3.py"
NUM_THREADS=4


ISOQUANT_INPUT_DIR="results/01_isoquant_raw"
GTF_INPUT_DIR="data/annotation/gencode.vM25.annotation.gtf.gz"

MERGED_READINFO_DIR="results/02_merged_readinfo"
REPLICATES_RESULTS_DIR="data/polya_results"

# --- Execution ---
cd "$PBS_O_WORKDIR"

mkdir -p "${MERGED_READINFO_DIR}"
mkdir -p "${REPLICATES_RESULTS_DIR}/${TISSUE_NAME}"

python3 "$UTIL_SCRIPT" \
    --replicates "$REP1" "$REP2" \
    --isoquant_dir "$ISOQUANT_INPUT_DIR" \
    --gtf_dir "${GTF_INPUT_DIR}" \
    --output_dir "${MERGED_READINFO_DIR}" \
    --results_dir "${REPLICATES_RESULTS_DIR}" \
    --num_processes "$NUM_THREADS"