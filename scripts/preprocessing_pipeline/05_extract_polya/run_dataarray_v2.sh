#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=99999:00:00
#PBS -l mem=100G
#PBS -j oe

set -euo pipefail

TISSUE="${TISSUE:-}"
REP1="${REP1:-}"
REP2="${REP2:-}"
GTF1="${GTF1:-}"
GTF2="${GTF2:-}"

UTIL_SCRIPT="scripts/utils/process_dataarray_v2.py"

NUM_THREADS=2 
BIN_WIDTH=5

OUTPUT_DIR_BASE="results/05_extracted_polya"
FINAL_OUTPUT_DIR="${OUTPUT_DIR_BASE}/${TISSUE}"

cd "$PBS_O_WORKDIR"
mkdir -p "$FINAL_OUTPUT_DIR"

python3 "$UTIL_SCRIPT" \
    --merged_data "$REP1" "$REP2" \
    --gtf_path "$GTF1" "$GTF2" \
    --label "$TISSUE" \
    --output_dir "$FINAL_OUTPUT_DIR" \
    --bin_width "$BIN_WIDTH" \
    --thread "$NUM_THREADS"