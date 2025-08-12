#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=99999:00:00
#PBS -l mem=100G
#PBS -j oe

set -euo pipefail

# --- Variables from qsub ---
TISSUE="${TISSUE:-}"
MERGED1="${MERGED1:-}" 
MERGED2="${MERGED2:-}" 
ISO1="${ISO1:-}"
ISO2="${ISO2:-}"

# --- Configuration ---
UTIL_PYTHON_SCRIPT="scripts/utils/gene_polya_matrix_v2.py"
UTIL_SHELL_SCRIPT="scripts/utils/concat_cpm_csv.sh"

OUTPUT_DIR_BASE="results/03_polya_matrix"
FINAL_OUTPUT_DIR="${OUTPUT_DIR_BASE}/${TISSUE}"

# --- Execution ---
cd "$PBS_O_WORKDIR"

mkdir -p "$FINAL_OUTPUT_DIR"

cp "$UTIL_SHELL_SCRIPT" "$FINAL_OUTPUT_DIR/"

python3 "$UTIL_PYTHON_SCRIPT" \
    --merged_data "$MERGED1" "$MERGED2" \
    --isoquant_path "$ISO1" "$ISO2" \
    --label "$TISSUE" \
    --output_dir "$FINAL_OUTPUT_DIR"