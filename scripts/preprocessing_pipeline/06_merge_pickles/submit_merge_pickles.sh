#!/bin/bash
#PBS -N merge_pickles
#PBS -l nodes=1:ppn=4,mem=50G
#PBS -j oe
#PBS -l walltime=24:00:00

set -euo pipefail
cd "$PBS_O_WORKDIR"

# --- Configuration ---
SCRIPT_PATH="scripts/preprocessing_pipeline/06_merge_data/merge_pickles.py"

INPUT_DIR="results/05_extracted_polya"

OUTPUT_INFO="results/06_merged_data/merged_info_df.pkl"
OUTPUT_ARRAY="results/06_merged_data/merged_array.npy"

# --- Execution ---
python3 "$SCRIPT_PATH" \
    --input_dir "$INPUT_DIR" \
    --output_info "$OUTPUT_INFO" \
    --output_array "$OUTPUT_ARRAY"
