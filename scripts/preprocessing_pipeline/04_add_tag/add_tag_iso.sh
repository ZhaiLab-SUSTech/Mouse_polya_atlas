#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=10G
#PBS -j oe
#PBS -l walltime=200:00:00

set -euo pipefail

SAMPLE_NAME="${SAMPLE_NAME:-}"
INPUT_BAM="${INPUT_BAM:-}"
READ_INFO_PARQUET="${READ_INFO_PARQUET:-}"

# --- Configuration ---
ADD_TAG_SCRIPT="scripts/utils/add_tag_to_bam_iso.py"
GET_POLYA_SCRIPT="scripts/utils/get_polyadenylated_reads.py"
GET_ELONG_SCRIPT="scripts/utils/get_elongating_reads.py"

ALIGNED_DATA_DIR="data/aligned_data"

OUTPUT_DIR_BASE="results/04_tagged_bams"
FINAL_OUTPUT_DIR="${OUTPUT_DIR_BASE}/${SAMPLE_NAME}"

TAGGED_BAM="${FINAL_OUTPUT_DIR}/${SAMPLE_NAME}.tagged.bam"
POLYA_BAM="${FINAL_OUTPUT_DIR}/${SAMPLE_NAME}.polyadenylated.bam"
ELONG_BAM="${FINAL_OUTPUT_DIR}/${SAMPLE_NAME}.elongating.bam"

cd "$PBS_O_WORKDIR"
mkdir -p "$FINAL_OUTPUT_DIR"

python3 "$ADD_TAG_SCRIPT" \
    -i "$INPUT_BAM" \
    --read_info "$READ_INFO_PARQUET" \
    --adapter_info "${ALIGNED_DATA_DIR}/${SAMPLE_NAME}.adapter.result.txt" \
    --polya_info "${ALIGNED_DATA_DIR}/${SAMPLE_NAME}.polyA_tail.result.txt" \
    -o "$TAGGED_BAM"
samtools index "$TAGGED_BAM"

python3 "$GET_POLYA_SCRIPT" -i "$TAGGED_BAM" -o "$POLYA_BAM"
samtools index "$POLYA_BAM"

python3 "$GET_ELONG_SCRIPT" -i "$TAGGED_BAM" -o "$ELONG_BAM"
samtools index "$ELONG_BAM"