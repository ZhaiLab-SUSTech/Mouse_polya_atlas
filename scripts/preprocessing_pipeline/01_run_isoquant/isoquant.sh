#!/bin/bash

set -euo pipefail

# --- Variables from qsub ---
BAM_FILE="${bam:-}"
SAMPLE_ID="${p1:-}"

# --- Configuration ---
REF_GENOME="data/annotation/mm10.noScaffold.fa"
GENE_DB="data/annotation/gencode.vM25.annotation.gtf.gz"


OUTPUT_DIR_BASE="results/01_isoquant_raw"
FINAL_OUTPUT_DIR="${OUTPUT_DIR_BASE}/${SAMPLE_ID}"

# --- Execution ---
cd "$PBS_O_WORKDIR"

mkdir -p "$FINAL_OUTPUT_DIR"

isoquant.py \
    -d nanopore \
    --bam "${BAM_FILE}" \
    --reference "${REF_GENOME}" \
    --genedb "${GENE_DB}" \
    --complete_genedb \
    --output "${FINAL_OUTPUT_DIR}" \
    --prefix "${SAMPLE_ID}" \
    --labels "${SAMPLE_ID}" \
    --force \
    --report_novel_unspliced true \
    --check_canonical \
    --sqanti_output