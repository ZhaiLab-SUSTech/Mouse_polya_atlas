#!/bin/bash

set -euo pipefail

# --- Configuration ---
MAX_JOBS=5
POLLING_INTERVAL=100
JOB_NAME_PREFIX="extract_polya_"
JOB_SCRIPT="scripts/preprocessing_pipeline/05_extract_polya/run_extract_polya.sh"

ISOQUANT_DIR="results/01_isoquant_raw"
MERGED_DIR="results/02_merged_readinfo"

TISSUE_LIST=(
    brain
    thyroid
    thymus
    heart
    lung
    liver
    spleen
    pancreas
    stomach
    small # short name for small intestine
    large # short name for large intestine
    adrenal
    kidney
    muscle
    adipose
    bone
    testis
    sperm
)

for tissue in "${TISSUE_LIST[@]}"; do
    while true; do
        running_jobs=$(qstat -u "$USER" | grep "$JOB_NAME_PREFIX" | grep -c '[QR]' || true)

        if [[ $running_jobs -lt $MAX_JOBS ]]; then
            rep1="${MERGED_DIR}/${tissue}/${tissue}_rep1.parquet"
            rep2="${MERGED_DIR}/${tissue}/${tissue}_rep2.parquet"
            gtf1="${ISOQUANT_DIR}/${tissue}_rep1/${tissue}_rep1.extended_annotation.gtf"
            gtf2="${ISOQUANT_DIR}/${tissue}_rep2/${tissue}_rep2.extended_annotation.gtf"

            echo "Slot available. Submitting job for ${tissue}."

            qsub -N "${JOB_NAME_PREFIX}${tissue}" \
                 -v TISSUE="$tissue",REP1="$rep1",REP2="$rep2",GTF1="$gtf1",GTF2="$gtf2" \
                 "$JOB_SCRIPT"

            sleep 1
            break
        else
            echo "Queue is full (${running_jobs}/${MAX_JOBS}). Waiting..."
            sleep "$POLLING_INTERVAL"
        fi
    done
done