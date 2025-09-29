#!/bin/bash

set -euo pipefail

# --- Configuration ---
MAX_JOBS=8
POLLING_INTERVAL=100
JOB_NAME_PREFIX="polya_matrix_"
JOB_SCRIPT="scripts/preprocessing_pipeline/03_calculate_gene_polya_median/run_polya_matrix.sh"

ISOQUANT_DIR="results/01_isoquant_raw"
MERGED_DIR="results/02_merged_readinfo"

TISSUE_LIST=(
    # brain
    # thyroid
    # thymus
    # heart
    # lung
    # liver
    # spleen
    # pancreas
    # stomach
    # small
    # large
    # adrenal
    # kidney
    # muscle
    # adipose
    # bone
    # testis
    # sperm
    test
)


for tissue in "${TISSUE_LIST[@]}"; do
    while true; do
        running_jobs=$(qstat -u "$USER" | grep "$JOB_NAME_PREFIX" | grep -c '[QR]' || true)

        if [[ $running_jobs -lt $MAX_JOBS ]]; then
            merged1="${MERGED_DIR}/${tissue}/${tissue}_rep1.parquet"
            merged2="${MERGED_DIR}/${tissue}/${tissue}_rep2.parquet"
            iso1="${ISOQUANT_DIR}/${tissue}_rep1/${tissue}_rep1"
            iso2="${ISOQUANT_DIR}/${tissue}_rep2/${tissue}_rep2"

            echo "Slot available. Submitting job for ${tissue}."

            qsub -N "${JOB_NAME_PREFIX}${tissue}" \
                 -v TISSUE="$tissue",MERGED1="$merged1",MERGED2="$merged2",ISO1="$iso1",ISO2="$iso2" \
                 "$JOB_SCRIPT"

            sleep 1
            break
        else
            echo "Queue is full (${running_jobs}/${MAX_JOBS}). Waiting..."
            sleep "$POLLING_INTERVAL"
        fi
    done
done