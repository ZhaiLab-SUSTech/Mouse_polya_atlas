#!/bin/bash

set -euo pipefail

# --- Configuration ---
MAX_JOBS=5
POLLING_INTERVAL=100
JOB_NAME_PREFIX="tissue_job_"
JOB_SCRIPT="scripts/preprocessing_pipeline/02_merge_isoquant2polya/run_tissue_job.sh"

# Define the list of tissues to be processed.
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
    # small # short name for small intestine
    # large # short name for large intestine
    # adrenal # short name for adrenal gland
    # kidney
    # muscle
    # adipose
    # bone # short name for bonemarrow
    # testis
    # sperm
    test
)

for tissue in "${TISSUE_LIST[@]}"; do
    while true; do
        running_jobs=$(qstat -u "$USER" | grep "$JOB_NAME_PREFIX" | grep -c '[QR]' || true)

        if [[ $running_jobs -lt $MAX_JOBS ]]; then
            rep1="${tissue}_rep1"
            rep2="${tissue}_rep2"

            echo "Slot available. Submitting job for ${tissue}."

            qsub -N "${JOB_NAME_PREFIX}${tissue}" \
                 -v TISSUE_NAME="$tissue",REP1="$rep1",REP2="$rep2" \
                 "$JOB_SCRIPT"

            sleep 1 # Stagger submissions slightly.
            break
        else
            echo "Queue is full (${running_jobs}/${MAX_JOBS}). Waiting..."
            sleep "$POLLING_INTERVAL"
        fi
    done
done