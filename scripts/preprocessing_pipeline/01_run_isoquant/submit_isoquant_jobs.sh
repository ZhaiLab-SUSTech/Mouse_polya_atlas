#!/bin/bash

set -euo pipefail

# --- Configuration ---
MAX_JOBS=3
POLLING_INTERVAL=100
JOB_NAME_PREFIX="isoquant_"


BAM_DIR="data/bams"
JOB_SCRIPT="scripts/preprocessing_pipeline/01_run_isoquant/isoquant.sh"

for f in "${BAM_DIR}"/*.rmRNA.bam; do
    [ -e "$f" ] || { echo "Warning: No .rmRNA.bam files found in ${BAM_DIR}. Exiting."; break; }

    root_name=$(basename "$f")
    label=$(echo "$root_name" | sed 's/\.rmRNA\.bam$//')

    while true; do
        running_jobs=$(qstat -u "$USER" | grep "$JOB_NAME_PREFIX" | grep -c '[QR]' || true)

        if [[ $running_jobs -lt $MAX_JOBS ]]; then
            echo "Slot available. Submitting job for ${label}."
            qsub -N "${JOB_NAME_PREFIX}${label}" -v bam="$f",p1="$label" "$JOB_SCRIPT"
            sleep 1
            break
        else
            echo "Queue is full (${running_jobs}/${MAX_JOBS}). Waiting..."
            sleep "$POLLING_INTERVAL"
        fi
    done
done