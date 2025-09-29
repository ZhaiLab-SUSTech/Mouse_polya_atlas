#!/bin/bash
###
 # @Description: Add tag for final samples; control submitted tasks num
 # @Date: 2024-06-30 21:34:57
 # @LastEditTime: 2024-10-30 15:53:19
 # @author: leihj
###

MAX_JOBS=10
POLLING_INTERVAL=100

MERGED_DIR="results/02_merged_readinfo"

for f in bamfile/*.rmRNA.bam; do
    root=$(basename $f)
    label=$(echo $root | sed 's/.rmRNA.bam//g')
    tissue=$(echo $label | sed 's/_rep.*//')

    while true; do

        running_jobs=$(qstat -a | grep leihj | grep 'add_tag_' | grep -c '[QR]')
        # running_jobs=$(qstat -a | grep leihj | grep 'test' | grep -c '[QR]') # for test

        if [[ $running_jobs -lt $MAX_JOBS ]]; then
            read_info="${MERGED_DIR}/${tissue}/${label}.parquet"

            echo "Curren00t running jobs: $running_jobs. Submitting a new job for ${label}."
            
            qsub -N "add_tag_${label}" \
                 -v "SAMPLE_NAME=${label},INPUT_BAM=${f},READ_INFO_PARQUET=${read_info}" \
                 "$EXECUTION_SCRIPT"

            sleep 1
            break
        else
            echo "Queue is full (${running_jobs}/${MAX_JOBS}). Waiting..."
            sleep $POLLING_INTERVAL
            # sleep 5
        fi
    done
done