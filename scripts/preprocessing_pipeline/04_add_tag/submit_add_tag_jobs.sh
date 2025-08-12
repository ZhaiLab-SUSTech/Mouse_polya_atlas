#!/bin/bash
###
 # @Description: Add tag for final samples; control submitted tasks num
 # @Date: 2024-06-30 21:34:57
 # @LastEditTime: 2024-10-30 15:53:19
 # @author: leihj
###

max_jobs=10
count=1
for f in bamfile/*.rmRNA.bam; do
    root=$(basename $f)
    label=$(echo $root | sed 's/.rmRNA.bam//g')


    while true; do

        running_jobs=$(qstat -a | grep leihj | grep 'add_tag_' | grep -c '[QR]')
        # running_jobs=$(qstat -a | grep leihj | grep 'test' | grep -c '[QR]') # for test


        if [[ $running_jobs -lt $max_jobs ]]; then
            echo $running_jobs
            qsub -N add_tag_$label -v inbam=$root add_tag_iso.sh
        #     qsub -N test sleep.sh # for test
            echo "Submitted job for $label"
            echo "$count jobs have been submitted"
            ((count++))
            sleep 1
            break
        else
            echo "Waiting for available slots..."
            sleep 100
            # sleep 5
        fi
    done
done