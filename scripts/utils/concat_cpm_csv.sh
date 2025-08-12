#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=100G
#PBS -j oe

cd $PBS_O_WORKDIR

python3 concat_cpm.py
