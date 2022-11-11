#!/bin/bash

## run trimming

## activate the `master` conda environment before running this script
# conda activate master

i=$1
echo "Handling $i"
j=${i/R1.fastq.gz/R2.fastq.gz}
log=$(basename $i)_trimming.log
trim_galore --cores 4 --paired --clip_R1 3 -q 20 --fastqc_args "-noextract" $i $j > $log 2>&1

echo "Ending $i at $(date)"
