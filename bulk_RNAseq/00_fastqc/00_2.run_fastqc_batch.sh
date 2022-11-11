#!/bin/bash

## run fastqc

## activate the `master` conda environment before running this script
# conda activate master

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "fastqc path: $(command -v fastqc)"

raw_fastq_path="/cndd2/junhao/ALS_FTD_bulk/rna_seq/ALS_FTD_batch/raw"

find -L ${raw_fastq_path} -name "*.fastq.gz" | sort | xargs -n 1 -P 16 -I {} ./00_1.fastqc_single_sample.sh {}

echo "Ending at $(date)"
