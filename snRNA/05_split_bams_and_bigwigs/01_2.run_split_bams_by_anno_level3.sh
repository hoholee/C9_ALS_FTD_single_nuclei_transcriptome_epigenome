#!/bin/bash

## split bam files by cluster assigned by cell ids

## activate the `master` conda environment before running this script
# conda activate master

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "samtools path: $(command -v samtools)"

cat sample_list.txt | \
	xargs -P 3 -I {} ./01_0.split_bam_by_cluster.sh {} rna_anno_level_3 level_3 4

echo "Ending at $(date)"
