#!/bin/bash

## convert bam to bigwig

## activate the `master` conda environment before running this script
# conda activate master

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "samtools path: $(command -v samtools)"
echo "bamCoverage path: $(command -v bamCoverage)"

for i in ./level_2/*.bam;
do
	echo "Processing ${i}..."
	j=${i/.bam/.bw}
	samtools index ${i}
	bamCoverage \
		-b ${i} \
		-o ${j} \
		--binSize 50 \
		--numberOfProcessors 16 \
		--normalizeUsing RPKM \
		--ignoreDuplicates
done

echo "Ending at $(date)"
