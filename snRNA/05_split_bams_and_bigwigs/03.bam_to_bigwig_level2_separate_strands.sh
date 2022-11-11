#!/bin/bash

## convert bam to bigwig, separating plus and minus strands

## activate the `master` conda environment before running this script
# conda activate master

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "bamCoverage path: $(command -v bamCoverage)"

for i in ./level_2/*.bam;
do
	echo "Processing ${i}..."
	j=${i/.bam/_forward.bw}
	k=${i/.bam/_reverse.bw}

	echo "Processing forward strand..."
	bamCoverage \
		-b ${i} \
		-o ${j} \
		--binSize 50 \
		--numberOfProcessors 16 \
		--normalizeUsing RPKM \
		--filterRNAstrand forward \
		--ignoreDuplicates

	echo "Processing reverse strand..."
	bamCoverage \
		-b ${i} \
		-o ${k} \
		--binSize 50 \
		--numberOfProcessors 16 \
		--normalizeUsing RPKM \
		--filterRNAstrand reverse \
		--ignoreDuplicates
done

echo "Ending at $(date)"
