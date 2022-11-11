#!/bin/bash

## use wiggletools to compute the mean signal across individuals
## and convert output wiggle files to bigwigs

## activate the `wiggletools` conda environment before running this script
# conda activate wiggletools

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "bamCoverage path: $(command -v wiggletools)"

for cell_type in Neurons Oligodendrocytes Other_glias
do
	echo "Processing cell type: ${cell_type}..."
	for region in mid_frontal_cortex motor_cortex
	do
		echo "Region: ${region}..."
		for diagnosis in Control ALS FTD
		do
			echo "Taking average of diagnosis group: ${diagnosis}..."
			wiggletools write ${region}_${diagnosis}_${cell_type}_avg_RPKM.wig \
				mean ./*_${region}_${cell_type}_${diagnosis}.bw
			./wigToBigWig ${region}_${diagnosis}_${cell_type}_avg_RPKM.wig chrom_size.txt ${region}_${diagnosis}_${cell_type}_avg_RPKM.bw
		done
	done
done

echo "Ending at $(date)"
