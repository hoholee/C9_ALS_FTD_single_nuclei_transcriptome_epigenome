#!/bin/bash

## get the TSS from a gene annotation file in BED format  and then get flanking regions as the promoter regions
## activate the `master` conda environment before running this script
# conda activate master

## logging which server and which conda environment this scirpt is running in
echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "bedtools path: $(command -v bedtools)"

# starts with the bed annotation file used in the 10X Cellranger snRNA-seq pipeline (Ensembl V93 filtered)
ref_bed="Homo_sapiens.GRCh38.93.filtered.bed"

for i in diffBind*Term_res*.bed
do
	j=${i/.bed/_closestGenes.txt}
	bedtools closest \
		-a <(bedtools sort -i "${i}") \
		-b <(bedtools sort -i "${ref_bed}") \
		-D ref > "${j}"

done

echo "Ending at $(date)"
