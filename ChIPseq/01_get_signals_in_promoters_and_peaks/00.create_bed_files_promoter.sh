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
input_bed="Homo_sapiens.GRCh38.93.filtered.bed"
chrom_size="/cndd/junhao/genomes/hg38/chrom_size.txt"
upstream_bps=2000
downstream_bps=1000
output_bed="promoter_TSS_upstream${upstream_bps}_downstream${downstream_bps}.bed"

awk 'BEGIN{FS = "\t"; OFS = "\t"} $6 == "+"{print $1, $2, $2 + 1, $4, $5, $6} $6 == "-"{print $1, $3 - 1, $3, $4, $5, $6}' ${input_bed} | \
	bedtools slop \
		-l ${upstream_bps} \
		-r ${downstream_bps} \
		-s \
		-i - \
		-g ${chrom_size} \
			> ${output_bed}

echo "Ending at $(date)"
