#!/bin/sh

# Generate data matrix using "expected counts" from RSEM results for EBSeq

path_genes=""

for i in `find . -name "*.genes.results" | sort`
do
	temp=$path_genes""$i" "
	path_genes=$temp
done


#echo $path_genes

/home/jul307/software/RSEM-1.2.30/rsem-generate-data-matrix $path_genes > allSamples_rsem_genes_results.txt
./generate-data-matrix-tpm $path_genes > allSamples_rsem_genes_results_TPM.txt
