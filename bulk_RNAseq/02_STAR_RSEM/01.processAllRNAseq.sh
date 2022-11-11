#!/bin/bash

# Run STAR and RSEM for each sample

## activate the `master` conda environment before running this script
# conda activate master

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "STAR path: $(command -v STAR)"

rnaseq_path="/cndd2/junhao/ALS_FTD_bulk/rna_seq/ALS_FTD_batch/clean_reads"

STARgenomeDir="/cndd2/junhao/genome/hg38/star_rsem_idx/star_idx"
RSEMrefDir="/cndd2/junhao/genome/hg38/star_rsem_idx/rsem_idx/RSEMref_hg38_gencodeV35"
dataType="str_PE" # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
nThreadsSTAR=16 # number of threads for STAR
nThreadsRSEM=16 # number of threads for RSEM

for i in `cat $rnaseq_path/sample_list.txt`
do
        echo "Handling "$i
        rm -rf $i
        mkdir $i && cd $i
        j=$rnaseq_path/combine_lanes/merged_${i}_R1.fastq.gz
        k=$rnaseq_path/combine_lanes/merged_${i}_R2.fastq.gz
        ../00.STAR_RSEM.sh $j $k $STARgenomeDir $RSEMrefDir $dataType $nThreadsSTAR $nThreadsRSEM
        cd ..
done

STAR --genomeDir $STARgenomeDir --genomeLoad Remove
echo "All Done!"
