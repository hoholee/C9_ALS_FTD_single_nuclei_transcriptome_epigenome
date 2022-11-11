#!/bin/bash

## summarized QC stats with multiqc

## activate the `multiqc` conda environment before running this script
# conda activate multiqc

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "multiqc path: $(command -v multiqc)"

## run multiqc
multiqc -n ALS_FTD_bulk_RNAseq_hg38_mapping_fastqc_report .

echo "Ending at $(date)"
