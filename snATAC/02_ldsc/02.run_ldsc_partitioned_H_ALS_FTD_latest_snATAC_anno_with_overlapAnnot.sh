#!/bin/bash

## Create annot files
## activate the `ldsc` conda environment before running this script
# conda activate ldsc

## logging which server and which conda environment this scirpt is running in
echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "python path: $(which python)"
LDSC_PATH="/home/jul307/software/ldsc"
echo "ldsc path: $LDSC_PATH"

## bypass the "LinAlgError: Singular matrix" error by adding the `--n-blocks 1000` parameter

cts_name="ALS_FTD_latest_snATAC_anno"
GWAS_trait_path="/cndd2/junhao/ldsc_tutorial/GWAS_source/sumstats"

for i in ${GWAS_trait_path}/*.sumstats.gz
do
	GWAS_trait=$(basename -s .sumstats.gz $i)
	echo ${GWAS_trait}

	python ${LDSC_PATH}/ldsc.py \
		--h2-cts ${GWAS_trait_path}/${GWAS_trait}.sumstats.gz \
		--out ./ldsc_results/${GWAS_trait}_${cts_name}_withOverlapAnnot \
		--ref-ld-chr ./1000G_EUR_Phase3_baseline/baseline. \
		--overlap-annot \
		--frqfile-chr ./1000G_Phase3_frq/1000G.EUR.QC. \
		--ref-ld-chr-cts ${cts_name}.ldcts \
		--w-ld-chr ./1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
		#--n-blocks 1000
done

echo "Ending at $(date)"
