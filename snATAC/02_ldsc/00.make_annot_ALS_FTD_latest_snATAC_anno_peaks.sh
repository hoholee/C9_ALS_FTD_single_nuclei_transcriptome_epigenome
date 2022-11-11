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

bed_path="/cndd3/junhao/ALS_FTD_singleCell/snATAC_ArchR/ALS_FTD_clean/ldsc/peaks"

for bed_file in ${bed_path}/*Hg19.bed
do

	sample=$(basename -s _liftToHg19.bed ${bed_file})
	echo ${sample}

	for chr in {1..22}
	do
		python ${LDSC_PATH}/make_annot.py \
			--bed-file ${bed_file} \
			--bimfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
			--annot-file ./annot_files/${sample}.${chr}.annot.gz
	done
done

echo "Ending at $(date)"
