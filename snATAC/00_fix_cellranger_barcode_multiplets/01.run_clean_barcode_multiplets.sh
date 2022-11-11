#!/bin/bash

## Run `clean_barcode_multiplets_1.1.py` provided by 10X genomics support team to remove barcode multiplets
## activate the `master` conda environment before running this script
# conda activate master

## logging which server and which conda environment this scirpt is running in
echo "Running this on ### "`hostname`" ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at "`date`
echo "cellranger-atac path: "`which cellranger-atac`

# Provide the cleanup tool access to the Cell Ranger ATAC python libraries
source /home/jul307/software/cellranger_atac/cellranger-atac-1.1.0/sourceme.bash
echo "python path: "`which python`

clean_sh="/home/jul307/software/cellranger_atac/cellranger-atac-1.1.0/clean_barcode_multiplets_1.1.py"
old_res_path="/cndd/junhao/ALS_FTD_singleCell/atac_seq_analysis/fix_cellranger_barcode_multiplets/cellranger_results"

for i in `ls -d ${old_res_path}/ATAC*`;
do
	echo "Running cleaning process for sample: "`basename $i`;
	python ${clean_sh} \
		--output_path . \
		--prefix `basename $i`_ \
		$i/outs
done

echo "All Done!"
echo "Ending at "`date`
