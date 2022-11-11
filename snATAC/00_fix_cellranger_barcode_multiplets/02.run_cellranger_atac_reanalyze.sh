#!/bin/bash

## Run cellranger-atac reanalyze, only used the valid barcodes 
## activate the `master` conda environment before running this script
# conda activate master

## logging which server and which conda environment this scirpt is running in
echo "Running this on ### "`hostname`" ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at "`date`
echo "cellranger-atac path: "`which cellranger-atac`

old_res_path="/cndd/junhao/ALS_FTD_singleCell/atac_seq_analysis/fix_cellranger_barcode_multiplets/cellranger_results"
reference_path="/cndd/junhao/genomes/cellRanger_references/refdata-cellranger-atac-GRCh38-1.1.0"

for i in `ls -d ${old_res_path}/ATAC*`;
do
	echo "Running reanalyze process for sample: "`basename $i`" at "`date`;
	cellranger-atac reanalyze \
		--id=`basename $i`"_reanalyze" \
		--peaks=${i}/outs/peaks.bed \
		--fragments=${i}/outs/fragments.tsv.gz \
		--reference=${reference_path} \
		--barcodes=`basename $i`"_cell_barcodes.csv" \
		--localcores=32 \
		--localmem=128
done

echo "All Done!"
echo "Ending at "`date`
