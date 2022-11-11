#!/bin/bash

## split cell meta data by sample

## activate the `master` conda environment before running this script
# conda activate master

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"

awk 'BEGIN {FS="\t";OFS="\t"} \
	NR>1 && $23!="Ambiguous" && $23!="NK_cell" && $23!="Exc_unknown"{ \
		split($1, cell_id, "_"); \
		if(!seen[$2]){ \
			seen[$2] = 1; \
			print "barcode", "rna_anno_level_1", "rna_anno_level_2", "rna_anno_level_3" > "cell_meta_data_"$2".txt" \
		};  \
		print cell_id[4], $22, $23, $24 >> "cell_meta_data_"$2".txt" \
	}' metadata_all_cells_2nd_round_annotations.txt

echo "Ending at $(date)"
