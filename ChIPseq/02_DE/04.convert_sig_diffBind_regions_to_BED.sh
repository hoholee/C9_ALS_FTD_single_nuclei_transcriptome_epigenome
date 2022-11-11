#!/bin/bash

## logging which server and which conda environment this scirpt is running in
echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"

for i in diffBind_H3K27ac_ALS_vs_Contro*Term_res*.txt
do
	j=${i/.txt/.bed}
	awk 'BEGIN{FS="\t";OFS="\t"} \
		NR>1 && $11<0.05 { \
			if($9>0){color = "225,49,94"}else{color = "32,104,144"}; \
			name=FILENAME; sub(/\.txt/, "", name); sub(/^diffBind_H3K27ac_/, "", name); \
			score = -log($11)/log(10); \
			printf "%s\t%d\t%d\t%s_sig_%d|%.2f|%.2f|%.2f\t%.2f\t%s\t%d\t%d\t%s\n", $1, $2-1, $3, name, NR, $7, $8, $9, score, ".", $2-1, $2-1, color \
		}' $i > $j
done

echo "Ending at $(date)"
