#!/bin/bash

## combine files run in the two different lanes for each biological sample

## activate the `master` conda environment before running this script
# conda activate master

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"

for i in $(cat sample_list.txt)
do
	echo "Handling $i";
	R1_files=()
	while IFS=  read -r -d $'\0'; do
		R1_files+=("$REPLY")
	done < <(find . -name "*${i}*R1_val_1.fq.gz" -print0 | sort -z)
	R1_lane1=${R1_files[0]}
	R1_lane2=${R1_files[1]}
	R2_lane1=${R1_lane1/R1_val_1.fq.gz/R2_val_2.fq.gz}
	R2_lane2=${R1_lane2/R1_val_1.fq.gz/R2_val_2.fq.gz}
	cat ${R1_lane1} ${R1_lane2} > ./combine_lanes/merged_${i}_R1.fastq.gz
	cat ${R2_lane1} ${R2_lane2} > ./combine_lanes/merged_${i}_R2.fastq.gz
done

echo "Ending at $(date)"
