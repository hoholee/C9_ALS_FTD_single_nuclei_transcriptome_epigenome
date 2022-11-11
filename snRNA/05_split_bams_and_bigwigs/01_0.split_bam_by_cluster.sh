#!/bin/bash

## split bam files by cluster assigned by cell ids

## activate the `master` conda environment before running this script
# conda activate master

sample=$1
split_by_column=$2
save_to_path_prefix=$3
num_thread=$4

raw_bam_path="/cndd2/junhao/ALS_FTD_singleCell/snRNA_postCellBender/split_bams/raw_bam"
cell_meta="cell_meta_data_${sample}.txt"
input_bam="${raw_bam_path}/${sample}.bam"
output_header="${save_to_path_prefix}/${sample}.header"

# saving header
samtools view -H ${input_bam} > ${output_header}

# spliting file
awk -v split_by_column=$split_by_column -v save_to_path_prefix=$save_to_path_prefix -v sample=$sample ' \
	BEGIN {FS = "\t"; OFS = "\t"} \
	NR==FNR && FNR==1 {for (i==1; i <=NF; i++){idx[$i] = i}} \
	NR==FNR && FNR>1 {cell_to_cluster[$1] = $idx[split_by_column]} \
	NR>FNR {
		for(j=12; j<=NF; j++) { \
			if($j~/^CB/){ \
				split($j, read_cell_id, ":"); \
				if(read_cell_id[3] in cell_to_cluster) {print $0 >> save_to_path_prefix"/"sample"_cluster_"cell_to_cluster[read_cell_id[3]]"_noheader.sam"}; \
				break \
				} \
			} \
		} \
	' ${cell_meta} <(samtools view ${input_bam})

# adding back header and save to bam files
for i in ${save_to_path_prefix}/${sample}_cluster*_noheader.sam
do
	j=${i/_noheader.sam/.bam}
	samtools view -@ ${num_thread} -b -o $j <(cat ${save_to_path_prefix}/${sample}.header $i)
done

# cleaning up
rm ${output_header}
rm ${save_to_path_prefix}/${sample}_cluster*_noheader.sam
