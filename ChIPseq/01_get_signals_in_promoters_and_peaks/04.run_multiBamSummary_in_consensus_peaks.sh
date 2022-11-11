#!/bin/bash

## activate the `master` conda environment before running this script
# conda activate master

# usage: ./04.run_multiBamSummary_in_consensus_peaks.sh cell_type
echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "multiBamSummary path: $(command -v multiBamSummary)"

cell_type=$1
bam_path="/cndd3/Public_Datasets/Dracheva_ALS_ChIPseq/bam"
input_bed="consensus_peaks_${cell_type}.bed"
blacklist_bed="/cndd/junhao/genomes/encode_blacklist/hg38-blacklist.v2.bed.gz"
output_npz="H3K27ac_ChIPseq_raw_counts_in_${cell_type}_consensus_peaks.npz"
output_raw="H3K27ac_ChIPseq_raw_counts_in_${cell_type}_consensus_peaks.tsv"
output_scaling_factors="H3K27ac_ChIPseq_raw_counts_in_${cell_type}_consensus_peaks_scaling_factors.tsv"
min_map_Q=20
num_proc=16

multiBamSummary BED-file \
	--BED ${input_bed} \
	--bamfiles \
		${bam_path}/113-Astro.merged.bam \
		${bam_path}/113-MG.merged.bam \
		${bam_path}/113-NeuN.merged.bam \
		${bam_path}/113-Olig.merged.bam \
		${bam_path}/114-Astro.merged.bam \
		${bam_path}/114-MG.merged.bam \
		${bam_path}/114-NeuN.merged.bam \
		${bam_path}/114-Olig.merged.bam \
		${bam_path}/115-Astro.merged.bam \
		${bam_path}/115-MG.merged.bam \
		${bam_path}/115-NeuN.merged.bam \
		${bam_path}/115-Olig.merged.bam \
		${bam_path}/128-Astro.merged.bam \
		${bam_path}/128-MG.merged.bam \
		${bam_path}/128-NeuN.merged.bam \
		${bam_path}/128-Olig.merged.bam \
		${bam_path}/129-Astro.merged.bam \
		${bam_path}/129-MG.merged.bam \
		${bam_path}/129-NeuN.merged.bam \
		${bam_path}/129-Olig.merged.bam \
		${bam_path}/130-Astro.merged.bam \
		${bam_path}/130-MG.merged.bam \
		${bam_path}/130-NeuN.merged.bam \
		${bam_path}/130-Olig.merged.bam \
		${bam_path}/302-Astro.merged.bam \
		${bam_path}/302-MG.merged.bam \
		${bam_path}/302-NeuN.merged.bam \
		${bam_path}/302-Olig.merged.bam \
		${bam_path}/303-Astro.merged.bam \
		${bam_path}/303-MG.merged.bam \
		${bam_path}/303-NeuN.merged.bam \
		${bam_path}/303-Olig.merged.bam \
		${bam_path}/309-Astro.merged.bam \
		${bam_path}/309-MG.merged.bam \
		${bam_path}/309-NeuN.merged.bam \
		${bam_path}/309-Olig.merged.bam \
		${bam_path}/311-Astro.merged.bam \
		${bam_path}/311-MG.merged.bam \
		${bam_path}/311-NeuN.merged.bam \
		${bam_path}/311-Olig.merged.bam \
		${bam_path}/312-Astro.merged.bam \
		${bam_path}/312-MG.merged.bam \
		${bam_path}/312-NeuN.merged.bam \
		${bam_path}/312-Olig.merged.bam \
		${bam_path}/320-Astro.merged.bam \
		${bam_path}/320-MG.merged.bam \
		${bam_path}/320-NeuN.merged.bam \
		${bam_path}/320-Olig.merged.bam \
		${bam_path}/322-Astro.merged.bam \
		${bam_path}/322-MG.merged.bam \
		${bam_path}/322-NeuN.merged.bam \
		${bam_path}/325-Astro.merged.bam \
		${bam_path}/325-MG.merged.bam \
		${bam_path}/325-NeuN.merged.bam \
		${bam_path}/325-Olig.merged.bam \
		${bam_path}/328-Astro.merged.bam \
		${bam_path}/328-MG.merged.bam \
		${bam_path}/328-NeuN.merged.bam \
		${bam_path}/328-Olig.merged.bam \
	--outFileName ${output_npz} \
	--outRawCounts ${output_raw} \
	--scalingFactors ${output_scaling_factors} \
	--labels \
		ALS_113_Astro \
		ALS_113_MG \
		ALS_113_NeuN \
		ALS_113_Olig \
		ALS_114_Astro \
		ALS_114_MG \
		ALS_114_NeuN \
		ALS_114_Olig \
		ALS_115_Astro \
		ALS_115_MG \
		ALS_115_NeuN \
		ALS_115_Olig \
		ALS_128_Astro \
		ALS_128_MG \
		ALS_128_NeuN \
		ALS_128_Olig \
		ALS_129_Astro \
		ALS_129_MG \
		ALS_129_NeuN \
		ALS_129_Olig \
		ALS_130_Astro \
		ALS_130_MG \
		ALS_130_NeuN \
		ALS_130_Olig \
		Control_302_Astro \
		Control_302_MG \
		Control_302_NeuN \
		Control_302_Olig \
		Control_303_Astro \
		Control_303_MG \
		Control_303_NeuN \
		Control_303_Olig \
		Control_309_Astro \
		Control_309_MG \
		Control_309_NeuN \
		Control_309_Olig \
		Control_311_Astro \
		Control_311_MG \
		Control_311_NeuN \
		Control_311_Olig \
		Control_312_Astro \
		Control_312_MG \
		Control_312_NeuN \
		Control_312_Olig \
		Control_320_Astro \
		Control_320_MG \
		Control_320_NeuN \
		Control_320_Olig \
		Control_322_Astro \
		Control_322_MG \
		Control_322_NeuN \
		Control_325_Astro \
		Control_325_MG \
		Control_325_NeuN \
		Control_325_Olig \
		Control_328_Astro \
		Control_328_MG \
		Control_328_NeuN \
		Control_328_Olig \
	--blackListFileName ${blacklist_bed} \
	--extendReads \
	--ignoreDuplicates \
	--minMappingQuality ${min_map_Q} \
	--centerReads \
	--numberOfProcessors ${num_proc}

echo "Ending at $(date)"
