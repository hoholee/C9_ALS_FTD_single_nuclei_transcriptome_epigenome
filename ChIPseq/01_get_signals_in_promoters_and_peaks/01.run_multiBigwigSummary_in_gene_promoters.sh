#!/bin/bash

## activate the `master` conda environment before running this script
# conda activate master

echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "multiBigwigSummary path: $(command -v multiBigwigSummary)"

bigwig_path="/cndd2/junhao/ALS_FTD_singleCell/ChIPseq/bigwig"
input_bed="promoter_TSS_upstream2000_downstream1000.bed"
blacklist_bed="/cndd/junhao/genomes/encode_blacklist/hg38-blacklist.v2.bed.gz"
output_npz="H3K27ac_signal_in_promoter_TSS_upstream2000_downstream1000.npz"
output_raw="H3K27ac_signal_in_promoter_TSS_upstream2000_downstream1000.tsv"
num_proc=16

multiBigwigSummary BED-file \
	--BED ${input_bed} \
	--bwfiles \
		${bigwig_path}/113-Astro.1M.bw \
		${bigwig_path}/113-MG.1M.bw \
		${bigwig_path}/113-NeuN.1M.bw \
		${bigwig_path}/113-Olig.1M.bw \
		${bigwig_path}/114-Astro.1M.bw \
		${bigwig_path}/114-MG.1M.bw \
		${bigwig_path}/114-NeuN.1M.bw \
		${bigwig_path}/114-Olig.1M.bw \
		${bigwig_path}/115-Astro.1M.bw \
		${bigwig_path}/115-MG.1M.bw \
		${bigwig_path}/115-NeuN.1M.bw \
		${bigwig_path}/115-Olig.1M.bw \
		${bigwig_path}/128-Astro.1M.bw \
		${bigwig_path}/128-MG.1M.bw \
		${bigwig_path}/128-NeuN.1M.bw \
		${bigwig_path}/128-Olig.1M.bw \
		${bigwig_path}/129-Astro.1M.bw \
		${bigwig_path}/129-MG.1M.bw \
		${bigwig_path}/129-NeuN.1M.bw \
		${bigwig_path}/129-Olig.1M.bw \
		${bigwig_path}/130-Astro.1M.bw \
		${bigwig_path}/130-MG.1M.bw \
		${bigwig_path}/130-NeuN.1M.bw \
		${bigwig_path}/130-Olig.1M.bw \
		${bigwig_path}/302-Astro.1M.bw \
		${bigwig_path}/302-MG.1M.bw \
		${bigwig_path}/302-NeuN.1M.bw \
		${bigwig_path}/302-Olig.1M.bw \
		${bigwig_path}/303-Astro.1M.bw \
		${bigwig_path}/303-MG.1M.bw \
		${bigwig_path}/303-NeuN.1M.bw \
		${bigwig_path}/303-Olig.1M.bw \
		${bigwig_path}/309-Astro.1M.bw \
		${bigwig_path}/309-MG.1M.bw \
		${bigwig_path}/309-NeuN.1M.bw \
		${bigwig_path}/309-Olig.1M.bw \
		${bigwig_path}/311-Astro.1M.bw \
		${bigwig_path}/311-MG.1M.bw \
		${bigwig_path}/311-NeuN.1M.bw \
		${bigwig_path}/311-Olig.1M.bw \
		${bigwig_path}/312-Astro.1M.bw \
		${bigwig_path}/312-MG.1M.bw \
		${bigwig_path}/312-NeuN.1M.bw \
		${bigwig_path}/312-Olig.1M.bw \
		${bigwig_path}/320-Astro.1M.bw \
		${bigwig_path}/320-MG.1M.bw \
		${bigwig_path}/320-NeuN.1M.bw \
		${bigwig_path}/320-Olig.1M.bw \
		${bigwig_path}/322-Astro.1M.bw \
		${bigwig_path}/322-MG.1M.bw \
		${bigwig_path}/322-NeuN.1M.bw \
		${bigwig_path}/325-Astro.1M.bw \
		${bigwig_path}/325-MG.1M.bw \
		${bigwig_path}/325-NeuN.1M.bw \
		${bigwig_path}/325-Olig.1M.bw \
		${bigwig_path}/328-Astro.1M.bw \
		${bigwig_path}/328-MG.1M.bw \
		${bigwig_path}/328-NeuN.1M.bw \
		${bigwig_path}/328-Olig.1M.bw \
	--outFileName ${output_npz} \
	--outRawCounts ${output_raw} \
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
	--numberOfProcessors ${num_proc}

echo "Ending at $(date)"
