#!/bin/bash

## run cellBender to remove ambient RNA background
## run this on djembe/tabla where GPU and CUDA were setup correctly
## activate the `CellBender_GPU` conda environment before running this script
# conda activate /home/AD/jchien/.conda/envs/CellBender_GPU

## logging which server and which conda environment this scirpt is running in
echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "cellbender path: $(command -v cellbender)"

# input/output settings
raw_h5_path="/cndd2/junhao/ALS_FTD_singleCell/run_cellBender_on_raw_snRNA/raw_cellRanger_output"

# model training parameters
expected_cells=1500
total_droplets_included=18000
low_count_threshold=15
num_epochs=150
learning_rate=0.0001
false_positive_rate=0.01
# lower the following two if GPU is running out of memory
posterior_batch_size=5
cells_posterior_reg_calc=50

for sample_id in $(cat sample_id_mFCX.txt)
do
	input_h5="${raw_h5_path}/${sample_id}_raw_feature_bc_matrix.h5"
	output_h5="${sample_id}_cellBender_corrected.h5"
	echo "Running sample: ${sample_id} in ${input_h5}"

	cellbender remove-background \
		--input ${input_h5} \
		--output ${output_h5} \
		--expected-cells ${expected_cells} \
		--total-droplets-included ${total_droplets_included} \
		--model full \
		--epochs ${num_epochs} \
		--cuda \
		--low-count-threshold ${low_count_threshold} \
		--fpr ${false_positive_rate} \
		--learning-rate ${learning_rate} \
		--posterior-batch-size ${posterior_batch_size} \
		--cells-posterior-reg-calc ${cells_posterior_reg_calc}
done

echo "Ending at $(date)"
