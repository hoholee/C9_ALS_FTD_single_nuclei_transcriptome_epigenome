# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(Augur)
library(tictoc)

tic()
data_obj_sub <- readRDS("../seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_SeuratV4_object.rds")
toc()

meta_data <- data_obj_sub@meta.data %>% as_tibble()

selected_n_subsamples <- 50
selected_subsample_size <- 200
selected_folds <- 3

all_diseases <- c("ALS", "FTD", "Control")
excluded_cell_types <- c("NK_cell", "Exc_unknown", "Ambiguous")

run_Augur <- function(selected_region, selected_disease_1, selected_disease_2) {
  message(
    str_glue(
      "Running Augur for {selected_region} {selected_disease_1} vs {selected_disease_2}..."
    )
  )

  selected_diseases <- c(selected_disease_1, selected_disease_2)
  excluded_disease <- all_diseases[!all_diseases %in% selected_diseases]

  # subset Seurat object
  sub_idx <- which(
    !(meta_data$rna_anno_2ndRound_level_2 %in% excluded_cell_types) &
      meta_data$region == selected_region &
      meta_data$disease != excluded_disease
  )
  data_obj_used <- data_obj_sub[, sub_idx]

  # run Augur
  tic()
  auc_out <- calculate_auc(
    data_obj_used,
    label_col = "disease",
    cell_type_col = "rna_anno_2ndRound_level_2",
    n_subsamples = selected_n_subsamples, # default 50
    subsample_size = selected_subsample_size, # default 20
    folds = selected_folds, # default 3
    min_cells = NULL,
    var_quantile = 0.5,
    feature_perc = 0.5,
    n_threads = 16,
    show_progress = TRUE,
    augur_mode = "default",
    classifier = "rf",
    rf_params = list(trees = 100, mtry = 2, min_n = NULL, importance = "accuracy")
  )
  toc()

  # save output
  saveRDS(
    auc_out,
    str_glue(
      "augur_out_obj_",
      "{selected_disease_1}_vs_{selected_disease_2}_",
      "{selected_region}_level2_",
      "nSubsamples_{selected_n_subsamples}_",
      "subsampleSize_{selected_subsample_size}_",
      "folds_{selected_folds}.rds"
    )
  )
  write_tsv(
    auc_out$AUC,
    str_glue(
      "augur_out_obj_",
      "{selected_disease_1}_vs_{selected_disease_2}_",
      "{selected_region}_level2_",
      "nSubsamples_{selected_n_subsamples}_",
      "subsampleSize_{selected_subsample_size}_",
      "folds_{selected_folds}.tsv"
    )
  )
  rm(data_obj_used)
}

# run scDist in groups
params <- expand_grid(
  selected_region = c("MCX", "mFCX"),
  selected_disease_1 = c("ALS", "FTD"),
  selected_disease_2 = c("Control")
)

pwalk(params, run_Augur)

## log sessionInfo
sessionInfo()
