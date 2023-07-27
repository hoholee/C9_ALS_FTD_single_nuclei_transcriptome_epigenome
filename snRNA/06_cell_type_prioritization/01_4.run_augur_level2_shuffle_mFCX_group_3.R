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

# set up random group by shuffle
indi_ALS <- c("332", "111", "113", "52", "388", "110")
indi_Control <- c("945", "902", "91", "1069", "904", "906")
indi_FTD <- c("36", "54", "55", "674", "61", "908") # subject 908 was excluded from the seurat object

indi_dict <- list(
  "ALS" = indi_ALS,
  "Control" = indi_Control,
  "FTD" = indi_FTD
)

random_group_dict <- list(
  "ALS_vs_Control" = list(
    "half" = combn(indi_ALS, 3, simplify = TRUE)[, c(1:10)],
    "full" = combn(indi_Control, 3, simplify = TRUE)
  ),
  "FTD_vs_Control" = list(
    "half" = combn(indi_FTD, 3, simplify = TRUE)[, c(1:10)],
    "full" = combn(indi_Control, 3, simplify = TRUE)
  ),
  "ALS_vs_FTD" = list(
    "half" = combn(indi_ALS, 3, simplify = TRUE)[, c(1:10)],
    "full" = combn(indi_FTD, 3, simplify = TRUE)
  )
)

combin_idx <- expand.grid(idx_1 = c(1:10), idx_2 = c(1:20)) %>% as_tibble()

# get seed for random groups that contain same number of subjects from both sexes
# !!!! only works for ALS vs Control for now !!!!
subject_groups <- read_tsv("../DE_MAST/second_round_level_2/shuffle_seed_subject_groups.tsv", col_types = "icc")
num_subject_per_group <- read_tsv("../DE_MAST/second_round_level_2/shuffle_seed_subject_groups_num_subject_per_sex_group1.tsv", col_types = "icii")

seed_list <- num_subject_per_group %>%
  filter(F == 3) %>%
  pull(seed)

# Augur parameters
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

  # get the meta data in the subset cells
  data_cell_meta <- data_obj_used@meta.data

  # shuffle the disease labels
  subject_disease_match <- data_cell_meta %>% distinct(disease, subject)

  run_shuffle <- function(selected_seed) {
    message(paste0("Iteration: ", selected_seed))

    comp <- paste0(selected_disease_1, "_vs_", selected_disease_2)
    group_A_1 <- random_group_dict[[comp]][["half"]][, combin_idx$idx_1[selected_seed]]
    group_A_2 <- random_group_dict[[comp]][["full"]][, combin_idx$idx_2[selected_seed]]
    group_B_1 <- indi_dict[[selected_disease_1]][!indi_dict[[selected_disease_1]] %in% group_A_1]
    group_B_2 <- indi_dict[[selected_disease_2]][!indi_dict[[selected_disease_2]] %in% group_A_2]

    full_permuted_df <- tibble(
      subject = c(group_A_1, group_A_2, group_B_1, group_B_2),
      disease_permute = rep(c(selected_disease_1, selected_disease_2), each = 6)
    )

    subject_disease_match_permuted <- subject_disease_match %>%
      left_join(full_permuted_df, by = "subject")

    disease_permute_dict <- subject_disease_match_permuted$disease_permute
    names(disease_permute_dict) <- subject_disease_match_permuted$subject

    data_obj_used$disease_permute <- disease_permute_dict[data_obj_used$subject]

    # run Augur
    tic()
    auc_out <- calculate_auc(
      data_obj_used,
      label_col = "disease_permute",
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
        "./Augur_shuffle/augur_out_obj_",
        "{selected_disease_1}_vs_{selected_disease_2}_",
        "{selected_region}_level2_",
        "nSubsamples_{selected_n_subsamples}_",
        "subsampleSize_{selected_subsample_size}_",
        "folds_{selected_folds}_",
        "seed_{selected_seed}.rds"
      )
    )
    write_tsv(
      auc_out$AUC,
      str_glue(
        "./Augur_shuffle/augur_out_obj_",
        "{selected_disease_1}_vs_{selected_disease_2}_",
        "{selected_region}_level2_",
        "nSubsamples_{selected_n_subsamples}_",
        "subsampleSize_{selected_subsample_size}_",
        "folds_{selected_folds}_",
        "seed_{selected_seed}.tsv"
      )
    )
    rm(data_obj_used)
  }

  walk(seed_list[43:63], run_shuffle)
}

# run scDist in groups
params <- expand_grid(
  selected_region = c("mFCX"),
  selected_disease_1 = c("ALS"),
  selected_disease_2 = c("Control")
)

pwalk(params, run_Augur)

## log sessionInfo
sessionInfo()
