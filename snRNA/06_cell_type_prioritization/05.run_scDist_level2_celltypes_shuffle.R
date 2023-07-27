# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(scDist)

data_obj_sub <- readRDS("../seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_SeuratV4_object.rds")

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

# use the "RNA" assay and the raw counts in the `counts` slot
# selected_assay <- "RNA"
# selected_slot <- "counts"

# use the "SCT" assay and the residuals from the negative binomial regression in scTransform
# as recommended by the scDist authors
# this "normalized" data is saved in the `scale.data` slot
selected_assay <- "SCT"
selected_slot <- "scale.data"

all_diseases <- c("ALS", "FTD", "Control")
excluded_cell_types <- c("NK_cell", "Exc_unknown", "Ambiguous")

run_scDist <- function(selected_region, selected_disease_1, selected_disease_2) {
  message(
    str_glue(
      "Running scDist for {selected_region} {selected_disease_1} vs {selected_disease_2}..."
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

  # extract counts
  data_used <- GetAssayData(data_obj_used[[selected_assay]], slot = selected_slot)

  # the following is not needed if we use the `SCT` assay and the `scale.data` slot
  # convert counts into log2(CPM+1)
  # data_used@x <- log2((10^6 * data_used@x / rep.int(colSums(data_used), diff(data_used@p))) + 1)

  # convert to full matrix
  data_used_mat <- as.matrix(data_used)

  # get the meta data in the subset cells
  data_cell_meta <- data_obj_used@meta.data

  # scale continuous variables
  # data_cell_meta$cn_nFeature_RNA <- scale(data_cell_meta$nFeature_RNA)
  # data_cell_meta$cn_nCount_RNA <- scale(data_cell_meta$nCount_RNA)
  data_cell_meta$cn_nFeature_SCT <- scale(data_cell_meta$nFeature_SCT)
  data_cell_meta$cn_nCount_SCT <- scale(data_cell_meta$nCount_SCT)
  data_cell_meta$cn_percent_mt <- scale(data_cell_meta$percent_mt)
  data_cell_meta$cn_age <- scale(data_cell_meta$age)

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

    data_cell_meta$disease_permute <- disease_permute_dict[data_cell_meta$subject]

    # run scDist
    set.seed(666)
    res <- scDist(data_used_mat, data_cell_meta,
      fixed.effects = c(
        "disease_permute",
        # "cn_nFeature_RNA", "cn_nCount_RNA",
        "cn_nFeature_SCT", "cn_nCount_SCT",
        "cn_percent_mt", "cn_age",
        "sex", "seq_batch"
      ),
      random.effects = "subject",
      clusters = "rna_anno_2ndRound_level_2",
      d = 20,
      truncate = FALSE,
      min.counts.per.cell = 20
    )

    # save results as RDS
    saveRDS(
      res,
      str_glue(
        "./scDist_shuffle/scDist_results_SCT_scaled_",
        "{selected_region}_",
        "{selected_disease_1}_vs_{selected_disease_2}_",
        "permuted_seed_{selected_seed}.rds"
      )
    )

    # save results as CSV
    res_df <- res$results %>%
      rownames_to_column("cell_type") %>%
      as_tibble() %>%
      mutate(
        seed = selected_seed,
        permuted_subject_disease_1 = str_flatten(c(group_A_1, group_A_2), collapse = ","),
        permuted_subject_disease_2 = str_flatten(c(group_B_1, group_B_2), collapse = ","),
      )
    write_tsv(
      res_df,
      str_glue(
        "./scDist_shuffle/scDist_results_SCT_scaled_",
        "{selected_region}_",
        "{selected_disease_1}_vs_{selected_disease_2}_",
        "permuted_seed_{selected_seed}.tsv"
      )
    )

    # plot results
    p <- DistPlot(res, return.plot = TRUE) +
      ylab(
        str_glue(
          "{selected_disease_1} vs {selected_disease_2} difference ",
          "in {selected_region} (scDist distance)"
        )
      ) +
      xlab("Cell type") +
      theme_bw(base_size = 8, base_family = "Helvetica") +
      theme(panel.grid.minor = element_blank())

    ggsave(
      filename = str_glue(
        "./scDist_shuffle/scDist_results_SCT_scaled_",
        "{selected_region}_",
        "{selected_disease_1}_vs_{selected_disease_2}_",
        "permuted_seed_{selected_seed}.pdf"
      ),
      plot = p,
      device = cairo_pdf(),
      width = 4,
      height = 3,
      useDingbats = FALSE
    )
    dev.off()
  }

  walk(seed_list, run_shuffle)
}

# run scDist in groups
params <- expand_grid(
  selected_region = c("MCX", "mFCX"),
  selected_disease_1 = c("ALS"),
  selected_disease_2 = c("Control")
)

pwalk(params, run_scDist)

## log sessionInfo
sessionInfo()
