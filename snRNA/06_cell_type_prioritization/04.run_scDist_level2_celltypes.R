# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(scDist)

data_obj_sub <- readRDS("../seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_SeuratV4_object.rds")

meta_data <- data_obj_sub@meta.data %>% as_tibble()

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

  # run scDist
  set.seed(666)
  res <- scDist(data_used_mat, data_cell_meta,
    fixed.effects = c(
      "disease",
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
      "scDist_results_SCT_scaled_",
      "{selected_region}_",
      "{selected_disease_1}_vs_{selected_disease_2}.rds"
    )
  )

  # save results as CSV
  res_df <- res$results %>%
    rownames_to_column("cell_type") %>%
    as_tibble()
  write_tsv(
    res_df,
    str_glue(
      "scDist_results_SCT_scaled_",
      "{selected_region}_",
      "{selected_disease_1}_vs_{selected_disease_2}.tsv"
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
      "scDist_results_SCT_scaled_",
      "{selected_region}_",
      "{selected_disease_1}_vs_{selected_disease_2}.pdf"
    ),
    plot = p,
    device = cairo_pdf(),
    width = 4,
    height = 3,
    useDingbats = FALSE
  )
}


# run scDist in groups
params <- expand_grid(
  selected_region = c("MCX", "mFCX"),
  selected_disease_1 = c("ALS", "FTD"),
  selected_disease_2 = c("Control")
)

pwalk(params, run_scDist)

## log sessionInfo
sessionInfo()
