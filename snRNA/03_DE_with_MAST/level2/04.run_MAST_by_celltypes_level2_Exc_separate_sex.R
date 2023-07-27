# run this dobro w/o activated conda env, with `MAST` v1.18.0 installed
# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(data.table)
library(NMF)
library(rsvd)
library(lme4)
library(Seurat)
library(sctransform)
library(future)
library(scales)
library(MAST)
library(magrittr)
library(tictoc)


options(mc.cores = 8)
# read in annotated data
data_obj_sub <- readRDS("../../seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_SeuratV4_object.rds")

meta_data <- data_obj_sub@meta.data %>% as_tibble()

# convert the Seurat object to the SummarizedExperiment class
# subset by cluster and by region
# using the "RNA" assay
# and use the raw counts in the `counts` slot
selected_assay <- "RNA"
selected_slot <- "counts"

all_diseases <- c("ALS", "FTD", "Control")

# threshold to keep genes that are expressed in certain percents of the chosen cells in comparison
freq_expressed <- 0.1
# fold-change threshold to call DE
# fc_threshold <- log2(1.1)

# test in female and male separately
run_MAST <- function(selected_cluster, selected_region, selected_disease_1, selected_disease_2, selected_sex) {
  message(str_glue("Processing: {selected_cluster}, {selected_region}, {selected_sex}, {selected_disease_1} vs {selected_disease_2}"))
  selected_diseases <- c(selected_disease_1, selected_disease_2)
  excluded_disease <- all_diseases[!all_diseases %in% selected_diseases]

  # subset Seurat object
  tic("Subseting seurat object")
  sub_idx <- which(meta_data$rna_anno_2ndRound_level_2 == selected_cluster &
    meta_data$region == selected_region &
    meta_data$sex == selected_sex &
    meta_data$disease != excluded_disease)
  data_obj_used <- data_obj_sub[, sub_idx]
  toc()

  tic("Extracting counts and normalized")
  # extract counts
  data_used <- GetAssayData(data_obj_used[[selected_assay]], slot = selected_slot)
  # convert counts into log2(CPM+1)
  data_used@x <- log2((10^6 * data_used@x / rep.int(colSums(data_used), diff(data_used@p))) + 1)
  toc()

  tic("Createing sca object")
  data_cell_meta <- data_obj_used@meta.data
  data_feature_meta <- tibble(gene_name = row.names(data_used))

  scaRaw <- FromMatrix(
    exprsArray = Matrix::as.matrix(data_used),
    cData = data_cell_meta,
    fData = data_feature_meta
  )

  # Remove genes are not expressed in any of the chosen cells
  sca <- scaRaw[which(freq(scaRaw) > 0), ]

  # calculate cellular detection rate (nGenesOn)
  # cdr <-colSums(assay(sca)>0)
  # this turns out to be excatly the same as the `nFeature_SCT` column...
  # all.equal(unname(cdr), colData(sca)$nFeature_SCT)
  # qplot(x=cdr, y=colData(sca)$nFeature_SCT) + xlab('CDR') + ylab('nFeature_SCT')

  # scale every covariant in consideration
  colData(sca)$cn_genes_on <- scale(colData(sca)$nFeature_SCT)
  colData(sca)$cn_age <- scale(colData(sca)$age)
  colData(sca)$cn_nCount_SCT <- scale(colData(sca)$nCount_SCT)
  colData(sca)$cn_percent_mt <- scale(colData(sca)$percent_mt)

  # only keep genes that found in at least `freq_expressed` of the chosen cells
  expressed_genes <- freq(sca) > freq_expressed
  sca <- sca[expressed_genes, ]
  toc()

  # differential expression using a Hurdle model
  cond <- factor(colData(sca)$disease)
  # use the 2nd condition in `selected_diseases` as base level
  cond <- relevel(cond, selected_diseases[2])
  colData(sca)$disease <- cond

  # fitting with a mixed model
  tic("Fitting ZLM")
  zlmCond <- zlm(~ disease + cn_genes_on + cn_nCount_SCT + cn_percent_mt + cn_age + seq_batch + (1 | subject),
    sca,
    method = "glmer",
    ebayes = FALSE,
    strictConvergence = FALSE,
    fitArgsD = list(nAGQ = 0)
  )
  toc()

  saveRDS(zlmCond, paste0(
    "./zlm_rds_separate_sex/MAST_zlmFit_level2_", selected_cluster, "_", selected_region, "_", selected_sex, "_",
    selected_diseases[1], "_vs_", selected_diseases[2], ".rds"
  ))
  # save convergence
  zlmCond_conv <- zlmCond@converged %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    dplyr::rename(conv_C = C, conv_D = D)

  lrt_term <- paste0("disease", selected_diseases[1])

  tic("Calculating LRT")
  summary_cond <- summary(zlmCond, doLRT = lrt_term, fitArgsD = list(nAGQ = 0))
  toc()
  saveRDS(summary_cond, paste0(
    "./zlm_rds_separate_sex/MAST_zlmSummary_level_2_", selected_cluster, "_", selected_region, "_", selected_sex, "_",
    selected_diseases[1], "_vs_", selected_diseases[2], ".rds"
  ))

  tic("Cleaning up results and output")
  summary_Dt <- summary_cond$datatable

  fcHurdle <- merge(summary_Dt[contrast == lrt_term & component == "H", .(primerid, `Pr(>Chisq)`)], # hurdle P values
    summary_Dt[contrast == lrt_term & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
    by = "primerid"
  ) # logFC coefficients
  fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, "fdr")]

  fcHurdle_df <- fcHurdle %>%
    as_tibble() %>%
    arrange(fdr) %>%
    dplyr::rename(
      gene = primerid,
      p_value = `Pr(>Chisq)`,
      model_log2FC = coef
    )

  # compute avg log2FC by hand
  data_cpm <- t(assay(sca)) %>%
    as_tibble() %>%
    mutate(disease = colData(sca)$disease, subject = colData(sca)$subject) %>%
    pivot_longer(names_to = "gene", values_to = "log2_CPM", cols = -c(disease, subject))
  data_fc <- data_cpm %>%
    group_by(gene, disease) %>%
    summarise(mean_log2_CPM = mean(log2_CPM)) %>%
    ungroup() %>%
    mutate(disease_group = case_when(
      disease == selected_disease_1 ~ "avg_CPM_disease_1",
      disease == selected_disease_2 ~ "avg_CPM_disease_2"
    )) %>%
    dplyr::select(-disease) %>%
    pivot_wider(names_from = "disease_group", values_from = "mean_log2_CPM") %>%
    mutate(avg_log2FC = avg_CPM_disease_1 - avg_CPM_disease_2)

  final_res <- fcHurdle_df %>%
    left_join(zlmCond_conv) %>%
    left_join(data_fc)

  write_tsv(final_res, paste0(
    "./results_separate_sex/MAST_res_level2_", selected_cluster, "_", selected_region, "_", selected_sex, "_",
    selected_diseases[1], "_vs_", selected_diseases[2], ".txt"
  ))

  toc()
  rm(sca, scaRaw, data_obj_used, zlmCond, summary_cond, summary_Dt, data_cpm, data_fc, final_res, fcHurdle, fcHurdle_df)
}

# test run
# run_MAST("Inh_VIP", "MCX", "ALS", "Control")

param <- expand_grid(
  selected_cluster = c("Exc_superficial", "Exc_intermediate", "Exc_deep"),
  selected_region = c("MCX", "mFCX"),
  selected_sex = c("M", "F"),
  selected_disease_1 = c("ALS", "FTD"),
  selected_disease_2 = c("Control", "FTD")
) %>%
  filter(selected_disease_1 != selected_disease_2)

pwalk(param, run_MAST)

## log sessionInfo
sessionInfo()
