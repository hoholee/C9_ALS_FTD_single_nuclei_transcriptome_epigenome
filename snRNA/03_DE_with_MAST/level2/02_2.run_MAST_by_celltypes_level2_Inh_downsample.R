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
library(pheatmap)
library(scico)
library(viridis)
library(tictoc)
library(RColorBrewer)
library(shades)
library(ggpubr)
library(ggbeeswarm)
library(GGally)


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
disease_color_panel <- c("#f2a034", "#bf36ff", "#36bcff")
names(disease_color_panel) <- all_diseases

# threshold to keep genes that are expressed in certain percents of the chosen cells in comparison
freq_expressed <- 0.1
# fold-change threshold to call DE
# fc_threshold <- log2(1.1)

run_MAST_downsample <- function(selected_cluster, selected_region, selected_disease_1, selected_disease_2, selected_seed, n_cells_per_subject){

  message(paste0("Processing: ", selected_cluster, ", ", selected_region, ", ",
                 selected_disease_1, " vs ", selected_disease_2))
  selected_diseases <- c(selected_disease_1, selected_disease_2)
  excluded_disease <- all_diseases[!all_diseases %in% selected_diseases]

  # subset Seurat object
  tic("Subseting seurat object")
  sub_idx <- which(meta_data$rna_anno_2ndRound_level_2 == selected_cluster &
                     meta_data$region == selected_region &
                     meta_data$disease != excluded_disease)
  data_obj_subset <- data_obj_sub[, sub_idx]
  toc()

  # downsample using a fixed number of cells in each individual (`n_cells_per_subject`)
  # note that the selected cells will fallback to the number of cells in that individual if it's smaller than `n_cells_per_subject`

  meta_data_subset <- data_obj_subset@meta.data %>%
    as_tibble() %>%
    rowid_to_column("row_id")

  set.seed(selected_seed)
  meta_data_downsample <- meta_data_subset %>%
    group_by(subject) %>%
    slice_sample(n = n_cells_per_subject) %>%
    ungroup()

  n_selected_cells <- meta_data_downsample %>%
    dplyr::count(rna_anno_2ndRound_level_2, region, disease, subject) %>%
    mutate(seed = selected_seed)

  write_tsv(n_selected_cells, paste0("./results_downsample/num_selected_downsample_cells_per_subject_level2_",
                                     selected_cluster, "_", selected_region, "_",
                              selected_diseases[1], "_vs_", selected_diseases[2], "_seed_", selected_seed, ".txt"))

  idx_downsample <- meta_data_downsample$row_id %>% sort()

  # subset again to get the downsampled object
  tic("Downsampling seurat object")
  data_obj_used <- data_obj_subset[, idx_downsample]
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
  sca <- scaRaw[which(freq(scaRaw)>0), ]

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
  zlmCond <- zlm(~disease + cn_genes_on + cn_nCount_SCT + cn_percent_mt + cn_age + sex + seq_batch + (1|subject),
                 sca,
                 method = 'glmer',
                 ebayes = FALSE,
                 strictConvergence = FALSE,
                 fitArgsD = list(nAGQ = 0)
  )
  toc()

  saveRDS(zlmCond, paste0("./zlm_rds_downsample/MAST_zlmFit_level2_", selected_cluster, "_", selected_region, "_",
                          selected_diseases[1], "_vs_", selected_diseases[2], "_seed_", selected_seed, ".rds"))
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
  saveRDS(summary_cond, paste0("./zlm_rds_downsample/MAST_zlmSummary_level_2_", selected_cluster, "_", selected_region, "_",
                               selected_diseases[1], "_vs_", selected_diseases[2], "_seed_", selected_seed, ".rds"))

  tic("Cleaning up results and output")
  summary_Dt <- summary_cond$datatable

  fcHurdle <- merge(summary_Dt[contrast==lrt_term & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summary_Dt[contrast==lrt_term & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

  fcHurdle_df <- fcHurdle %>%
    as_tibble() %>%
    arrange(fdr) %>%
    dplyr::rename(gene = primerid,
                  p_value = `Pr(>Chisq)`,
                  model_log2FC = coef)

  # compute avg log2FC by hand
  data_cpm <- t(assay(sca)) %>%
    as_tibble() %>%
    mutate(disease = colData(sca)$disease, subject = colData(sca)$subject) %>%
    pivot_longer(names_to = "gene", values_to = "log2_CPM", cols = -c(disease, subject))
  data_fc <- data_cpm %>%
    group_by(gene, disease) %>%
    summarise(mean_log2_CPM = mean(log2_CPM)) %>%
    ungroup() %>%
    mutate(disease_group = case_when(disease == selected_disease_1 ~ "avg_CPM_disease_1",
                                     disease == selected_disease_2 ~ "avg_CPM_disease_2")) %>%
    dplyr::select(-disease) %>%
    pivot_wider(names_from = "disease_group", values_from = "mean_log2_CPM") %>%
    mutate(avg_log2FC = avg_CPM_disease_1 - avg_CPM_disease_2)

  final_res <- fcHurdle_df %>%
    left_join(zlmCond_conv) %>%
    left_join(data_fc)

  write_tsv(final_res, paste0("./results_downsample/MAST_res_level2_", selected_cluster, "_", selected_region, "_",
                              selected_diseases[1], "_vs_", selected_diseases[2], "_seed_", selected_seed, ".txt"))

  toc()
  rm(sca, scaRaw, data_obj_used, zlmCond, summary_cond, summary_Dt, data_cpm, data_fc, final_res, fcHurdle, fcHurdle_df)
}


param <- expand_grid(
  selected_cluster = c("Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST", "Inh_ADARB2_Other"),
  selected_region = c("MCX", "mFCX"),
  selected_seed = c(4, 8, 15, 16, 23, 42, 78, 173, 666, 2021))

pwalk(param, run_MAST_downsample,
      selected_disease_1 = "ALS", selected_disease_2 = "Control", n_cells_per_subject = 30)
pwalk(param, run_MAST_downsample,
      selected_disease_1 = "FTD", selected_disease_2 = "Control", n_cells_per_subject = 30)
pwalk(param, run_MAST_downsample,
      selected_disease_1 = "ALS", selected_disease_2 = "FTD", n_cells_per_subject = 30)

## log sessionInfo
sessionInfo()
