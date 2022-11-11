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

run_MAST <- function(selected_cluster){

  message(paste0("Processing: ", selected_cluster))

  # subset Seurat object
  # comparing MCX vs. mFCX in Control samples
  tic("Subseting seurat object")
  sub_idx <- which(meta_data$rna_anno_2ndRound_level_2 == selected_cluster &
                     meta_data$disease == "Control")
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
  cond <- factor(colData(sca)$region)
  # use frontal cortex (mFCX) as base level
  cond <- relevel(cond, "mFCX")
  colData(sca)$region <- cond

  # fitting with a mixed model
  tic("Fitting ZLM")
  zlmCond <- zlm(~region + cn_genes_on + cn_nCount_SCT + cn_percent_mt + cn_age + sex + seq_batch + (1|subject),
                 sca,
                 method = 'glmer',
                 ebayes = FALSE,
                 strictConvergence = FALSE,
                 fitArgsD = list(nAGQ = 0)
  )
  toc()

  saveRDS(zlmCond, paste0("./zlm_rds_MCX_vs_mFCX/MAST_zlmFit_level2_", selected_cluster, "_Control_MCX_vs_mFCX.rds"))
  # save convergence
  zlmCond_conv <- zlmCond@converged %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    dplyr::rename(conv_C = C, conv_D = D)

  lrt_term <- "regionMCX"

  tic("Calculating LRT")
  summary_cond <- summary(zlmCond, doLRT = lrt_term, fitArgsD = list(nAGQ = 0))
  toc()
  saveRDS(summary_cond, paste0("./zlm_rds_MCX_vs_mFCX/MAST_zlmSummary_level_2_", selected_cluster, "_Control_MCX_vs_mFCX.rds"))

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
    mutate(region = colData(sca)$region, subject = colData(sca)$subject) %>%
    pivot_longer(names_to = "gene", values_to = "log2_CPM", cols = -c(region, subject))
  data_fc <- data_cpm %>%
    group_by(gene, region) %>%
    summarise(mean_log2_CPM = mean(log2_CPM)) %>%
    ungroup() %>%
    pivot_wider(names_from = "region", values_from = "mean_log2_CPM") %>%
    dplyr::rename(avg_CPM_MCX = MCX, avg_CPM_mFCX = mFCX) %>%
    mutate(avg_log2FC = avg_CPM_MCX - avg_CPM_mFCX)

  final_res <- fcHurdle_df %>%
    left_join(zlmCond_conv) %>%
    left_join(data_fc)

  write_tsv(final_res, paste0("./results_MCX_vs_mFCX/MAST_res_level2_", selected_cluster, "_Control_MCX_vs_mFCX.txt"))

  toc()
  rm(sca, scaRaw, data_obj_used, zlmCond, summary_cond, summary_Dt, data_cpm, data_fc, final_res, fcHurdle, fcHurdle_df)
}


cluster_list = c("Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST", "Inh_ADARB2_Other")

walk(cluster_list, run_MAST)

## log sessionInfo
sessionInfo()
