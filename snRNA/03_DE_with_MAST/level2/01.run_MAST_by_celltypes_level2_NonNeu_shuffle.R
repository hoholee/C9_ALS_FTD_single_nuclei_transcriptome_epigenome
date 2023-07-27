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
  "ALS_vs_Control" = list("half" = combn(indi_ALS, 3, simplify = TRUE)[, c(1:10)],
                          "full" = combn(indi_Control, 3, simplify = TRUE)),
  "FTD_vs_Control" = list("half" = combn(indi_FTD, 3, simplify = TRUE)[, c(1:10)],
                          "full" = combn(indi_Control, 3, simplify = TRUE)),
  "ALS_vs_FTD" = list("half" = combn(indi_ALS, 3, simplify = TRUE)[, c(1:10)],
                      "full" = combn(indi_FTD, 3, simplify = TRUE))
)

combin_idx <- expand.grid(idx_1 = c(1:10), idx_2 = c(1:20)) %>% as_tibble()

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

run_MAST <- function(selected_cluster, selected_region, selected_disease_1, selected_disease_2){

  message(paste0("Processing: ", selected_cluster, ", ", selected_region, ", ",
                 selected_disease_1, " vs ", selected_disease_2))
  selected_diseases <- c(selected_disease_1, selected_disease_2)
  excluded_disease <- all_diseases[!all_diseases %in% selected_diseases]

  # subset Seurat object
  tic("Subseting seurat object")
  sub_idx <- which(meta_data$rna_anno_2ndRound_level_2 == selected_cluster &
                     meta_data$region == selected_region &
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

  # random shuffle of the labels of each subject
  
  meta <- colData(sca) %>% as_tibble()
  subject_disease_match <- meta %>% distinct(disease, subject)
  
  run_shuffle <- function(selected_seed){
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
     
    names(subject_disease_match_permuted$disease_permute) <- subject_disease_match_permuted$subject
    
    colData(sca)$disease <- subject_disease_match_permuted$disease_permute[colData(sca)$subject]
    
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

  saveRDS(zlmCond, paste0("./zlm_rds_shuffle/MAST_zlmFit_level2_", selected_cluster, "_", selected_region, "_",
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
  saveRDS(summary_cond, paste0("./zlm_rds_shuffle/MAST_zlmSummary_level_2_", selected_cluster, "_", selected_region, "_",
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

    fcHurdle_df %>% 
      left_join(zlmCond_conv) %>% 
      left_join(data_fc) %>% 
      mutate(seed = selected_seed,
             permuted_subject_disease_1 = str_flatten(c(group_A_1, group_A_2), collapse = ","),
             permuted_subject_disease_2 = str_flatten(c(group_B_1, group_B_2), collapse = ","),
             )
  }

  final_res <- map_dfr(1:200, run_shuffle)
  

  write_tsv(final_res, paste0("./results_shuffle/MAST_res_level2_", selected_cluster, "_", selected_region, "_",
                              selected_diseases[1], "_vs_", selected_diseases[2], ".txt"))

  toc()
  rm(sca, scaRaw, data_obj_used, zlmCond, summary_cond, summary_Dt, data_cpm, data_fc, final_res, fcHurdle, fcHurdle_df)
}


param <- expand_grid(
  selected_cluster = c("Astro", "Endo", "Micro", "OPC", "Oligo", "VLMC"),
  selected_region = c("MCX", "mFCX"))

walk2(param$selected_cluster, param$selected_region, run_MAST,
      selected_disease_1 = "ALS", selected_disease_2 = "Control")
# walk2(param$selected_cluster, param$selected_region, run_MAST,
#       selected_disease_1 = "FTD", selected_disease_2 = "Control")
# walk2(param$selected_cluster, param$selected_region, run_MAST,
#       selected_disease_1 = "ALS", selected_disease_2 = "FTD")

## log sessionInfo
sessionInfo()
