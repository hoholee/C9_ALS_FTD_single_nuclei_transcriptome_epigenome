library(tidyverse)
library(forcats)
library(pheatmap)
library(viridis)
library(scico)
library(RColorBrewer)
library(BuenColors)
library(qvalue)
library(dendsort)

res <- read_tsv("./ldsc_results/ldsc_res_ALS_FTD_latest_snATAC_anno_withOverlapAnnot_summary.txt")

group_order <- c(
  "Exc_superficial_MCX_ALS", "Exc_superficial_MCX_FTD", "Exc_superficial_MCX_Control",
  "Exc_superficial_mFCX_ALS", "Exc_superficial_mFCX_FTD", "Exc_superficial_mFCX_Control",
  "Exc_intermediate_MCX_ALS", "Exc_intermediate_MCX_FTD", "Exc_intermediate_MCX_Control",
  "Exc_intermediate_mFCX_ALS", "Exc_intermediate_mFCX_FTD", "Exc_intermediate_mFCX_Control",
  "Exc_deep_MCX_ALS", "Exc_deep_MCX_FTD", "Exc_deep_MCX_Control",
  "Exc_deep_mFCX_ALS", "Exc_deep_mFCX_FTD", "Exc_deep_mFCX_Control",
  "Inh_VIP_MCX_ALS", "Inh_VIP_MCX_FTD", "Inh_VIP_MCX_Control",
  "Inh_VIP_mFCX_ALS", "Inh_VIP_mFCX_FTD", "Inh_VIP_mFCX_Control",
  "Inh_LAMP5_MCX_ALS", "Inh_LAMP5_MCX_Control",
  "Inh_LAMP5_mFCX_ALS", "Inh_LAMP5_mFCX_FTD", "Inh_LAMP5_mFCX_Control",
  "Inh_PVALB_MCX_ALS", "Inh_PVALB_MCX_FTD", "Inh_PVALB_MCX_Control",
  "Inh_PVALB_mFCX_ALS", "Inh_PVALB_mFCX_FTD", "Inh_PVALB_mFCX_Control",
  "Inh_SST_MCX_ALS", "Inh_SST_MCX_FTD", "Inh_SST_MCX_Control",
  "Inh_SST_mFCX_ALS", "Inh_SST_mFCX_FTD", "Inh_SST_mFCX_Control",
  "Astro_MCX_ALS", "Astro_MCX_FTD", "Astro_MCX_Control",
  "Astro_mFCX_ALS", "Astro_mFCX_FTD", "Astro_mFCX_Control",
  "Micro_MCX_ALS", "Micro_MCX_FTD", "Micro_MCX_Control",
  "Micro_mFCX_ALS", "Micro_mFCX_FTD", "Micro_mFCX_Control",
  "Oligo_MCX_ALS", "Oligo_MCX_FTD", "Oligo_MCX_Control",
  "Oligo_mFCX_ALS", "Oligo_mFCX_FTD", "Oligo_mFCX_Control",
  "OPC_MCX_ALS", "OPC_MCX_FTD", "OPC_MCX_Control",
  "OPC_mFCX_ALS", "OPC_mFCX_FTD", "OPC_mFCX_Control"
)

res_mod <- res %>%
  group_by(group) %>%
  mutate(FDR = p.adjust(coefficient_p, method = "fdr")) %>%
  ungroup() %>%
  arrange(FDR)

mat <- res %>%
  mutate(enrichment = coefficient / coefficient_std_err) %>%
  dplyr::select(group, trait, enrichment) %>%
  spread(trait, enrichment) %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group) %>%
  column_to_rownames("group") %>%
  as.data.frame()

mat_fdr <- res_mod %>%
  dplyr::select(group, trait, FDR) %>%
  spread(trait, FDR) %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group) %>%
  column_to_rownames("group") %>%
  as.data.frame()

anno <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
colnames(anno) <- colnames(mat)
rownames(anno) <- rownames(mat)
# anno[mat_fdr < 0.1] <- "-"
anno[mat_fdr < 0.05] <- "*"
# anno[mat_fdr < 0.01] <- "**"
# anno[mat_fdr < 0.001] <- "***"

mat_p <- res_mod %>%
  dplyr::select(group, trait, coefficient_p) %>%
  spread(trait, coefficient_p) %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group) %>%
  column_to_rownames("group") %>%
  as.data.frame()

anno2 <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
colnames(anno2) <- colnames(mat)
rownames(anno2) <- rownames(mat)
# anno2[mat_p < 0.1] <- "-"
anno2[mat_p < 0.05] <- "*"
# anno2[mat_p < 0.01] <- "**"
# anno2[mat_p < 0.001] <- "***"

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
hclust_rows <- sort_hclust(hclust(dist(mat), method = "ward.D2"))
hclust_cols <- hclust(dist(t(mat)), method = "ward.D2")

pheatmap(
  mat,
  # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(80),
  color = scico(100, palette = "vik"),
  breaks = seq(-10, 10, length.out = 101),
  cellwidth = 8,
  cellheight = 8,
  fontsize_row = 8,
  fontsize_col = 8,
  display_numbers = anno,
  # display_numbers = anno2,
  fontsize_number = 6,
  cluster_cols = hclust_cols,
  # cluster_rows = hclust_rows,
  cluster_rows = FALSE,
  angle_col = 45,
  main = "snATAC peaks by group",
  filename = "ldsc_ALS_FTD_latest_snATAC_peaks_by_group_res_fdr.pdf"
)