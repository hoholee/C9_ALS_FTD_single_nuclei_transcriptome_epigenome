library(tidyverse)
library(viridis)
library(scico)
library(pheatmap)
library(plotly)
library(htmlwidgets)
library(ggrepel)
library(ggrastr)

# read DEG list
our_DE <- read_tsv("../MAST_res_level2_summary.txt", col_types = "cccciicdddddllddd")

Manolis_DE <- read_tsv("Manolis_ALS_single_cell_DE_genes_June2021_submission.tsv", col_types = "cccccddddddc")
Mathys_DE <- read_tsv("Mathys_Nat2019_AD_prefrontalCortex_single_cell_DE_noPathology_vs_Pathology.tsv", col_types = "ccddddddllc")
Morabito_DE <- read_tsv("Morabito_NatGenet2021_AD_prefrontalCortex_single_cell_cell_type_diagnosis_DE.tsv", col_types = "dddddcdcc")
Grubman_DE <- read_tsv("Grubman_Nat2019_AD_EntorhinalCortex_single_Cell_DE.tsv", col_types = "cddc")

# define DE threshold
FDR_threshold <- 0.05
FC_threshold <- 1.2

our_selected_cell_types <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

# subset significant DEGs in our DE list
our_DE_sig <- our_DE %>%
  filter(
    cell_type %in% our_selected_cell_types,
    FDR < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    abs(model_log2FC - avg_logFC) < 2
  )

# subset significant DEGs in published DE list using their reported threshold
# Manolis: p-adjust < 0.005, abs log2 FC > 1 standard deviation for each cell type
# (this SD range from 0.0236 to 0.1521, which is ~3-10% changes)
# set to a fixed values of 10% for all cell types for this comparison
Manolis_DE_sig <- Manolis_DE %>%
  filter(
    padj < 0.005,
    abs(log2FCMean) > log2(1.1)
  )

# Mathys: supported by both models using the criteria:
# FDR-corrected P < 0.01 in a two sided Wilcoxon-rank sum test,
# absolute log2(mean gene expression in AD category x/mean gene expression in AD category y) > 0.25,
# and FDR-corrected P < 0.05 in a Poisson mixed model
Mathys_DE_sig <- Mathys_DE %>%
  filter(
    DEGs.Ind.Model == TRUE,
    DEGs.Ind.Mix.models == TRUE
  )

# Morabito: MAST, already filtered by unknown threshold?
# appears to be abs log2 FC > 0.1, p-adjust < 0.05

# Grubman: edgeR, download from http://adsn.ddnetbio.com
# already filtered by FDR < 0.05, abs log2 FC > 1

n_all_genes <- 33538

# get Jaccard index for selected pair-wise comparisons
get_jaccard_Manolis <- function(
    selected_type_published,
    selected_comparison,
    selected_type_ours,
    selected_cond_1,
    selected_cond_2,
    selected_region) {
  published_DE_sub <- Manolis_DE_sig %>%
    filter(
      cell_type == selected_type_published,
      comparison == selected_comparison,
    ) %>%
    pull(Gene)

  our_DE_sub <- our_DE_sig %>%
    filter(
      cond_1 == selected_cond_1,
      cond_2 == selected_cond_2,
      cell_type == selected_type_ours,
      region == selected_region
    ) %>%
    pull(gene)

  tibble(
    type_published = selected_type_published,
    comparison = selected_comparison,
    type_ours = selected_type_ours,
    cond_1 = selected_cond_1,
    cond_2 = selected_cond_2,
    region = selected_region,
    n_DE_ours = length(our_DE_sub),
    n_DE_published = length(published_DE_sub),
    n_intersection = length(intersect(our_DE_sub, published_DE_sub)),
    n_union = length(union(our_DE_sub, published_DE_sub)),
    jaccard_index = n_intersection / n_union,
    p_hyper = sum(dhyper(n_intersection:n_DE_published, n_DE_ours, n_all_genes - n_DE_ours, n_DE_published)),
  )
}

# get_jaccard_Manolis(
#   selected_type_published = "OPC",
#   selected_comparison = "Manolis_C9ALS_vs_Control",
#   selected_type_ours = "OPC",
#   selected_cond_1 = "ALS",
#   selected_cond_2 = "Control",
#   selected_region = "MCX"
# )

params_Manolis <- expand_grid(
  selected_type_published = Manolis_DE_sig %>% distinct(cell_type) %>% pull(cell_type),
  selected_comparison = Manolis_DE_sig %>% distinct(comparison) %>% pull(comparison),
  selected_type_ours = our_selected_cell_types,
  selected_cond_1 = c("ALS", "FTD"),
  selected_cond_2 = "Control",
  selected_region = c("MCX", "mFCX")
)

res_Manolis <- pmap_df(params_Manolis, get_jaccard_Manolis)
write_tsv(res_Manolis, "DE_genes_overlaps_jaccard_Manolis.txt")


get_jaccard_Mathys <- function(
    selected_type_published,
    selected_comparison,
    selected_type_ours,
    selected_cond_1,
    selected_cond_2,
    selected_region) {
  published_DE_sub <- Mathys_DE_sig %>%
    filter(
      cell_type == selected_type_published,
      comparison == selected_comparison,
    ) %>%
    pull(gene)

  our_DE_sub <- our_DE_sig %>%
    filter(
      cond_1 == selected_cond_1,
      cond_2 == selected_cond_2,
      cell_type == selected_type_ours,
      region == selected_region
    ) %>%
    pull(gene)

  tibble(
    type_published = selected_type_published,
    comparison = selected_comparison,
    type_ours = selected_type_ours,
    cond_1 = selected_cond_1,
    cond_2 = selected_cond_2,
    region = selected_region,
    n_DE_ours = length(our_DE_sub),
    n_DE_published = length(published_DE_sub),
    n_intersection = length(intersect(our_DE_sub, published_DE_sub)),
    n_union = length(union(our_DE_sub, published_DE_sub)),
    jaccard_index = n_intersection / n_union,
    p_hyper = sum(dhyper(n_intersection:n_DE_published, n_DE_ours, n_all_genes - n_DE_ours, n_DE_published)),
  )
}

params_Mathys <- expand_grid(
  selected_type_published = Mathys_DE_sig %>% distinct(cell_type) %>% pull(cell_type),
  selected_comparison = Mathys_DE_sig %>% distinct(comparison) %>% pull(comparison),
  selected_type_ours = our_selected_cell_types,
  selected_cond_1 = c("ALS", "FTD"),
  selected_cond_2 = "Control",
  selected_region = c("MCX", "mFCX")
)

res_Mathys <- pmap_df(params_Mathys, get_jaccard_Mathys)
write_tsv(res_Mathys, "DE_genes_overlaps_jaccard_Mathys.txt")


get_jaccard_Morabito <- function(
    selected_type_published,
    selected_comparison,
    selected_type_ours,
    selected_cond_1,
    selected_cond_2,
    selected_region) {
  published_DE_sub <- Morabito_DE %>%
    filter(
      celltype == selected_type_published,
      comparison == selected_comparison,
    ) %>%
    pull(gene)

  our_DE_sub <- our_DE_sig %>%
    filter(
      cond_1 == selected_cond_1,
      cond_2 == selected_cond_2,
      cell_type == selected_type_ours,
      region == selected_region
    ) %>%
    pull(gene)

  tibble(
    type_published = selected_type_published,
    comparison = selected_comparison,
    type_ours = selected_type_ours,
    cond_1 = selected_cond_1,
    cond_2 = selected_cond_2,
    region = selected_region,
    n_DE_ours = length(our_DE_sub),
    n_DE_published = length(published_DE_sub),
    n_intersection = length(intersect(our_DE_sub, published_DE_sub)),
    n_union = length(union(our_DE_sub, published_DE_sub)),
    jaccard_index = n_intersection / n_union,
    p_hyper = sum(dhyper(n_intersection:n_DE_published, n_DE_ours, n_all_genes - n_DE_ours, n_DE_published)),
  )
}

params_Morabito <- expand_grid(
  selected_type_published = Morabito_DE %>% distinct(celltype) %>% pull(celltype),
  selected_comparison = Morabito_DE %>% distinct(comparison) %>% pull(comparison),
  selected_type_ours = our_selected_cell_types,
  selected_cond_1 = c("ALS", "FTD"),
  selected_cond_2 = "Control",
  selected_region = c("MCX", "mFCX")
)

res_Morabito <- pmap_df(params_Morabito, get_jaccard_Morabito)
write_tsv(res_Morabito, "DE_genes_overlaps_jaccard_Morabito.txt")


get_jaccard_Grubman <- function(
    selected_type_published,
    selected_comparison,
    selected_type_ours,
    selected_cond_1,
    selected_cond_2,
    selected_region) {
  published_DE_sub <- Grubman_DE %>%
    filter(
      cell_type == selected_type_published,
    ) %>%
    pull(geneName)

  our_DE_sub <- our_DE_sig %>%
    filter(
      cond_1 == selected_cond_1,
      cond_2 == selected_cond_2,
      cell_type == selected_type_ours,
      region == selected_region
    ) %>%
    pull(gene)

  tibble(
    type_published = selected_type_published,
    comparison = selected_comparison,
    type_ours = selected_type_ours,
    cond_1 = selected_cond_1,
    cond_2 = selected_cond_2,
    region = selected_region,
    n_DE_ours = length(our_DE_sub),
    n_DE_published = length(published_DE_sub),
    n_intersection = length(intersect(our_DE_sub, published_DE_sub)),
    n_union = length(union(our_DE_sub, published_DE_sub)),
    jaccard_index = n_intersection / n_union,
    p_hyper = sum(dhyper(n_intersection:n_DE_published, n_DE_ours, n_all_genes - n_DE_ours, n_DE_published)),
  )
}

params_Grubman <- expand_grid(
  selected_type_published = Grubman_DE %>% distinct(cell_type) %>% pull(cell_type),
  selected_comparison = "Grubman_Nat2019_AD_EntorhinalCortex",
  selected_type_ours = our_selected_cell_types,
  selected_cond_1 = c("ALS", "FTD"),
  selected_cond_2 = "Control",
  selected_region = c("MCX", "mFCX")
)

res_Grubman <- pmap_df(params_Grubman, get_jaccard_Grubman)
write_tsv(res_Grubman, "DE_genes_overlaps_jaccard_Grubman.txt")

# plot heatmap of Jaccard index
res_sig_Mathys <- res_Mathys %>%
  # filter(p_hyper < 0.05) %>%
  mutate(
    group_published = paste(comparison, type_published, sep = "_"),
    group_ours = paste(region, type_ours, cond_1, cond_2, sep = "_")
  ) %>%
  select(group_published, group_ours, jaccard_index)

res_mat_Mathys <- res_sig_Mathys %>%
  pivot_wider(names_from = "group_ours", values_from = "jaccard_index", values_fill = NA_real_) %>%
  as.data.frame() %>%
  column_to_rownames("group_published")

pheatmap(res_mat_Mathys,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  # color = colorRampPalette(scico(100, palette = "acton", direction = -1))(100),
  # color = colorRampPalette(viridis(100, option = "A", end = 0.97))(100),
  border_color = NA,
  cellwidth = 5.5,
  cellheight = 5.5,
  fontsize = 5.5,
  filename = "DE_overlaps_with_Mathys_jaccard.pdf",
  width = 6, height = 4
)

# remove neurons in FTD
# custom order of rows and columns
res_sig_Mathys <- res_Mathys %>%
  # filter(p_hyper < 0.05) %>%
  filter(!(cond_1 == "FTD" & grepl("Exc_|Inh_", type_ours))) %>%
  mutate(
    group_published = paste(comparison, type_published, sep = "_"),
    group_ours = paste(region, type_ours, cond_1, cond_2, sep = "_")
  ) %>%
  select(group_published, group_ours, jaccard_index) %>%
  mutate(
    group_published = str_replace(group_published, "Mathys_Nat2019_AD_prefrontalCortex_noPathology_vs_Pathology_", ""),
    group_published = factor(group_published, levels = c("Ast", "Mic", "Oli", "Opc", "Ex", "In")),
    region = if_else(
      grepl("MCX_", group_ours),
      "Motor",
      "Frontal"
    ),
    disease = if_else(
      grepl("_ALS_", group_ours),
      "ALS",
      "FTD"
    ),
    cell_type = group_ours,
    cell_type = str_replace(cell_type, "MCX_|mFCX_", ""),
    cell_type = str_replace(cell_type, "_ALS_Control|_FTD_Control", ""),
    cell_type = factor(
      cell_type,
      levels = c(
        "Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC",
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
      )
    ),
    disease = factor(disease, levels = c("ALS", "FTD")),
    region = factor(region, levels = c("Motor", "Frontal"))
  ) %>%
  arrange(group_published, cell_type, disease, region)

res_mat_Mathys <- res_sig_Mathys %>%
  select(group_published, group_ours, jaccard_index) %>%
  pivot_wider(names_from = "group_ours", values_from = "jaccard_index", values_fill = NA_real_) %>%
  as.data.frame() %>%
  column_to_rownames("group_published")

column_annotations <- res_sig_Mathys %>%
  distinct(group_ours, cell_type, region, disease) %>%
  as.data.frame() %>%
  column_to_rownames("group_ours")

column_annotations_colors <- list(
  disease = c(
    "ALS" = "#F6416C",
    "FTD" = "#FFDE7D"
  ),
  region = c(
    "Motor" = "#432266",
    "Frontal" = "#FAA51B"
  ),
  cell_type = c(
    "Astro" = "#ef6075",
    "Micro" = "#ec5b42",
    "Oligo" = "#89520d",
    "OPC" = "#b98735",
    "Endo" = "#ad312c",
    "VLMC" = "#7c5816",
    "Exc_superficial" = "#54c879",
    "Exc_intermediate" = "#bfe676",
    "Exc_deep" = "#c9a836",
    "Inh_VIP" = "#cb0077",
    "Inh_LAMP5" = "#016cf2",
    "Inh_ADARB2_Other" = "#02a7b8",
    "Inh_PVALB" = "#8e90d8",
    "Inh_SST" = "#721297"
  )
)

pheatmap(res_mat_Mathys,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  # color = colorRampPalette(scico(100, palette = "acton", direction = -1))(100),
  # color = colorRampPalette(viridis(100, option = "A", end = 0.97))(100),
  border_color = NA,
  cellwidth = 5.5,
  cellheight = 5.5,
  fontsize = 5.5,
  annotation_col = column_annotations,
  annotation_colors = column_annotations_colors,
  filename = "DE_overlaps_with_Mathys_jaccard_rmFTDneurons.pdf",
  width = 6, height = 4
)


res_sig_Morabito <- res_Morabito %>%
  # filter(p_hyper < 0.05) %>%
  mutate(
    group_published = paste(comparison, type_published, sep = "_"),
    group_ours = paste(region, type_ours, cond_1, cond_2, sep = "_")
  ) %>%
  select(group_published, group_ours, jaccard_index)

res_mat_Morabito <- res_sig_Morabito %>%
  pivot_wider(names_from = "group_ours", values_from = "jaccard_index", values_fill = NA_real_) %>%
  as.data.frame() %>%
  column_to_rownames("group_published")

pheatmap(res_mat_Morabito,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  # color = colorRampPalette(scico(100, palette = "acton", direction = -1))(100),
  # color = colorRampPalette(viridis(100, option = "A", end = 0.97))(100),
  border_color = NA,
  cellwidth = 5.5,
  cellheight = 5.5,
  fontsize = 5.5,
  filename = "DE_overlaps_with_Morabito_jaccard.pdf",
  width = 6, height = 4
)

# remove neurons in FTD
# custom order of rows and columns
res_sig_Morabito <- res_Morabito %>%
  # filter(p_hyper < 0.05) %>%
  filter(!(cond_1 == "FTD" & grepl("Exc_|Inh_", type_ours))) %>%
  mutate(
    group_published = paste(comparison, type_published, sep = "_"),
    group_ours = paste(region, type_ours, cond_1, cond_2, sep = "_")
  ) %>%
  select(group_published, group_ours, jaccard_index) %>%
  mutate(
    group_published = str_replace(group_published, "Morabito_NatGenet2021_AD_prefrontalCortex_cell_type_diagnosis_DE_", ""),
    group_published = factor(group_published, levels = c("ASC", "MG", "ODC", "OPC", "PER.END", "EX", "INH")),
    region = if_else(
      grepl("MCX_", group_ours),
      "Motor",
      "Frontal"
    ),
    disease = if_else(
      grepl("_ALS_", group_ours),
      "ALS",
      "FTD"
    ),
    cell_type = group_ours,
    cell_type = str_replace(cell_type, "MCX_|mFCX_", ""),
    cell_type = str_replace(cell_type, "_ALS_Control|_FTD_Control", ""),
    cell_type = factor(
      cell_type,
      levels = c(
        "Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC",
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
      )
    ),
    disease = factor(disease, levels = c("ALS", "FTD")),
    region = factor(region, levels = c("Motor", "Frontal"))
  ) %>%
  arrange(group_published, cell_type, disease, region)

res_mat_Morabito <- res_sig_Morabito %>%
  select(group_published, group_ours, jaccard_index) %>%
  pivot_wider(names_from = "group_ours", values_from = "jaccard_index", values_fill = NA_real_) %>%
  as.data.frame() %>%
  column_to_rownames("group_published")

column_annotations <- res_sig_Morabito %>%
  distinct(group_ours, cell_type, region, disease) %>%
  as.data.frame() %>%
  column_to_rownames("group_ours")

column_annotations_colors <- list(
  disease = c(
    "ALS" = "#F6416C",
    "FTD" = "#FFDE7D"
  ),
  region = c(
    "Motor" = "#432266",
    "Frontal" = "#FAA51B"
  ),
  cell_type = c(
    "Astro" = "#ef6075",
    "Micro" = "#ec5b42",
    "Oligo" = "#89520d",
    "OPC" = "#b98735",
    "Endo" = "#ad312c",
    "VLMC" = "#7c5816",
    "Exc_superficial" = "#54c879",
    "Exc_intermediate" = "#bfe676",
    "Exc_deep" = "#c9a836",
    "Inh_VIP" = "#cb0077",
    "Inh_LAMP5" = "#016cf2",
    "Inh_ADARB2_Other" = "#02a7b8",
    "Inh_PVALB" = "#8e90d8",
    "Inh_SST" = "#721297"
  )
)

pheatmap(res_mat_Morabito,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  # color = colorRampPalette(scico(100, palette = "acton", direction = -1))(100),
  # color = colorRampPalette(viridis(100, option = "A", end = 0.97))(100),
  border_color = NA,
  cellwidth = 5.5,
  cellheight = 5.5,
  fontsize = 5.5,
  annotation_col = column_annotations,
  annotation_colors = column_annotations_colors,
  filename = "DE_overlaps_with_Morabito_jaccard_rmFTDneurons.pdf",
  width = 6, height = 4
)

# remove all FTD
# custom order of rows and columns
res_sig_Morabito <- res_Morabito %>%
  filter(cond_1 != "FTD") %>%
  mutate(
    group_published = paste(comparison, type_published, sep = "_"),
    group_ours = paste(region, type_ours, cond_1, cond_2, sep = "_")
  ) %>%
  select(group_published, group_ours, jaccard_index) %>%
  mutate(
    group_published = str_replace(group_published, "Morabito_NatGenet2021_AD_prefrontalCortex_cell_type_diagnosis_DE_", ""),
    group_published = factor(group_published, levels = c("ASC", "MG", "ODC", "OPC", "PER.END", "EX", "INH")),
    region = if_else(
      grepl("MCX_", group_ours),
      "Motor",
      "Frontal"
    ),
    disease = if_else(
      grepl("_ALS_", group_ours),
      "ALS",
      "FTD"
    ),
    cell_type = group_ours,
    cell_type = str_replace(cell_type, "MCX_|mFCX_", ""),
    cell_type = str_replace(cell_type, "_ALS_Control|_FTD_Control", ""),
    cell_type = factor(
      cell_type,
      levels = c(
        "Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC",
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
      )
    ),
    disease = factor(disease, levels = c("ALS", "FTD")),
    region = factor(region, levels = c("Motor", "Frontal"))
  ) %>%
  arrange(group_published, cell_type, disease, region)

res_mat_Morabito <- res_sig_Morabito %>%
  select(group_published, group_ours, jaccard_index) %>%
  pivot_wider(names_from = "group_ours", values_from = "jaccard_index", values_fill = NA_real_) %>%
  as.data.frame() %>%
  column_to_rownames("group_published")

column_annotations <- res_sig_Morabito %>%
  distinct(group_ours, cell_type, region) %>%
  as.data.frame() %>%
  column_to_rownames("group_ours")

column_annotations_colors <- list(
  region = c(
    "Motor" = "#432266",
    "Frontal" = "#FAA51B"
  ),
  cell_type = c(
    "Astro" = "#ef6075",
    "Micro" = "#ec5b42",
    "Oligo" = "#89520d",
    "OPC" = "#b98735",
    "Endo" = "#ad312c",
    "VLMC" = "#7c5816",
    "Exc_superficial" = "#54c879",
    "Exc_intermediate" = "#bfe676",
    "Exc_deep" = "#c9a836",
    "Inh_VIP" = "#cb0077",
    "Inh_LAMP5" = "#016cf2",
    "Inh_ADARB2_Other" = "#02a7b8",
    "Inh_PVALB" = "#8e90d8",
    "Inh_SST" = "#721297"
  )
)

pheatmap(res_mat_Morabito,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  # color = colorRampPalette(scico(100, palette = "acton", direction = -1))(100),
  # color = colorRampPalette(viridis(100, option = "A", end = 0.97))(100),
  border_color = NA,
  cellwidth = 5.5,
  cellheight = 5.5,
  fontsize = 5.5,
  annotation_col = column_annotations,
  annotation_colors = column_annotations_colors,
  filename = "DE_overlaps_with_Morabito_jaccard_rmFTD.pdf",
  width = 6, height = 4
)


res_sig_Grubman <- res_Grubman %>%
  # filter(p_hyper < 0.05) %>%
  mutate(
    group_published = paste(comparison, type_published, sep = "_"),
    group_ours = paste(region, type_ours, cond_1, cond_2, sep = "_")
  ) %>%
  select(group_published, group_ours, jaccard_index)

res_mat_Grubman <- res_sig_Grubman %>%
  pivot_wider(names_from = "group_ours", values_from = "jaccard_index", values_fill = NA_real_) %>%
  as.data.frame() %>%
  column_to_rownames("group_published")

pheatmap(res_mat_Grubman,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  # color = colorRampPalette(scico(100, palette = "acton", direction = -1))(100),
  # color = colorRampPalette(viridis(100, option = "A", end = 0.97))(100),
  border_color = NA,
  cellwidth = 5.5,
  cellheight = 5.5,
  fontsize = 5.5,
  filename = "DE_overlaps_with_Grubman_jaccard.pdf",
  width = 6, height = 4
)


res_sig_Manolis <- res_Manolis %>%
  # filter(p_hyper < 0.05) %>%
  mutate(
    group_published = paste(comparison, type_published, sep = "_"),
    group_ours = paste(region, type_ours, cond_1, cond_2, sep = "_")
  ) %>%
  select(group_published, group_ours, jaccard_index)

res_mat_Manolis <- res_sig_Manolis %>%
  pivot_wider(names_from = "group_ours", values_from = "jaccard_index", values_fill = NA_real_) %>%
  as.data.frame() %>%
  column_to_rownames("group_published")

pheatmap(res_mat_Manolis,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  # color = colorRampPalette(scico(100, palette = "acton", direction = -1))(100),
  # color = colorRampPalette(viridis(100, option = "A", end = 0.97))(100),
  border_color = NA,
  cellwidth = 5.5,
  cellheight = 5.5,
  fontsize = 5.5,
  filename = "DE_overlaps_with_Manolis_jaccard.pdf",
  width = 6, height = 4
)

# get Pearson/Spearman correlation of log2FC between for our DE genes and published DE genes
# Morabito
plot_FC_cor_Morabito <- function(snRNA_cell_type, Morabito_cell_type,
                                 snRNA_cond_1, snRNA_cond_2,
                                 snRNA_region) {
  message(str_glue(
    "Processing snRNA {snRNA_cell_type} vs. Morabito {Morabito_cell_type} in {snRNA_region} ",
    "with FC > {FC_threshold} and FDR < {FDR_threshold}, {snRNA_cond_1} vs. {snRNA_cond_2}"
  ))
  snRNA_sub <- our_DE_sig %>%
    filter(
      cell_type == snRNA_cell_type,
      cond_1 == snRNA_cond_1,
      cond_2 == snRNA_cond_2,
      region == snRNA_region
    )
  if (nrow(snRNA_sub) == 0) {
    message(str_glue("No significant DE genes in snRNA {snRNA_cell_type}!"))
  } else {
    Morabito_DE_sub <- Morabito_DE %>%
      filter(
        celltype == Morabito_cell_type
      )
    res <- snRNA_sub %>%
      full_join(Morabito_DE_sub, by = c("gene" = "gene"))
    write_tsv(res, str_glue("./plot_vs_Morabito/{snRNA_cond_1}_vs_{snRNA_cond_2}_snRNA_{snRNA_cell_type}_expFC_vs_Morabito_{Morabito_cell_type}_in_{snRNA_region}_FC_{FC_threshold}_FDR_{FDR_threshold}_merged_table.tsv"))
    cor_spearman <- cor.test(res$model_log2FC, res$avg_logFC.y, method = "spearman")
    cor_pearson <- cor.test(res$model_log2FC, res$avg_logFC.y, method = "pearson")
    num_overlap_genes <- res %>%
      drop_na() %>%
      nrow()
    rho <- round(cor_spearman$estimate, 3)
    p_value_spearman <- signif(cor_spearman$p.value, 3)
    pearson_r <- round(cor_pearson$estimate, 3)
    p_value_pearson <- signif(cor_pearson$p.value, 3)
    out_stat_df <- tibble(
      snRNA_cond_1 = snRNA_cond_1,
      snRNA_cond_2 = snRNA_cond_2,
      snRNA_cell_type = snRNA_cell_type,
      Morabito_cell_type = Morabito_cell_type,
      region = snRNA_region,
      FC_threshold = FC_threshold,
      FDR_threshold = FDR_threshold,
      num_overlap_genes = num_overlap_genes,
      spearman_rho = cor_spearman$estimate,
      spearman_p_value = cor_spearman$p.value,
      pearson_r = cor_pearson$estimate,
      pearson_p_value = cor_pearson$p.value
    )
    write_tsv(out_stat_df, str_glue("./plot_vs_Morabito/{snRNA_cond_1}_vs_{snRNA_cond_2}_snRNA_{snRNA_cell_type}_expFC_vs_Morabito_{Morabito_cell_type}_in_{snRNA_region}_FC_{FC_threshold}_FDR_{FDR_threshold}_stat.txt"))
    p <- res %>% ggplot(aes(model_log2FC, avg_logFC.y))
    p1 <- p +
      geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
      geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
      geom_point(size = 0.1) +
      xlab(str_glue("snRNA gene expresssion log2FC\n({snRNA_cond_1} vs. {snRNA_cond_2})")) +
      ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
      ggtitle(str_glue(
        "snRNA {snRNA_cell_type} vs. Morabito {Morabito_cell_type} in {snRNA_region}\n",
        "FC > {FC_threshold}, FDR < {FDR_threshold}\n",
        "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}, r = {pearson_r}, p = {p_value_pearson}"
      )) +
      theme_bw(base_size = 8, base_family = "Helvetica") +
      theme(panel.grid.minor = element_blank())
    ggsave(str_glue("./plot_vs_Morabito/{snRNA_cond_1}_vs_{snRNA_cond_2}_snRNA_{snRNA_cell_type}_expFC_vs_Morabito_{Morabito_cell_type}_in_{snRNA_region}_FC_{FC_threshold}_FDR_{FDR_threshold}.pdf"),
      plot = p1, device = cairo_pdf(), width = 4, height = 3, useDingbats = FALSE
    )
    dev.off()
    p2 <- p +
      geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
      geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
      geom_point(size = 0.1, aes(label = gene)) +
      xlab(str_glue("snRNA gene expresssion log2FC\n({snRNA_cond_1} vs. {snRNA_cond_2})")) +
      ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
      ggtitle(str_glue(
        "snRNA {snRNA_cell_type} vs. Morabito {Morabito_cell_type} in {snRNA_region}\n",
        "FC > {FC_threshold}, FDR < {FDR_threshold}\n",
        "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}, r = {pearson_r}, p = {p_value_pearson}"
      )) +
      theme_bw(base_size = 8, base_family = "Helvetica") +
      theme(panel.grid.minor = element_blank())
    p2 <- ggplotly(p2)
    saveWidget(
      p2,
      str_glue("./plot_vs_Morabito/{snRNA_cond_1}_vs_{snRNA_cond_2}_snRNA_{snRNA_cell_type}_expFC_vs_Morabito_{Morabito_cell_type}_in_{snRNA_region}_FC_{FC_threshold}_FDR_{FDR_threshold}.html"),
      selfcontained = TRUE
    )
  }
}

plot_FC_cor_Morabito("Astro", "ASC", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Astro", "ASC", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Astro", "ASC", "FTD", "Control", "MCX")
plot_FC_cor_Morabito("Astro", "ASC", "FTD", "Control", "mFCX")

plot_FC_cor_Morabito("Oligo", "ODC", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Oligo", "ODC", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Oligo", "ODC", "FTD", "Control", "MCX")
plot_FC_cor_Morabito("Oligo", "ODC", "FTD", "Control", "mFCX")

plot_FC_cor_Morabito("Micro", "MG", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Micro", "MG", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Micro", "MG", "FTD", "Control", "MCX")
plot_FC_cor_Morabito("Micro", "MG", "FTD", "Control", "mFCX")

plot_FC_cor_Morabito("OPC", "OPC", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("OPC", "OPC", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("OPC", "OPC", "FTD", "Control", "MCX")
plot_FC_cor_Morabito("OPC", "OPC", "FTD", "Control", "mFCX")

plot_FC_cor_Morabito("Endo", "PER.END", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Endo", "PER.END", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Endo", "PER.END", "FTD", "Control", "MCX")
plot_FC_cor_Morabito("Endo", "PER.END", "FTD", "Control", "mFCX")

plot_FC_cor_Morabito("VLMC", "PER.END", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("VLMC", "PER.END", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("VLMC", "PER.END", "FTD", "Control", "MCX")
plot_FC_cor_Morabito("VLMC", "PER.END", "FTD", "Control", "mFCX")

plot_FC_cor_Morabito("Exc_superficial", "EX", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Exc_superficial", "EX", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Exc_intermediate", "EX", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Exc_intermediate", "EX", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Exc_deep", "EX", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Exc_deep", "EX", "ALS", "Control", "mFCX")

plot_FC_cor_Morabito("Inh_PVALB", "INH", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Inh_PVALB", "INH", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Inh_SST", "INH", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Inh_SST", "INH", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Inh_VIP", "INH", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Inh_VIP", "INH", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Inh_LAMP5", "INH", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Inh_LAMP5", "INH", "ALS", "Control", "mFCX")
plot_FC_cor_Morabito("Inh_ADARB2_Other", "INH", "ALS", "Control", "MCX")
plot_FC_cor_Morabito("Inh_ADARB2_Other", "INH", "ALS", "Control", "mFCX")

# read summarized FC correlation results
cor_summary <- read_tsv("./our_snRNA_vs_Morabito_FC_correlation_stat_summary.txt")
region_relabel <- c(
  "MCX" = "motor cortex",
  "mFCX" = "frontal cortex"
)

num_overlap_genes_cutoff <- 5

plot_FC_correlation_heatmap <- function(selected_disease,
                                        selected_region,
                                        selected_FC_threshold) {
  df_sub <- cor_summary %>%
    filter(
      FC_threshold == selected_FC_threshold,
      snRNA_cond_1 == selected_disease,
      region == selected_region,
      num_overlap_genes > num_overlap_genes_cutoff
    )

  if (selected_disease == "ALS") {
    rho_mat <- df_sub %>%
      select(snRNA_cell_type, Morabito_cell_type, spearman_rho) %>%
      pivot_wider(
        names_from = "Morabito_cell_type",
        values_from = "spearman_rho"
      ) %>%
      mutate(snRNA_cell_type = factor(snRNA_cell_type,
        levels = c(
          "Astro", "Micro", "Oligo", "OPC",
          "Exc_superficial", "Exc_intermediate", "Exc_deep",
          "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
        )
      )) %>%
      arrange(snRNA_cell_type) %>%
      select(snRNA_cell_type, ASC, MG, ODC, OPC, EX, INH) %>%
      as.data.frame() %>%
      column_to_rownames(var = "snRNA_cell_type")

    sig_mat <- df_sub %>%
      mutate(
        rho_FDR = p.adjust(spearman_p_value, method = "fdr"),
        r_FDR = p.adjust(pearson_p_value, method = "fdr"),
        rho_sig = case_when(
          rho_FDR < 0.05 ~ "*",
          rho_FDR < 0.01 ~ "**",
          rho_FDR < 0.001 ~ "***",
          TRUE ~ ""
        ),
        r_sig = case_when(
          r_FDR < 0.05 ~ "*",
          r_FDR < 0.01 ~ "**",
          r_FDR < 0.001 ~ "***",
          TRUE ~ ""
        )
      ) %>%
      select(snRNA_cell_type, Morabito_cell_type, rho_sig) %>%
      pivot_wider(
        names_from = "Morabito_cell_type",
        values_from = "rho_sig",
        values_fill = ""
      ) %>%
      mutate(snRNA_cell_type = factor(snRNA_cell_type,
        levels = c(
          "Astro", "Micro", "Oligo", "OPC",
          "Exc_superficial", "Exc_intermediate", "Exc_deep",
          "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
        )
      )) %>%
      arrange(snRNA_cell_type) %>%
      select(snRNA_cell_type, ASC, MG, ODC, OPC, EX, INH) %>%
      as.data.frame() %>%
      column_to_rownames(var = "snRNA_cell_type")

    r_mat <- df_sub %>%
      select(snRNA_cell_type, Morabito_cell_type, pearson_r) %>%
      pivot_wider(
        names_from = "Morabito_cell_type",
        values_from = "pearson_r"
      ) %>%
      mutate(snRNA_cell_type = factor(snRNA_cell_type,
        levels = c(
          "Astro", "Micro", "Oligo", "OPC",
          "Exc_superficial", "Exc_intermediate", "Exc_deep",
          "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
        )
      )) %>%
      arrange(snRNA_cell_type) %>%
      select(snRNA_cell_type, ASC, MG, ODC, OPC, EX, INH) %>%
      as.data.frame() %>%
      column_to_rownames(var = "snRNA_cell_type")

    r_sig_mat <- df_sub %>%
      mutate(
        rho_FDR = p.adjust(spearman_p_value, method = "fdr"),
        r_FDR = p.adjust(pearson_p_value, method = "fdr"),
        rho_sig = case_when(
          rho_FDR < 0.05 ~ "*",
          rho_FDR < 0.01 ~ "**",
          rho_FDR < 0.001 ~ "***",
          TRUE ~ ""
        ),
        r_sig = case_when(
          r_FDR < 0.05 ~ "*",
          r_FDR < 0.01 ~ "**",
          r_FDR < 0.001 ~ "***",
          TRUE ~ ""
        )
      ) %>%
      select(snRNA_cell_type, Morabito_cell_type, r_sig) %>%
      pivot_wider(
        names_from = "Morabito_cell_type",
        values_from = "r_sig",
        values_fill = ""
      ) %>%
      mutate(snRNA_cell_type = factor(snRNA_cell_type,
        levels = c(
          "Astro", "Micro", "Oligo", "OPC",
          "Exc_superficial", "Exc_intermediate", "Exc_deep",
          "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
        )
      )) %>%
      arrange(snRNA_cell_type) %>%
      select(snRNA_cell_type, ASC, MG, ODC, OPC, EX, INH) %>%
      as.data.frame() %>%
      column_to_rownames(var = "snRNA_cell_type")
  } else {
    rho_mat <- df_sub %>%
      select(snRNA_cell_type, Morabito_cell_type, spearman_rho) %>%
      pivot_wider(
        names_from = "Morabito_cell_type",
        values_from = "spearman_rho"
      ) %>%
      mutate(snRNA_cell_type = factor(snRNA_cell_type,
        levels = c(
          "Astro", "Micro", "Oligo", "OPC"
        )
      )) %>%
      arrange(snRNA_cell_type) %>%
      select(snRNA_cell_type, ASC, MG, ODC, OPC) %>%
      as.data.frame() %>%
      column_to_rownames(var = "snRNA_cell_type")

    sig_mat <- df_sub %>%
      mutate(
        rho_FDR = p.adjust(spearman_p_value, method = "fdr"),
        r_FDR = p.adjust(pearson_p_value, method = "fdr"),
        rho_sig = case_when(
          rho_FDR < 0.05 ~ "*",
          rho_FDR < 0.01 ~ "**",
          rho_FDR < 0.001 ~ "***",
          TRUE ~ ""
        ),
        r_sig = case_when(
          r_FDR < 0.05 ~ "*",
          r_FDR < 0.01 ~ "**",
          r_FDR < 0.001 ~ "***",
          TRUE ~ ""
        )
      ) %>%
      select(snRNA_cell_type, Morabito_cell_type, rho_sig) %>%
      pivot_wider(
        names_from = "Morabito_cell_type",
        values_from = "rho_sig",
        values_fill = ""
      ) %>%
      mutate(snRNA_cell_type = factor(snRNA_cell_type,
        levels = c(
          "Astro", "Micro", "Oligo", "OPC"
        )
      )) %>%
      arrange(snRNA_cell_type) %>%
      select(snRNA_cell_type, ASC, MG, ODC, OPC) %>%
      as.data.frame() %>%
      column_to_rownames(var = "snRNA_cell_type")

    r_mat <- df_sub %>%
      select(snRNA_cell_type, Morabito_cell_type, pearson_r) %>%
      pivot_wider(
        names_from = "Morabito_cell_type",
        values_from = "pearson_r"
      ) %>%
      mutate(snRNA_cell_type = factor(snRNA_cell_type,
        levels = c(
          "Astro", "Micro", "Oligo", "OPC"
        )
      )) %>%
      arrange(snRNA_cell_type) %>%
      select(snRNA_cell_type, ASC, MG, ODC, OPC) %>%
      as.data.frame() %>%
      column_to_rownames(var = "snRNA_cell_type")

    r_sig_mat <- df_sub %>%
      mutate(
        rho_FDR = p.adjust(spearman_p_value, method = "fdr"),
        r_FDR = p.adjust(pearson_p_value, method = "fdr"),
        rho_sig = case_when(
          rho_FDR < 0.05 ~ "*",
          rho_FDR < 0.01 ~ "**",
          rho_FDR < 0.001 ~ "***",
          TRUE ~ ""
        ),
        r_sig = case_when(
          r_FDR < 0.05 ~ "*",
          r_FDR < 0.01 ~ "**",
          r_FDR < 0.001 ~ "***",
          TRUE ~ ""
        )
      ) %>%
      select(snRNA_cell_type, Morabito_cell_type, r_sig) %>%
      pivot_wider(
        names_from = "Morabito_cell_type",
        values_from = "r_sig",
        values_fill = ""
      ) %>%
      mutate(snRNA_cell_type = factor(snRNA_cell_type,
        levels = c(
          "Astro", "Micro", "Oligo", "OPC"
        )
      )) %>%
      arrange(snRNA_cell_type) %>%
      select(snRNA_cell_type, ASC, MG, ODC, OPC) %>%
      as.data.frame() %>%
      column_to_rownames(var = "snRNA_cell_type")
  }
  pheatmap(rho_mat,
    color = scico(100, palette = "vikO", direction = 1),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 8,
    cellheight = 8,
    fontsize = 8,
    show_rownames = TRUE,
    show_colnames = TRUE,
    border_color = NA,
    scale = "none",
    display_numbers = sig_mat,
    number_color = "#000000",
    main = str_glue(
      "Spearman correlation of our_snRNA {selected_disease} vs. Control {region_relabel[selected_region]} \n",
      "vs. Morabito AD vs. Control expression FC\n",
      "in significnat DE genes (FDR < 0.05, FC > {selected_FC_threshold})"
    ),
    filename = str_glue(
      "./plot_vs_Morabito/our_snRNA_vs_Morabito_AD_FC_spearman_cor_heatmap_",
      "{selected_disease}_vs_Control_{selected_region}_FC_{selected_FC_threshold}.pdf"
    ),
    width = 6,
    height = 4
  )
  pheatmap(r_mat,
    color = scico(100, palette = "vikO", direction = 1),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 8,
    cellheight = 8,
    fontsize = 8,
    show_rownames = TRUE,
    show_colnames = TRUE,
    border_color = NA,
    scale = "none",
    display_numbers = r_sig_mat,
    number_color = "#000000",
    main = str_glue(
      "Pearson correlation of our_snRNA {selected_disease} vs. Control {region_relabel[selected_region]} \n",
      "vs. Morabito AD vs. Control expression FC\n",
      "in significnat DE genes (FDR < 0.05, FC > {selected_FC_threshold})"
    ),
    filename = str_glue(
      "./plot_vs_Morabito/our_snRNA_vs_Morabito_AD_FC_pearson_cor_heatmap_",
      "{selected_disease}_vs_Control_{selected_region}_FC_{selected_FC_threshold}.pdf"
    ),
    width = 6,
    height = 4
  )
}

plot_FC_correlation_heatmap("ALS", "MCX", 1.2)
plot_FC_correlation_heatmap("ALS", "mFCX", 1.2)
plot_FC_correlation_heatmap("FTD", "MCX", 1.2)
plot_FC_correlation_heatmap("FTD", "mFCX", 1.2)


# scatter plots of FC correalation for selected cell types, highlighting example genes
# ALS vs. Control, Astro, MCX
df <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Astro_expFC_vs_Morabito_ASC_in_MCX_FC_1.2_FDR_0.05_merged_table.tsv")
test_stats <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Astro_expFC_vs_Morabito_ASC_in_MCX_FC_1.2_FDR_0.05_stat.txt")

# df %>% drop_na() %>% filter(model_log2FC < 0, avg_logFC.y < 0) %>% pull(gene)
highlighted_genes <- c(
  "COLEC12", "MAOB", "FAS", "NEAT1", "CD44",
  "GRIA2", "FOXG1",
  "APOE",
  "RANBP3L",
  "PCLO", "LINC00499"
)

# only keep genes that are singificant DE in both datasets
df_sub <- df %>% drop_na()
num_overlap_genes <- test_stats$num_overlap_genes
rho <- test_stats$spearman_rho
p_value_spearman <- test_stats$spearman_p_value

p <- df_sub %>% ggplot(aes(model_log2FC, avg_logFC.y))
p +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
  geom_point_rast(
    color = "#ef6075",
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    fill = "#ef6075",
    color = "black",
    data = df_sub %>% filter(gene %in% highlighted_genes),
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% highlighted_genes),
    size = 2
  ) +
  xlab(str_glue("snRNA gene expresssion log2FC\n(ALS vs. Control)")) +
  ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
  ggtitle(
    str_glue(
      "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}"
    )
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())

ggsave(
  "./plot_vs_Morabito/ALS_vs_Control_snRNA_Astro_expFC_vs_Morabito_ASC_in_MCX_FC_1.2_FDR_0.05_add_example_genes.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 2,
  useDingbats = FALSE
)
dev.off()

# ALS vs. Control, Astro, mFCX
df <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Astro_expFC_vs_Morabito_ASC_in_mFCX_FC_1.2_FDR_0.05_merged_table.tsv")
test_stats <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Astro_expFC_vs_Morabito_ASC_in_mFCX_FC_1.2_FDR_0.05_stat.txt")

# df %>% drop_na() %>% filter(model_log2FC < 0, avg_logFC.y < 0) %>% pull(gene)
highlighted_genes <- c(
  "CHI3L1", "MAOB", "FAS", "NEAT1", "CD44",
  "GRIA2", "TET1",
  "KCNIP1", "KNDC1",
  "NEBL", "NFIB",
  "PCLO", "LINC00499"
)

# only keep genes that are singificant DE in both datasets
df_sub <- df %>% drop_na()
num_overlap_genes <- test_stats$num_overlap_genes
rho <- test_stats$spearman_rho
p_value_spearman <- test_stats$spearman_p_value

p <- df_sub %>% ggplot(aes(model_log2FC, avg_logFC.y))
p +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
  geom_point_rast(
    color = "#ef6075",
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    fill = "#ef6075",
    color = "black",
    data = df_sub %>% filter(gene %in% highlighted_genes),
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% highlighted_genes),
    size = 2
  ) +
  xlab(str_glue("snRNA gene expresssion log2FC\n(ALS vs. Control)")) +
  ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
  ggtitle(str_glue(
    "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}"
  )) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())

ggsave(
  "./plot_vs_Morabito/ALS_vs_Control_snRNA_Astro_expFC_vs_Morabito_ASC_in_mFCX_FC_1.2_FDR_0.05_add_example_genes.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 2,
  useDingbats = FALSE
)
dev.off()

# ALS vs. Control, Oligo, MCX
df <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Oligo_expFC_vs_Morabito_ODC_in_MCX_FC_1.2_FDR_0.05_merged_table.tsv")
test_stats <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Oligo_expFC_vs_Morabito_ODC_in_MCX_FC_1.2_FDR_0.05_stat.txt")

# df %>% drop_na() %>% filter(model_log2FC < 0, avg_logFC.y < 0) %>% pull(gene)
highlighted_genes <- c(
  "PLP1", "UNC5C", "NRXN3",
  "FXR1", "STAT3",
  "COL18A1", "RANBP17", "RASGRF2", "ARHGAP22",
  "ENO4", "MAP4K5", "SLC38A2",
  "FCHSD2", "LAMA2", "NRXN3", "XIST", "NEAT1", "NAV2"
)

# only keep genes that are singificant DE in both datasets
df_sub <- df %>% drop_na()
num_overlap_genes <- test_stats$num_overlap_genes
rho <- test_stats$spearman_rho
p_value_spearman <- test_stats$spearman_p_value

p <- df_sub %>% ggplot(aes(model_log2FC, avg_logFC.y))
p +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
  geom_point_rast(
    color = "#89520d",
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    fill = "#89520d",
    color = "black",
    data = df_sub %>% filter(gene %in% highlighted_genes),
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% highlighted_genes),
    size = 2
  ) +
  xlab(str_glue("snRNA gene expresssion log2FC\n(ALS vs. Control)")) +
  ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
  ggtitle(
    str_glue(
      "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}"
    )
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())

ggsave(
  "./plot_vs_Morabito/ALS_vs_Control_snRNA_Oligo_expFC_vs_Morabito_ODC_in_MCX_FC_1.2_FDR_0.05_add_example_genes.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 2,
  useDingbats = FALSE
)
dev.off()

# ALS vs. Control, Oligo, mFCX
df <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Oligo_expFC_vs_Morabito_ODC_in_mFCX_FC_1.2_FDR_0.05_merged_table.tsv")
test_stats <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Oligo_expFC_vs_Morabito_ODC_in_mFCX_FC_1.2_FDR_0.05_stat.txt")

# df %>% drop_na() %>% filter(model_log2FC < 0, avg_logFC.y < 0) %>% pull(gene)
highlighted_genes <- c(
  "PLP1", "UNC5C", "NRXN3", "CTNND2", "DNM1L",
  "FXR1", "STAT3", "HDAC4", "SEMA5A", "FOXO1", "HSPA1A",
  "COL18A1", "RANBP17", "RASGRF2", "ARHGAP22",
  "ENO4", "MAP4K5", "SLC38A2",
  "FCHSD2", "LAMA2", "NRXN3", "XIST", "NEAT1", "NAV2", "LINC01099"
)

# only keep genes that are singificant DE in both datasets
df_sub <- df %>% drop_na()
num_overlap_genes <- test_stats$num_overlap_genes
rho <- test_stats$spearman_rho
p_value_spearman <- test_stats$spearman_p_value

p <- df_sub %>% ggplot(aes(model_log2FC, avg_logFC.y))
p +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
  geom_point_rast(
    color = "#89520d",
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    fill = "#89520d",
    color = "black",
    data = df_sub %>% filter(gene %in% highlighted_genes),
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% highlighted_genes),
    size = 2
  ) +
  xlab(str_glue("snRNA gene expresssion log2FC\n(ALS vs. Control)")) +
  ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
  ggtitle(
    str_glue(
      "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}"
    )
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())

ggsave(
  "./plot_vs_Morabito/ALS_vs_Control_snRNA_Oligo_expFC_vs_Morabito_ODC_in_mFCX_FC_1.2_FDR_0.05_add_example_genes.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 2,
  useDingbats = FALSE
)
dev.off()

# ALS vs. Control, Exc upper, MCX
df <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_superficial_expFC_vs_Morabito_EX_in_MCX_FC_1.2_FDR_0.05_merged_table.tsv")
test_stats <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_superficial_expFC_vs_Morabito_EX_in_MCX_FC_1.2_FDR_0.05_stat.txt")

# df %>% drop_na() %>% filter(model_log2FC < 0, avg_logFC.y < 0) %>% pull(gene)
highlighted_genes <- c(
  "TET3", "RBFOX3", "DNMT3A", "SLIT1", "MEIS3",
  "SLIT2", "FXR1", "HSPA1A", "HSPH1", "SNAP23",
  "GRID2", "ROBO2", "MOXD1",
  "KIF1A", "NTRK2", "HDAC4",
  "SNX31"
)

# only keep genes that are singificant DE in both datasets
df_sub <- df %>% drop_na()
num_overlap_genes <- test_stats$num_overlap_genes
rho <- test_stats$spearman_rho
p_value_spearman <- test_stats$spearman_p_value

p <- df_sub %>% ggplot(aes(model_log2FC, avg_logFC.y))
p +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
  geom_point_rast(
    color = "#54c879",
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    fill = "#54c879",
    color = "black",
    data = df_sub %>% filter(gene %in% highlighted_genes),
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% highlighted_genes),
    size = 2,
    max.overlaps = 15
  ) +
  xlab(str_glue("snRNA gene expresssion log2FC\n(ALS vs. Control)")) +
  ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
  ggtitle(
    str_glue(
      "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}"
    )
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())

ggsave(
  "./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_superficial_expFC_vs_Morabito_EX_in_MCX_FC_1.2_FDR_0.05_add_example_genes.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 2,
  useDingbats = FALSE
)
dev.off()

# ALS vs. Control, Exc upper, mFCX
df <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_superficial_expFC_vs_Morabito_EX_in_mFCX_FC_1.2_FDR_0.05_merged_table.tsv")
test_stats <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_superficial_expFC_vs_Morabito_EX_in_mFCX_FC_1.2_FDR_0.05_stat.txt")

# df %>% drop_na() %>% filter(model_log2FC < 0, avg_logFC.y < 0) %>% pull(gene)
highlighted_genes <- c(
  "TET3", "RBFOX3", "DNMT3A", "BACH2", "BCL2", "SMAD3", "SGK1",
  "SLC26A3", "FXR1", "HSPD1", "HSPH1", "LAMA2",
  "LINC02296", "TMSB4X", "SPARCL1",
  "PDE1A", "CHD20", "TRPC6"
)

# only keep genes that are singificant DE in both datasets
df_sub <- df %>% drop_na()
num_overlap_genes <- test_stats$num_overlap_genes
rho <- test_stats$spearman_rho
p_value_spearman <- test_stats$spearman_p_value

p <- df_sub %>% ggplot(aes(model_log2FC, avg_logFC.y))
p +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
  geom_point_rast(
    color = "#54c879",
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    fill = "#54c879",
    color = "black",
    data = df_sub %>% filter(gene %in% highlighted_genes),
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% highlighted_genes),
    size = 2,
    max.overlaps = 15
  ) +
  xlab(str_glue("snRNA gene expresssion log2FC\n(ALS vs. Control)")) +
  ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
  ggtitle(
    str_glue(
      "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}"
    )
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())

ggsave(
  "./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_superficial_expFC_vs_Morabito_EX_in_mFCX_FC_1.2_FDR_0.05_add_example_genes.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 2,
  useDingbats = FALSE
)
dev.off()

# ALS vs. Control, Exc deep, MCX
df <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_deep_expFC_vs_Morabito_EX_in_MCX_FC_1.2_FDR_0.05_merged_table.tsv")
test_stats <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_deep_expFC_vs_Morabito_EX_in_MCX_FC_1.2_FDR_0.05_stat.txt")

# df %>% drop_na() %>% filter(model_log2FC < 0, avg_logFC.y < 0) %>% pull(gene)
highlighted_genes <- c(
  "CAMK2A", "POLR1A", "DNMT3A", "NAV1", "PIGL", "SLC26A3", "RIT2",
  "HSPH1", "RSPO2", "SNX31",
  "CALM3", "DLGAP1-AS4",
  "ZDHHC11B", "SPHKAP", "PTPRT"
)

# only keep genes that are singificant DE in both datasets
df_sub <- df %>% drop_na()
num_overlap_genes <- test_stats$num_overlap_genes
rho <- test_stats$spearman_rho
p_value_spearman <- test_stats$spearman_p_value

p <- df_sub %>% ggplot(aes(model_log2FC, avg_logFC.y))
p +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
  geom_point_rast(
    color = "#c9a836",
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    fill = "#c9a836",
    color = "black",
    data = df_sub %>% filter(gene %in% highlighted_genes),
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% highlighted_genes),
    size = 2,
    max.overlaps = 15
  ) +
  xlab(str_glue("snRNA gene expresssion log2FC\n(ALS vs. Control)")) +
  ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
  ggtitle(
    str_glue(
      "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}"
    )
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())

ggsave(
  "./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_deep_expFC_vs_Morabito_EX_in_MCX_FC_1.2_FDR_0.05_add_example_genes.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 2,
  useDingbats = FALSE
)
dev.off()

# ALS vs. Control, Exc deep, mFCX
df <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_deep_expFC_vs_Morabito_EX_in_mFCX_FC_1.2_FDR_0.05_merged_table.tsv")
test_stats <- read_tsv("./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_deep_expFC_vs_Morabito_EX_in_mFCX_FC_1.2_FDR_0.05_stat.txt")

# df %>% drop_na() %>% filter(model_log2FC < 0, avg_logFC.y < 0) %>% pull(gene)
highlighted_genes <- c(
  "DNMT3A", "TET3", "RBFOX3", "SLC6A17", "SEPT9",
  "HSPH1", "SNX31", "SLC26A3", "TTN", "SNAP23", "HSPA1A",
  "TMSB4X", "SREBF2", "TMEM165",
  "PTPRT", "HTR1E"
)

# only keep genes that are singificant DE in both datasets
df_sub <- df %>% drop_na()
num_overlap_genes <- test_stats$num_overlap_genes
rho <- test_stats$spearman_rho
p_value_spearman <- test_stats$spearman_p_value

p <- df_sub %>% ggplot(aes(model_log2FC, avg_logFC.y))
p +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
  geom_point_rast(
    color = "#c9a836",
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    fill = "#c9a836",
    color = "black",
    data = df_sub %>% filter(gene %in% highlighted_genes),
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% highlighted_genes),
    size = 2,
    max.overlaps = 15
  ) +
  xlab(str_glue("snRNA gene expresssion log2FC\n(ALS vs. Control)")) +
  ylab(str_glue("Morabito gene expression log2FC\n(AD vs. Control)")) +
  ggtitle(
    str_glue(
      "n = {num_overlap_genes}, rho = {rho}, p = {p_value_spearman}"
    )
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())

ggsave(
  "./plot_vs_Morabito/ALS_vs_Control_snRNA_Exc_deep_expFC_vs_Morabito_EX_in_mFCX_FC_1.2_FDR_0.05_add_example_genes.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 2,
  useDingbats = FALSE
)
dev.off()
