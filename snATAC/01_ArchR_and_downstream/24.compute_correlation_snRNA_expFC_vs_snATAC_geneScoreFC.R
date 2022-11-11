# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(pheatmap)
library(viridis)
library(scico)

rna <- read_tsv("MAST_res_level2_summary.txt")

get_FC_cor <- function(selected_cell_type_rna,
                       selected_cell_type_atac,
                       selected_region,
                       selected_cond_1,
                       selected_cond_2,
                       FDR_threshold,
                       FC_threshold) {
  message(
    paste0(
      "Processing ", selected_cell_type_rna, " vs. ", selected_cell_type_atac,
      " in ", selected_region, " ",
      selected_cond_1, " vs. ", selected_cond_2,
      " with FC: ", FC_threshold, " and FDR: ", FDR_threshold
    )
  )
  rna_sub <- rna %>%
    filter(
      cell_type == selected_cell_type_rna,
      cond_1 == selected_cond_1,
      cond_2 == selected_cond_2,
      region == selected_region,
      FDR < FDR_threshold,
      abs(model_log2FC) > log2(FC_threshold),
      conv_C == TRUE,
      conv_D == TRUE,
      model_log2FC_ci_hi * model_log2FC_ci_low > 0,
      abs(model_log2FC - avg_logFC) < 2
    )
  if (nrow(rna_sub) == 0) {
    message("No significant DE genes in snRNA!")
  } else {
    num_sig_DE <- nrow(rna_sub)
    atac <- read_tsv(
      paste0(
        "DE_geneScores_", selected_cell_type_atac, "_", selected_region,
        "_", selected_cond_1, "_vs_", selected_cell_type_atac, "_",
        selected_region, "_", selected_cond_2, "_res.txt"
      )
    )
    res <- rna_sub %>%
      left_join(atac, by = c("gene" = "name")) %>%
      mutate(atac_sig = if_else(FDR.y < 0.05, "yes", "no")) %>%
      arrange(atac_sig)
    cor <- cor.test(res$model_log2FC, res$Log2FC, method = "spearman")
    tibble(
      cell_type_rna = selected_cell_type_rna,
      cell_type_atac = selected_cell_type_atac,
      cond_1 = selected_cond_1,
      cond_2 = selected_cond_2,
      region = selected_region,
      FDR_thres_snRNA_DE = FDR_threshold,
      FC_thres_snRNA_DE = FC_threshold,
      num_sig_snRNA_DEG = num_sig_DE,
      rho = cor$estimate,
      p = cor$p.value
    )
  }
}

# get_FC_cor("Astro", "Astro", "MCX", "ALS", "Control", 0.05, 2)

cell_types <- rna %>%
  distinct(cell_type) %>%
  filter(
    !cell_type %in% c("Exc_unknown", "Inh_ADARB2_Other", "Endo", "VLMC")
  ) %>%
  pull(cell_type)

params <- expand_grid(
  selected_cell_type_rna = cell_types,
  selected_cell_type_atac = cell_types,
  selected_region = c("MCX", "mFCX"),
  # selected_cond_1 = c("ALS", "FTD"),
  selected_cond_1 = "ALS",
  selected_cond_2 = "Control",
  FDR_threshold = 0.05,
  FC_threshold = c(1.2, 1.5, 2)
)

cor_res <- pmap_dfr(params, get_FC_cor)

write_tsv(cor_res, "spearman_cor_snRNA_expFC_vs_snATAC_geneScoreFC.txt")

# plot heatmap of rho
region_relabel <- c(
  "MCX" = "motor cortex",
  "mFCX" = "frontal cortex"
)

plot_cor_heatmap <- function(to_plot_cond_1,
                             to_plot_cond_2,
                             to_plot_region,
                             to_plot_FC_thres) {
  cor_res_sub <- cor_res %>%
    filter(
      cond_1 == to_plot_cond_1,
      cond_2 == to_plot_cond_2,
      region == to_plot_region,
      FDR_thres_snRNA_DE == 0.05,
      FC_thres_snRNA_DE == to_plot_FC_thres
    )

  mat_rho <- cor_res_sub %>%
    select(cell_type_rna, num_sig_snRNA_DEG, cell_type_atac, rho) %>%
    pivot_wider(
      names_from = cell_type_atac,
      values_from = rho,
      values_fill = NA_real_
    ) %>%
    mutate(cell_type_rna = factor(
      cell_type_rna,
      levels = c(
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST",
        "Astro", "Micro", "Oligo", "OPC"
      )
    )) %>%
    arrange(cell_type_rna) %>%
    mutate(
      row_label = paste0(cell_type_rna, " (", num_sig_snRNA_DEG, ")")
    ) %>%
    select(
      row_label,
      Exc_superficial, Exc_intermediate, Exc_deep,
      Inh_LAMP5, Inh_VIP, Inh_PVALB, Inh_SST,
      Astro, Micro, Oligo, OPC
    ) %>%
    as.data.frame() %>%
    column_to_rownames("row_label")

  mat_sig <- cor_res_sub %>%
    # mutate(
    #   FDR = p.adjust(p, method = "fdr"),
    #   sig = case_when(
    #     p < 0.05 & p >= 0.01 ~ "*",
    #     p < 0.01 & p >= 0.001 ~ "**",
    #     p < 0.001 ~ "***",
    #     TRUE ~ ""
    #   )
    # ) %>%
    mutate(
      sig = case_when(
        p < 0.05 & p >= 0.01 ~ "*",
        p < 0.01 & p >= 0.001 ~ "**",
        p < 0.001 ~ "***",
        TRUE ~ ""
      )
    ) %>%
    select(cell_type_rna, num_sig_snRNA_DEG, cell_type_atac, sig) %>%
    pivot_wider(
      names_from = cell_type_atac,
      values_from = sig,
      values_fill = NA_character_
    ) %>%
    mutate(cell_type_rna = factor(
      cell_type_rna,
      levels = c(
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST",
        "Astro", "Micro", "Oligo", "OPC"
      )
    )) %>%
    arrange(cell_type_rna) %>%
    mutate(
      row_label = paste0(cell_type_rna, " (", num_sig_snRNA_DEG, ")")
    ) %>%
    select(
      row_label,
      Exc_superficial, Exc_intermediate, Exc_deep,
      Inh_LAMP5, Inh_VIP, Inh_PVALB, Inh_SST,
      Astro, Micro, Oligo, OPC
    ) %>%
    as.data.frame() %>%
    column_to_rownames("row_label")

  pheatmap(
    mat_rho,
    color = scico(100, palette = "vik", direction = 1),
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
    display_numbers = mat_sig,
    number_color = "#000000",
    main = str_glue(
      "Spearman correlation of snRNA expression FC and snATAC gene score FC\n",
      "{to_plot_cond_1} vs. {to_plot_cond_2}, {region_relabel[to_plot_region]} \n",
      "in significnat snRNA DE genes (FDR < 0.05, FC > {to_plot_FC_thres})"
    ),
    filename = str_glue(
      "./plot/snRNA_expFC_vs_snATAC_geneScore_FC_spearman_cor_heatmap_",
      "{to_plot_cond_1}_vs_{to_plot_cond_2}_{to_plot_region}_FC_{to_plot_FC_thres}.pdf"
    ),
    width = 6,
    height = 4
  )
}

# plot_cor_heatmap("ALS", "Control", "MCX", 2)

params_plot <- expand_grid(
  to_plot_cond_1 = "ALS",
  to_plot_cond_2 = "Control",
  to_plot_region = c("MCX", "mFCX"),
  to_plot_FC_thres = c(1.2, 1.5, 2)
)

pwalk(params_plot, plot_cor_heatmap)

# log session info
sessionInfo()