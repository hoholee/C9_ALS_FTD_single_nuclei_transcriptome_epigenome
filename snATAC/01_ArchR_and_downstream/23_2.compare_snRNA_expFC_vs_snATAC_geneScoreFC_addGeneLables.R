# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(plotly)
library(htmlwidgets)

rna <- read_tsv("MAST_res_level2_summary.txt")

cell_types <- rna %>%
  distinct(cell_type) %>%
  filter(!cell_type %in% c("Exc_unknown", "Inh_ADARB2_Other", "Endo", "VLMC")) %>%
  pull(cell_type)

save_res_table_and_html <- function(selected_cell_type, selected_region, FDR_threshold, FC_threshold) {
  message(paste0("Processing ", selected_cell_type, " in ", selected_region, " with FC: ", FC_threshold, " and FDR: ", FDR_threshold))

  rna_sub <- rna %>%
    filter(
      cell_type == selected_cell_type,
      cond_1 == "ALS",
      cond_2 == "Control",
      region == selected_region,
      FDR < FDR_threshold,
      abs(model_log2FC) > log2(FC_threshold),
      conv_C == TRUE,
      conv_D == TRUE,
      model_log2FC_ci_hi * model_log2FC_ci_low > 0,
      abs(model_log2FC - avg_logFC) < 2
    ) %>%
    rename(
      rna_num_cells_cond_1 = num_cells_cond_1,
      rna_num_cells_cond_2 = num_cells_cond_2,
      rna_p_value = p_value,
      rna_model_log2FC = model_log2FC,
      rna_model_log2FC_ci_hi = model_log2FC_ci_hi,
      rna_model_log2FC_ci_low = model_log2FC_ci_low,
      rna_FDR = FDR,
      rna_conv_C = conv_C,
      rna_conv_D = conv_D,
      rna_avg_log2CPM_cond_1 = avg_log2CPM_cond_1,
      rna_avg_log2CPM_cond_2 = avg_log2CPM_cond_2,
      rna_avg_log2FC = avg_logFC,
    )
  if (nrow(rna_sub) == 0) {
    message("No significant DE genes in snRNA!")
  } else {
    atac <- read_tsv(paste0(
      "DE_geneScores_", selected_cell_type, "_", selected_region,
      "_ALS_vs_", selected_cell_type, "_", selected_region, "_Control_res.txt"
    )) %>%
      select(-idx, -group_1, -group_2) %>%
      rename(
        chrom = seqnames,
        atac_log2FC = Log2FC,
        atac_FDR = FDR,
        atac_meanDiff = MeanDiff
      )
    res <- rna_sub %>%
      left_join(atac, by = c("gene" = "name")) %>%
      mutate(atac_sig = if_else(atac_FDR < 0.05, "yes", "no"))
    output_table <- str_glue(
      "./formatted_table_and_html/",
      "ALS_vs_Control_snRNA_expFC_vs_snATAC_genenScoreFC_in_",
      "{selected_cell_type}_{selected_region}_",
      "FC_{FC_threshold}_FDR_{FDR_threshold}.tsv"
    )
    write_tsv(res, output_table)

    cor <- cor.test(res$rna_model_log2FC, res$atac_log2FC, method = "spearman")
    p <- res %>% ggplot(aes(rna_model_log2FC, atac_log2FC))
    p <- p +
      geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
      geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
      geom_point(aes(color = atac_sig, label = gene), size = 0.5) +
      scale_color_manual(values = c("yes" = "#CD2526", "no" = "#000000")) +
      xlab("snRNA gene expresssion log2FC\n(ALS vs. Control)") +
      ylab("snATAC gene score log2FC\n(ALS vs. Control)") +
      ggtitle(paste0(
        selected_cell_type, ", ", selected_region, "\n",
        "FC > ", FC_threshold, ", FDR < ", FDR_threshold, "\n",
        "rho = ", round(cor$estimate, 3), ", p = ", signif(cor$p.value, 3)
      )) +
      theme_bw(base_size = 8, base_family = "Helvetica") +
      theme(panel.grid.minor = element_blank())

    p1 <- ggplotly(p)
    output_html <- str_glue(
      "./formatted_table_and_html/",
      "ALS_vs_Control_snRNA_expFC_vs_snATAC_genenScoreFC_in_",
      "{selected_cell_type}_{selected_region}_",
      "FC_{FC_threshold}_FDR_{FDR_threshold}.html"
    )
    saveWidget(p1, file = output_html, selfcontained = TRUE)
  }
}

walk(cell_types, save_res_table_and_html, selected_region = "MCX", FDR_threshold = 0.05, FC_threshold = 1.2)
walk(cell_types, save_res_table_and_html, selected_region = "mFCX", FDR_threshold = 0.05, FC_threshold = 1.2)

sessionInfo()
