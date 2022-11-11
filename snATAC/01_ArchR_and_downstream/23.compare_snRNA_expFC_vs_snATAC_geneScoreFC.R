# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)

rna <- read_tsv("MAST_res_level2_summary.txt")

plot_FC_cor <- function(selected_cell_type, selected_region, FDR_threshold, FC_threshold){
  message(paste0("Processing ", selected_cell_type, " in ", selected_region, " with FC: ", FC_threshold, " and FDR: ", FDR_threshold ))
  rna_sub <- rna %>%
    filter(cell_type == selected_cell_type,
           cond_1 == "ALS",
           cond_2 == "Control",
           region == selected_region,
           FDR < FDR_threshold,
           abs(model_log2FC) > log2(FC_threshold),
           conv_C == TRUE,
           conv_D == TRUE,
           model_log2FC_ci_hi * model_log2FC_ci_low > 0,
           abs(model_log2FC - avg_logFC) < 2
    )
  if(nrow(rna_sub) == 0){
    message("No significant DE genes in snRNA!")
  } else {
    atac <- read_tsv(paste0("DE_geneScores_", selected_cell_type, "_", selected_region,
                            "_ALS_vs_", selected_cell_type, "_", selected_region, "_Control_res.txt"))
    res <- rna_sub %>%
      left_join(atac, by = c("gene" = "name")) %>%
      mutate(atac_sig = if_else(FDR.y < 0.05, "yes", "no")) %>%
      arrange(atac_sig)
    cor <- cor.test(res$model_log2FC, res$Log2FC, method = "spearman")
    p <- res %>% ggplot(aes(model_log2FC, Log2FC))
    p1 <- p +
      geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
      geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
      geom_point(aes(color = atac_sig), size = 0.5) +
      scale_color_manual(values = c("yes" = "#CD2526", "no" = "#000000")) +
      xlab("snRNA gene expresssion log2FC\n(ALS vs. Control)") +
      ylab("snATAC gene score log2FC\n(ALS vs. Control)") +
      ggtitle(paste0(selected_cell_type, ", ", selected_region, "\n",
                     "FC > ", FC_threshold, ", FDR < ", FDR_threshold, "\n",
                     "rho = ", round(cor$estimate, 3), ", p = ", signif(cor$p.value, 3))) +
      theme_bw(base_size = 8, base_family = "Helvetica") +
      theme(panel.grid.minor = element_blank())
    ggsave(paste0("./plot/ALS_vs_Control_snRNA_expFC_vs_snATAC_genenScoreFC_in_", selected_cell_type, "_",
                  selected_region, "_FC_", FC_threshold, "_FDR_", FDR_threshold,".pdf"),
           plot = p1, device = cairo_pdf(), width = 4, height = 3, useDingbats = FALSE)
    dev.off()
  }
}

# plot_FC_cor("Astro", "MCX", 0.05, 2)

cell_types <- rna %>% distinct(cell_type) %>% filter(!cell_type %in% c("Exc_unknown", "Inh_ADARB2_Other", "Endo", "VLMC")) %>% pull(cell_type)
walk(cell_types, plot_FC_cor, selected_region = "MCX", FDR_threshold = 0.05, FC_threshold = 2)
walk(cell_types, plot_FC_cor, selected_region = "mFCX", FDR_threshold = 0.05, FC_threshold = 2)

walk(cell_types, plot_FC_cor, selected_region = "MCX", FDR_threshold = 0.05, FC_threshold = 1.5)
walk(cell_types, plot_FC_cor, selected_region = "mFCX", FDR_threshold = 0.05, FC_threshold = 1.5)

walk(cell_types, plot_FC_cor, selected_region = "MCX", FDR_threshold = 0.05, FC_threshold = 1.2)
walk(cell_types, plot_FC_cor, selected_region = "mFCX", FDR_threshold = 0.05, FC_threshold = 1.2)

sessionInfo()
