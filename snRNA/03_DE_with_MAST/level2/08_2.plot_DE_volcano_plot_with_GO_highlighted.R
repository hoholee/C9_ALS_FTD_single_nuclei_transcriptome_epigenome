library(tidyverse)

DE <- read_tsv("../MAST_res_level2_summary.txt")

FDR_threshold <- 0.05
FC_threshold <- 1.2

genes_in_GO <- read_tsv("Astr.actin.overlap_genes.txt", col_names = "gene")
# genes_in_GO <- read_tsv("Astr.potassium_ion_transport.overlap_genes.txt", col_names = "gene")

DE_sub <- DE %>%
  filter(cond_1 == "ALS",
         cond_2 == "Control",
         region == "MCX",
         cell_type != "Exc_unknown") %>%
  mutate(sig = if_else(
    (FDR < FDR_threshold &
       abs(model_log2FC) > log2(FC_threshold) &
       conv_C == TRUE &
       conv_D == TRUE &
       model_log2FC_ci_hi * model_log2FC_ci_low > 0 &
       abs(model_log2FC - avg_logFC) < 2),
    "DE",
    "nonDE"
  ),
  in_GO = if_else(gene %in% genes_in_GO$gene, "in_GO", "not_in_GO")
  ) %>%
  arrange(sig)

p <- DE_sub %>% ggplot(aes(model_log2FC, -log10(FDR)))
p +
  geom_point(aes(color = sig, size = sig)) +
  geom_point(aes(size = sig), data = DE_sub %>% filter(in_GO == "in_GO"), color = "#3F619F") +
  facet_wrap(~cell_type) +
  scale_color_manual(values = c("#9C302F", "#B3AAAA")) +
  scale_size_manual(values = c(1, 0.2)) +
  coord_cartesian(xlim = c(-5, 5)) +
  xlab("log2FC") +
  ggtitle("MCX, ALS vs. Control, GO:Actin") +
  # ggtitle("MCX, ALS vs. Control, GO:K+ transport") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("volcano_plots_ALS_vs_Control_MCX_show_GO_actin.pdf",
# ggsave("volcano_plots_ALS_vs_Control_MCX_show_GO_Ktransport.pdf",
       device = cairo_pdf(),
       width = 6, height = 4, useDingbats = FALSE)
