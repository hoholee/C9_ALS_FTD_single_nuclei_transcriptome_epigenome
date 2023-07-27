# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(scico)
library(ggrepel)
library(plotly)
library(htmlwidgets)

# read DE list
DE_level2 <- read_tsv("/cndd2/junhao/ALS_FTD_singleCell/snRNA_postCellBender/DE_MAST/second_round_level_2/MAST_res_level2_summary.txt")

DE_level2_Astro <- DE_level2 %>%
  filter(
    cond_1 == "ALS",
    cond_2 == "Control",
    cell_type == "Astro",
    region == "MCX"
  )

DE_level2_Micro <- DE_level2 %>%
  filter(
    cond_1 == "ALS",
    cond_2 == "Control",
    cell_type == "Micro",
    region == "MCX"
  )

DE_level2_Oligo <- DE_level2 %>%
  filter(
    cond_1 == "ALS",
    cond_2 == "Control",
    cell_type == "Oligo",
    region == "MCX"
  )

DE_level2_Exc_upper <- DE_level2 %>%
  filter(
    cond_1 == "ALS",
    cond_2 == "Control",
    cell_type == "Exc_superficial",
    region == "MCX"
  )

DE_level2_Exc_intermediate <- DE_level2 %>%
  filter(
    cond_1 == "ALS",
    cond_2 == "Control",
    cell_type == "Exc_intermediate",
    region == "MCX"
  )

DE_level2_Exc_deep <- DE_level2 %>%
  filter(
    cond_1 == "ALS",
    cond_2 == "Control",
    cell_type == "Exc_deep",
    region == "MCX"
  )


# read gene promoter H3K27ac Chip-seq signal
promoter_ChIP <- read_tsv("H3K27ac_signal_in_promoter_TSS_upstream2000_downstream1000.tsv") %>%
  rename_all(function(x) (str_replace_all(x, "'", ""))) %>%
  rename(chr = `#chr`) %>%
  pivot_longer(cols = matches("^(ALS_)|(Control_)"), names_to = "sample", values_to = "value") %>%
  separate(sample, c("disease", "subject", "cell_type"), sep = "_") %>%
  mutate(value = if_else(is.na(value), 0, value))

gene_promoter_bed <- read_tsv("promoter_TSS_upstream2000_downstream1000.bed",
  col_names = c("chr", "start", "end", "id", "score", "strand"),
  col_types = "cddcic"
) %>%
  separate(id, c("gene_id", "gene_name"), sep = "_") %>%
  select(-score)

ChIP_res <- promoter_ChIP %>%
  left_join(gene_promoter_bed) %>%
  group_by(gene_id, gene_name, disease, cell_type) %>%
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = "disease", values_from = "mean_value")

# correlates gene expression FC in snRNA and delta H3K27ac signal in ChIP-seq
pseudo_count <- 0.001

# Astro
astro_diff <- ChIP_res %>%
  filter(cell_type == "Astro") %>%
  mutate(delta = ALS - Control) %>%
  mutate(log2FC = log2(ALS + pseudo_count) - log2(Control + pseudo_count)) %>%
  left_join(DE_level2_Astro, by = c("gene_name" = "gene")) %>%
  rename(
    FANS_cell_type = cell_type.x,
    promoter_H3K27ac_signal_ALS = ALS,
    promoter_H3K27ac_signal_Control = Control,
    promoter_H3K27ac_diff = delta,
    promoter_H3K27ac_log2FC = log2FC,
    snRNA_cell_type = cell_type.y,
  )

astro_diff_sig <- astro_diff %>%
  filter(
    FDR < 0.05,
    abs(model_log2FC) > log2(1.2),
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    conv_C == TRUE,
    conv_D == TRUE,
    abs(model_log2FC - avg_logFC) < 2
  )

write_tsv(astro_diff_sig, "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Astro_DE_genes.tsv")

p <- astro_diff_sig %>%
  ggplot(aes(model_log2FC, promoter_H3K27ac_log2FC))

p <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(aes(label = gene_name)) +
  xlab("snRNA gene expression log2FC") +
  ylab("H3K27ac signal at promoter log2FC") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )
p1 <- ggplotly(p)
saveWidget(
  p1,
  "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Astro_DE_genes.html",
  selfcontained = TRUE
)


# Micro
micro_diff <- ChIP_res %>%
  filter(cell_type == "MG") %>%
  mutate(delta = ALS - Control) %>%
  mutate(log2FC = log2(ALS + pseudo_count) - log2(Control + pseudo_count)) %>%
  left_join(DE_level2_Micro, by = c("gene_name" = "gene")) %>%
  rename(
    FANS_cell_type = cell_type.x,
    promoter_H3K27ac_signal_ALS = ALS,
    promoter_H3K27ac_signal_Control = Control,
    promoter_H3K27ac_diff = delta,
    promoter_H3K27ac_log2FC = log2FC,
    snRNA_cell_type = cell_type.y,
  )

micro_diff_sig <- micro_diff %>%
  filter(
    FDR < 0.05,
    abs(model_log2FC) > log2(1.2),
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    conv_C == TRUE,
    conv_D == TRUE,
    abs(model_log2FC - avg_logFC) < 2
  )
write_tsv(micro_diff_sig, "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Micro_DE_genes.tsv")

p <- micro_diff_sig %>%
  ggplot(aes(model_log2FC, promoter_H3K27ac_log2FC))

p <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(aes(label = gene_name)) +
  xlab("snRNA gene expression log2FC") +
  ylab("H3K27ac signal at promoter log2FC") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )
p1 <- ggplotly(p)
saveWidget(
  p1,
  "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Micro_DE_genes.html",
  selfcontained = TRUE
)

# Oligo
oligo_diff <- ChIP_res %>%
  filter(cell_type == "Olig") %>%
  mutate(delta = ALS - Control) %>%
  mutate(log2FC = log2(ALS + pseudo_count) - log2(Control + pseudo_count)) %>%
  left_join(DE_level2_Oligo, by = c("gene_name" = "gene")) %>%
  rename(
    FANS_cell_type = cell_type.x,
    promoter_H3K27ac_signal_ALS = ALS,
    promoter_H3K27ac_signal_Control = Control,
    promoter_H3K27ac_diff = delta,
    promoter_H3K27ac_log2FC = log2FC,
    snRNA_cell_type = cell_type.y,
  )

oligo_diff_sig <- oligo_diff %>%
  filter(
    FDR < 0.05,
    abs(model_log2FC) > log2(1.2),
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    conv_C == TRUE,
    conv_D == TRUE,
    abs(model_log2FC - avg_logFC) < 2
  )

write_tsv(oligo_diff_sig, "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Oligo_DE_genes.tsv")

p <- oligo_diff_sig %>%
  ggplot(aes(model_log2FC, promoter_H3K27ac_log2FC))

p <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(aes(label = gene_name)) +
  xlab("snRNA gene expression log2FC") +
  ylab("H3K27ac signal at promoter log2FC") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )
p1 <- ggplotly(p)
saveWidget(
  p1,
  "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Oligo_DE_genes.html",
  selfcontained = TRUE
)



# Exc upper
exc_upper_diff <- ChIP_res %>%
  filter(cell_type == "NeuN") %>%
  mutate(delta = ALS - Control) %>%
  mutate(log2FC = log2(ALS + pseudo_count) - log2(Control + pseudo_count)) %>%
  left_join(DE_level2_Exc_upper, by = c("gene_name" = "gene")) %>%
  rename(
    FANS_cell_type = cell_type.x,
    promoter_H3K27ac_signal_ALS = ALS,
    promoter_H3K27ac_signal_Control = Control,
    promoter_H3K27ac_diff = delta,
    promoter_H3K27ac_log2FC = log2FC,
    snRNA_cell_type = cell_type.y,
  )

exc_upper_diff_sig <- exc_upper_diff %>%
  filter(
    FDR < 0.05,
    abs(model_log2FC) > log2(1.2),
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    conv_C == TRUE,
    conv_D == TRUE,
    abs(model_log2FC - avg_logFC) < 2
  )

write_tsv(exc_upper_diff_sig, "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Exc_upper_DE_genes.tsv")

p <- exc_upper_diff_sig %>%
  ggplot(aes(model_log2FC, promoter_H3K27ac_log2FC))

p <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(aes(label = gene_name)) +
  xlab("snRNA gene expression log2FC") +
  ylab("H3K27ac signal at promoter log2FC") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )
p1 <- ggplotly(p)
saveWidget(
  p1,
  "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Exc_upper_DE_genes.html",
  selfcontained = TRUE
)


# Exc intermediate
exc_intermediate_diff <- ChIP_res %>%
  filter(cell_type == "NeuN") %>%
  mutate(delta = ALS - Control) %>%
  mutate(log2FC = log2(ALS + pseudo_count) - log2(Control + pseudo_count)) %>%
  left_join(DE_level2_Exc_intermediate, by = c("gene_name" = "gene")) %>%
  rename(
    FANS_cell_type = cell_type.x,
    promoter_H3K27ac_signal_ALS = ALS,
    promoter_H3K27ac_signal_Control = Control,
    promoter_H3K27ac_diff = delta,
    promoter_H3K27ac_log2FC = log2FC,
    snRNA_cell_type = cell_type.y,
  )

exc_intermediate_diff_sig <- exc_intermediate_diff %>%
  filter(
    FDR < 0.05,
    abs(model_log2FC) > log2(1.2),
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    conv_C == TRUE,
    conv_D == TRUE,
    abs(model_log2FC - avg_logFC) < 2
  )

write_tsv(exc_intermediate_diff_sig, "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Exc_intermediate_DE_genes.tsv")

p <- exc_intermediate_diff_sig %>%
  ggplot(aes(model_log2FC, promoter_H3K27ac_log2FC))

p <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(aes(label = gene_name)) +
  xlab("snRNA gene expression log2FC") +
  ylab("H3K27ac signal at promoter log2FC") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )
p1 <- ggplotly(p)
saveWidget(
  p1,
  "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Exc_intermediate_DE_genes.html",
  selfcontained = TRUE
)

# Exc deep
exc_deep_diff <- ChIP_res %>%
  filter(cell_type == "NeuN") %>%
  mutate(delta = ALS - Control) %>%
  mutate(log2FC = log2(ALS + pseudo_count) - log2(Control + pseudo_count)) %>%
  left_join(DE_level2_Exc_deep, by = c("gene_name" = "gene")) %>%
  rename(
    FANS_cell_type = cell_type.x,
    promoter_H3K27ac_signal_ALS = ALS,
    promoter_H3K27ac_signal_Control = Control,
    promoter_H3K27ac_diff = delta,
    promoter_H3K27ac_log2FC = log2FC,
    snRNA_cell_type = cell_type.y,
  )

exc_deep_diff_sig <- exc_deep_diff %>%
  filter(
    FDR < 0.05,
    abs(model_log2FC) > log2(1.2),
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    conv_C == TRUE,
    conv_D == TRUE,
    abs(model_log2FC - avg_logFC) < 2
  )

write_tsv(exc_deep_diff_sig, "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Exc_deep_DE_genes.tsv")

p <- exc_deep_diff_sig %>%
  ggplot(aes(model_log2FC, promoter_H3K27ac_log2FC))

p <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(aes(label = gene_name)) +
  xlab("snRNA gene expression log2FC") +
  ylab("H3K27ac signal at promoter log2FC") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )
p1 <- ggplotly(p)
saveWidget(
  p1,
  "./snRNA_geneFC_vs_promoter_H3K27ac_signal/snRNA_FC_vs_promoter_H3K27ac_FC_Exc_deep_DE_genes.html",
  selfcontained = TRUE
)
