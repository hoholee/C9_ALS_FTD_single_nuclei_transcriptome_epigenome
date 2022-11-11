# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(scico)
library(ggrepel)

# read DE list
DE_level2 <- read_tsv("MAST_DE_res_level2_summary.txt")

DE_level2_Astro <- DE_level2 %>% 
  filter(cond_1 == "ALS",
         cond_2 == "Control",
         cell_type == "Astro",
         region == "MCX")

DE_level2_Micro <- DE_level2 %>% 
  filter(cond_1 == "ALS",
         cond_2 == "Control",
         cell_type == "Micro",
         region == "MCX")

DE_level2_Oligo <- DE_level2 %>% 
  filter(cond_1 == "ALS",
         cond_2 == "Control",
         cell_type == "Oligo",
         region == "MCX")

# read gene promoter H3K27ac Chip-seq signal
promoter_ChIP <- read_tsv("H3K27ac_signal_in_promoter_TSS_upstream2000_downstream1000.tsv") %>% 
  rename_all(function(x)(str_replace_all(x, "'", ""))) %>% 
  rename(chr = `#chr`) %>% 
  pivot_longer(cols = matches("^(ALS_)|(Control_)"), names_to = "sample", values_to = "value") %>% 
  separate(sample, c("disease", "subject", "cell_type"), sep = "_") %>% 
  mutate(value = if_else(is.na(value), 0, value))

gene_promoter_bed <- read_tsv("promoter_TSS_upstream2000_downstream1000.bed",
                              col_names = c("chr", "start", "end", "id", "score", "strand"),
                              col_types = "cddcic") %>% 
  separate(id, c("gene_id", "gene_name"), sep = "_") %>% 
  select(-score)

ChIP_res <- promoter_ChIP %>% 
  left_join(gene_promoter_bed) %>% 
  group_by(gene_id, gene_name, disease, cell_type) %>% 
  summarise(mean_value = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "disease", values_from = "mean_value")

# sanity check: make a figure of some cell-type specific genes
# box plot with individuals
# SNAP25, RBFOX3, APOE, GFAP, AQP4, C1QB, P2RY12, MOG, MOBP

gene_list <- c("SNAP25", "RBFOX3", "ETNPPL", "AQP4", "C1QB", "P2RY12", "MOG", "MOBP")

marker_genes <- promoter_ChIP %>% 
  left_join(gene_promoter_bed) %>% 
  filter(gene_name %in% gene_list) %>% 
  mutate(gene_name = factor(gene_name, levels = gene_list))

p <- marker_genes %>% ggplot(aes(cell_type, value))
p + geom_boxplot(aes(fill = disease), outlier.size = 0.5, size = 0.5) +
  facet_wrap(~gene_name, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("ALS" = "#EE436D", Control = "#1BB7A8")) +
  xlab("FANS-sorted cell type") +
  ylab("H3K27ac ChIP-seq signal (RPM) in promoter\n[TSS-2kb, TSS+1kb]") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("H3K27ac_ChIPseq_signal_in_cellType_marker_genes.pdf", device = cairo_pdf(), width = 6.5, height = 3.5, useDingbats = FALSE)

# scatter plot with mean, highlight some gene
control_res <- promoter_ChIP %>% 
  filter(disease == "Control") %>% 
  left_join(gene_promoter_bed) %>% 
  group_by(gene_id, gene_name, disease, cell_type) %>% 
  summarise(mean_value = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "cell_type", values_from = "mean_value")

# astro vs micro
gene_list1 <- c("ETNPPL", "AQP4", "C1QB", "P2RY12")
p <- control_res %>% ggplot(aes(Astro, MG))
p + 
  geom_hex(bins = 100) +
  geom_point(data = . %>% filter(gene_name %in% gene_list1), shape = 1, color = "#C91878") +
  geom_text_repel(data = . %>% filter(gene_name %in% gene_list1), color = "#C91878", aes(label = gene_name), size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey", linetype = 2) +
  scale_fill_viridis_c(trans = "log10", name = "# genes", option = "D") +
  coord_cartesian(xlim = c(0, 2.1), ylim = c(0, 2.1)) +
  xlab("Mean H3K27ac ChIP-seq signal (RPM) in Astrocytes of Control samples\nin promoter [TSS-2kb, TSS+1kb]") +
  ylab("Mean H3K27ac ChIP-seq signal (RPM)\nin Microglia of Control samples\nin promoter [TSS-2kb, TSS+1kb]") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("H3K27ac_ChIPseq_signal_in_Control_Astro_vs_Micro.pdf", device = cairo_pdf(), width = 5, height = 3, useDingbats = FALSE)  


# neuron vs oligo
gene_list2 <- c("SNAP25", "RBFOX3", "MOG", "MOBP")
p <- control_res %>% ggplot(aes(NeuN, Olig))
p + 
  geom_hex(bins = 100) +
  geom_point(data = . %>% filter(gene_name %in% gene_list2), shape = 1, color = "#C91878") +
  geom_text_repel(data = . %>% filter(gene_name %in% gene_list2), color = "#C91878", aes(label = gene_name), size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey", linetype = 2) +
  scale_fill_viridis_c(trans = "log10", name = "# genes", option = "D") +
  coord_cartesian(xlim = c(0, 1.8), ylim = c(0, 1.8)) +
  xlab("Mean H3K27ac ChIP-seq signal (RPM) in Neurons of Control samples\nin promoter [TSS-2kb, TSS+1kb]") +
  ylab("Mean H3K27ac ChIP-seq signal (RPM)\nin Oligodendrocytes of Control samples\nin promoter [TSS-2kb, TSS+1kb]") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("H3K27ac_ChIPseq_signal_in_Control_Neuron_vs_Oligo.pdf", device = cairo_pdf(), width = 5, height = 3, useDingbats = FALSE)  

# disease effect
p <- ChIP_res %>% ggplot(aes(Control, ALS))
p + 
  geom_hex(bins = 100) +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey", linetype = 2) +
  facet_wrap(~cell_type, scales = "free") +
  scale_fill_viridis_c(trans = "log10", name = "# genes", option = "D") +
  xlab("Mean H3K27ac ChIP-seq signal (RPM) in Control samples\nin promoter [TSS-2kb, TSS+1kb]") +
  ylab("Mean H3K27ac ChIP-seq signal (RPM)\nin ALS samples\nin promoter [TSS-2kb, TSS+1kb]") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("H3K27ac_ChIPseq_signal_Control_vs_ALS.pdf", device = cairo_pdf(), width = 6, height = 5, useDingbats = FALSE)  


# correlates gene expression FC in snRNA and delta H3K27ac signal in ChIP-seq
# Astro
astro_diff <- ChIP_res %>% 
  filter(cell_type == "Astro") %>% 
  mutate(delta = ALS - Control) %>% 
  left_join(DE_level2_Astro, by = c("gene_name" = "gene"))

astro_diff_sig <- astro_diff %>%
  filter(FDR < 0.05,
         abs(model_log2FC) > log2(1.5), 
         model_log2FC_ci_hi *  model_log2FC_ci_low > 0,
         conv_C == TRUE,
         conv_D == TRUE,
         abs(model_log2FC - avg_log2FC) < 2)

p <- astro_diff_sig %>% 
  ggplot(aes(model_log2FC, delta))

p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point() +
  geom_point(data = . %>% filter(abs(delta) > 0.05, model_log2FC * delta > 0), color = "#C91878") + 
  # geom_text_repel(data = . %>% filter(abs(model_log2FC) > 1.8, abs(delta) > 0.05, model_log2FC * delta > 0),
  geom_text_repel(data = . %>% filter(abs(delta) > 0.05, model_log2FC * delta > 0),
                  aes(label = gene_name), color = "#C91878") + 
  xlab("snRNA gene expression log2FC") +
  ylab("delta H3K27ac signal at promoter (RPM)") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("snRNA_FC_vs_delta_H3K27ac_at_promoter_Astro_DE_genes.pdf", device = cairo_pdf(), width = 6, height = 5, useDingbats = FALSE)  


gene_list <- c("AQP4", "ITGB4", "ARHGEF10L", "CD44", "AQP1", "LINC00499", "C9orf72", "ARHGAP31", "GRM3", "GRID1")

marker_genes <- promoter_ChIP %>% 
  filter(cell_type == "Astro") %>% 
  left_join(gene_promoter_bed) %>% 
  filter(gene_name %in% gene_list) %>% 
  mutate(gene_name = factor(gene_name, levels = gene_list))

p <- marker_genes %>% ggplot(aes(gene_name, value))
p + geom_boxplot(aes(fill = disease), outlier.size = 0.5, size = 0.5, width = 0.6) +
  scale_fill_manual(values = c("ALS" = "#EE436D", Control = "#1BB7A8")) +
  xlab("Gene") +
  ylab("H3K27ac ChIP-seq signal (RPM) in promoter\n[TSS-2kb, TSS+1kb]") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("boxplot_H3K27ac_at_promoter_selected_Astro_DE_genes.pdf", device = cairo_pdf(), width = 6, height = 5, useDingbats = FALSE)  

# Micro
micro_diff <- ChIP_res %>% 
  filter(cell_type == "MG") %>% 
  mutate(delta = ALS - Control) %>% 
  left_join(DE_level2_Micro, by = c("gene_name" = "gene"))

micro_diff_sig <- micro_diff %>%
  filter(FDR < 0.05,
         abs(model_log2FC) > log2(1.5), 
         model_log2FC_ci_hi *  model_log2FC_ci_low > 0,
         conv_C == TRUE,
         conv_D == TRUE,
         abs(model_log2FC - avg_log2FC) < 2)

p <- micro_diff_sig %>% 
  ggplot(aes(model_log2FC, delta))

p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point() +
  geom_point(data = . %>% filter(abs(delta) > 0.05, model_log2FC * delta > 0), color = "#C91878") + 
  # geom_text_repel(data = . %>% filter(abs(model_log2FC) > 1.8, abs(delta) > 0.05, model_log2FC * delta > 0),
  geom_text_repel(data = . %>% filter(abs(delta) > 0.05, model_log2FC * delta > 0),
                  aes(label = gene_name), color = "#C91878") + 
  xlab("snRNA gene expression log2FC") +
  ylab("delta H3K27ac signal at promoter (RPM)") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("snRNA_FC_vs_delta_H3K27ac_at_promoter_Micro_DE_genes.pdf", device = cairo_pdf(), width = 6, height = 5, useDingbats = FALSE)  


gene_list <- c("DSCAM", "MITF", "MYO1E", "PPARG", "SPP1", "MRC1", "RNF150", "LRRFIP1", "CD300A", "LINC01684")

marker_genes <- promoter_ChIP %>% 
  filter(cell_type == "MG") %>% 
  left_join(gene_promoter_bed) %>% 
  filter(gene_name %in% gene_list) %>% 
  mutate(gene_name = factor(gene_name, levels = gene_list))

p <- marker_genes %>% ggplot(aes(gene_name, value))
p + geom_boxplot(aes(fill = disease), outlier.size = 0.5, size = 0.5, width = 0.6) +
  scale_fill_manual(values = c("ALS" = "#EE436D", Control = "#1BB7A8")) +
  xlab("Gene") +
  ylab("H3K27ac ChIP-seq signal (RPM) in promoter\n[TSS-2kb, TSS+1kb]") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("boxplot_H3K27ac_at_promoter_selected_Micro_DE_genes.pdf", device = cairo_pdf(), width = 6, height = 5, useDingbats = FALSE)  

# Astro
oligo_diff <- ChIP_res %>% 
  filter(cell_type == "Olig") %>% 
  mutate(delta = ALS - Control) %>% 
  left_join(DE_level2_Astro, by = c("gene_name" = "gene"))

oligo_diff_sig <- oligo_diff %>%
  filter(FDR < 0.05,
         abs(model_log2FC) > log2(1.5), 
         model_log2FC_ci_hi *  model_log2FC_ci_low > 0,
         conv_C == TRUE,
         conv_D == TRUE,
         abs(model_log2FC - avg_log2FC) < 2)

p <- oligo_diff_sig %>% 
  ggplot(aes(model_log2FC, delta))

p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point() +
  geom_point(data = . %>% filter(abs(delta) > 0.05, model_log2FC * delta > 0), color = "#C91878") + 
  # geom_text_repel(data = . %>% filter(abs(model_log2FC) > 1.8, abs(delta) > 0.05, model_log2FC * delta > 0),
  geom_text_repel(data = . %>% filter(abs(delta) > 0.05, model_log2FC * delta > 0),
                  aes(label = gene_name), color = "#C91878") + 
  xlab("snRNA gene expression log2FC") +
  ylab("delta H3K27ac signal at promoter (RPM)") +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("snRNA_FC_vs_delta_H3K27ac_at_promoter_Oligo_DE_genes.pdf", device = cairo_pdf(), width = 6, height = 5, useDingbats = FALSE)  


gene_list <- c("BCL6", "AR1D5B", "KLF6", "NEDD4", "CLK1", "HSP90AB1", "C9orf72")

marker_genes <- promoter_ChIP %>% 
  filter(cell_type == "Olig") %>% 
  left_join(gene_promoter_bed) %>% 
  filter(gene_name %in% gene_list) %>% 
  mutate(gene_name = factor(gene_name, levels = gene_list))

p <- marker_genes %>% ggplot(aes(gene_name, value))
p + geom_boxplot(aes(fill = disease), outlier.size = 0.5, size = 0.5, width = 0.6) +
  scale_fill_manual(values = c("ALS" = "#EE436D", Control = "#1BB7A8")) +
  xlab("Gene") +
  ylab("H3K27ac ChIP-seq signal (RPM) in promoter\n[TSS-2kb, TSS+1kb]") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave("boxplot_H3K27ac_at_promoter_selected_Oligo_DE_genes.pdf", device = cairo_pdf(), width = 6, height = 5, useDingbats = FALSE)  
