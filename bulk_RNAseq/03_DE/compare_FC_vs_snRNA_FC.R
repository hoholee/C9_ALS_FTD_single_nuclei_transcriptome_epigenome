library(tidyverse)
library(viridis)
library(scico)

bulk_DE <- read_tsv("edgeR_res_glmTest_allExpressedGenes_rmOutlierSamples_simplified.txt")

# MA plot
p <- bulk_DE %>% ggplot(aes(logCPM, logFC))
p + 
  geom_hex(bins = 100) + 
  theme_bw(base_size = 10, base_family = "Helvetica") + 
  facet_wrap(~comparison + region + cell_type, nrow = 3) +
  scale_fill_viridis_c(option = "B") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/DE_MA_plots.pdf", device = cairo_pdf(), width = 8, height = 6, useDingbats = FALSE)

# read snRNA DE

snRNA_DE_level2 <- read_tsv("/cndd2/junhao/ALS_FTD_singleCell/snRNA_clean/DE_MAST/level_2/MAST_DE_res_level2_summary_expr_frac_0.1.txt")

snRNA_DE_level2_sig <- snRNA_DE_level2 %>% 
  filter(conv_C == TRUE,
         conv_D == TRUE,
         FDR < 0.05,
         abs(model_log2FC) > log2(1.1),
         model_log2FC_ci_hi * model_log2FC_ci_low > 0,
         abs(model_log2FC - avg_lo2FC) < 2)

# Astro vs. OtherGlia
snRNA_DE_level2_sig_astro <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Astro") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_glia <- bulk_DE %>%
  filter(cell_type == "Other_glias") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

astro_merged <- snRNA_DE_level2_sig_astro %>%
  left_join(bulk_DE_glia, by = c("comparison", "region", "gene" = "geneName"))

p <- astro_merged %>% 
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Astrocyte, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Astrocyte + Microglia, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_astro_FC_vs_FANSsorted_otherGlia_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

# Micro vs. OtherGlia
snRNA_DE_level2_sig_micro <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Micro") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_glia <- bulk_DE %>%
  filter(cell_type == "Other_glias") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

micro_merged <- snRNA_DE_level2_sig_micro %>%
  left_join(bulk_DE_glia, by = c("comparison", "region", "gene" = "geneName"))

p <- micro_merged %>% 
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Microglia, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Astrocyte + Microglia, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_micro_FC_vs_FANSsorted_otherGlia_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)


# Oligo vs. Oligo
snRNA_DE_level2_sig_oligo <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Oligo") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_oligo <- bulk_DE %>%
  filter(cell_type == "Oligodendrocytes") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

oligo_merged <- snRNA_DE_level2_sig_oligo %>%
  left_join(bulk_DE_oligo, by = c("comparison", "region", "gene" = "geneName"))

p <- oligo_merged %>% 
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Oligodendrocyte, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Oligodendrocyte, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_oligo_FC_vs_FANSsorted_oligo_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)


# OPC vs. Oligo
snRNA_DE_level2_sig_OPC <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "OPC") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_oligo <- bulk_DE %>%
  filter(cell_type == "Oligodendrocytes") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

OPC_merged <- snRNA_DE_level2_sig_OPC %>%
  left_join(bulk_DE_oligo, by = c("comparison", "region", "gene" = "geneName"))

p <- OPC_merged %>% 
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, OPC, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Oligodendrocyte, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_OPC_FC_vs_FANSsorted_oligo_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

# Exc L2 vs Neurons
snRNA_DE_level2_sig_ExcL2 <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Exc_L2") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_neurons <- bulk_DE %>%
  filter(cell_type == "Neurons") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

ExcL2_merged <- snRNA_DE_level2_sig_ExcL2 %>%
  left_join(bulk_DE_neurons, by = c("comparison", "region", "gene" = "geneName"))

p <- ExcL2_merged %>%
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Exc L2, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Neurons, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_ExcL2_FC_vs_FANSsorted_neuron_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

# Exc L3-5 vs Neurons
snRNA_DE_level2_sig_ExcL35 <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Exc_L3-5") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_neurons <- bulk_DE %>%
  filter(cell_type == "Neurons") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

ExcL35_merged <- snRNA_DE_level2_sig_ExcL35 %>%
  left_join(bulk_DE_neurons, by = c("comparison", "region", "gene" = "geneName"))

p <- ExcL35_merged %>%
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Exc L3-5, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Neurons, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_ExcL35_FC_vs_FANSsorted_neuron_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)


# Exc L5-6 vs Neurons
snRNA_DE_level2_sig_ExcL56 <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Exc_L5-6") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_neurons <- bulk_DE %>%
  filter(cell_type == "Neurons") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

ExcL56_merged <- snRNA_DE_level2_sig_ExcL56 %>%
  left_join(bulk_DE_neurons, by = c("comparison", "region", "gene" = "geneName"))

p <- ExcL56_merged %>%
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Exc L5-6, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Neurons, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_ExcL56_FC_vs_FANSsorted_neuron_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

## Inh CGE LAMP5 vs. Neurons
snRNA_DE_level2_sig_InhCGELamp5 <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Inh_CGE_LAMP5") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_neurons <- bulk_DE %>%
  filter(cell_type == "Neurons") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

InhCGELamp5_merged <- snRNA_DE_level2_sig_InhCGELamp5 %>%
  left_join(bulk_DE_neurons, by = c("comparison", "region", "gene" = "geneName"))

p <- InhCGELamp5_merged %>%
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Inh CGE LAMP5, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Neurons, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_InhCGELamp5_FC_vs_FANSsorted_neuron_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

## Inh CGE VIP vs. Neurons
snRNA_DE_level2_sig_InhCGEVip <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Inh_CGE_VIP") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_neurons <- bulk_DE %>%
  filter(cell_type == "Neurons") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

InhCGEVip_merged <- snRNA_DE_level2_sig_InhCGEVip %>%
  left_join(bulk_DE_neurons, by = c("comparison", "region", "gene" = "geneName"))

p <- InhCGEVip_merged %>%
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Inh CGE VIP, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Neurons, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_InhCGEVip_FC_vs_FANSsorted_neuron_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

## Inh CGE Other vs. Neurons
snRNA_DE_level2_sig_InhCGEOther <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Inh_CGE_Other") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_neurons <- bulk_DE %>%
  filter(cell_type == "Neurons") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

InhCGEOther_merged <- snRNA_DE_level2_sig_InhCGEOther %>%
  left_join(bulk_DE_neurons, by = c("comparison", "region", "gene" = "geneName"))

p <- InhCGEOther_merged %>%
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Inh CGE Other, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Neurons, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_InhCGEOther_FC_vs_FANSsorted_neuron_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

## Inh MGE PVALB vs. Neurons
snRNA_DE_level2_sig_InhMGEPvalb <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Inh_MGE_PVALB") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_neurons <- bulk_DE %>%
  filter(cell_type == "Neurons") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

InhMGEPvalb_merged <- snRNA_DE_level2_sig_InhMGEPvalb %>%
  left_join(bulk_DE_neurons, by = c("comparison", "region", "gene" = "geneName"))

p <- InhMGEPvalb_merged %>%
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Inh MGE PVALB, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Neurons, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_InhMGEPvalb_FC_vs_FANSsorted_neuron_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

## Inh MGE SST vs. Neurons
snRNA_DE_level2_sig_InhMGESst <- snRNA_DE_level2_sig %>% 
  filter(cell_type == "Inh_MGE_SST") %>% 
  mutate(comparison = paste0(cond_1, "_vs_", cond_2))
bulk_DE_neurons <- bulk_DE %>%
  filter(cell_type == "Neurons") %>% 
  mutate(region = if_else(region == "motor_cortex", "MCX", "mFCX"))

InhMGESst_merged <- snRNA_DE_level2_sig_InhMGESst %>%
  left_join(bulk_DE_neurons, by = c("comparison", "region", "gene" = "geneName"))

p <- InhMGESst_merged %>%
  ggplot(aes(model_log2FC, logFC))

p + geom_hex(bins = 100) +
  facet_wrap(~comparison + region, dir = "v", nrow = 2) +
  scale_fill_viridis_c(option = "B") +
  xlab("snRNA, Inh MGE SST, MAST model log2FC") +
  ylab("FANS-sorted bulk RNA, Neurons, edgeR log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/snRNA_level2_InhMGESst_FC_vs_FANSsorted_neuron_FC.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)
