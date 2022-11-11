library(tidyverse)
library(ggrastr)

DE_level_2 <- read_tsv("MAST_res_level2_splitControls_summary.txt")

FDR_threshold <- 0.05
FC_threshold <- 1.2

cell_types <- c(
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

res <- DE_level_2 %>%
  mutate(
    sig_DE = if_else(
      (fdr < FDR_threshold &
        abs(model_log2FC) > log2(FC_threshold) &
        conv_C == TRUE &
        conv_D == TRUE &
        ci.hi * ci.lo > 0 &
        abs(model_log2FC - avg_log2FC) < 2
      ), "Sig", "nonSig"
    )
  ) %>%
  select(gene, model_log2FC, sig_DE, cell_type, region, diagnosis_1) %>%
  pivot_wider(names_from = "diagnosis_1", values_from = c("model_log2FC", "sig_DE")) %>%
  mutate(
    sig = case_when(
      sig_DE_ALS == "Sig" & sig_DE_FTD == "Sig" ~ "Both",
      sig_DE_ALS == "Sig" & (sig_DE_FTD == "nonSig" | is.na(sig_DE_FTD)) ~ "ALS only",
      (sig_DE_ALS == "nonSig" | is.na(sig_DE_ALS)) & sig_DE_FTD == "Sig" ~ "FTD only",
      (sig_DE_ALS == "nonSig" | is.na(sig_DE_ALS)) & (sig_DE_FTD == "nonSig" | is.na(sig_DE_FTD)) ~ "Neither"
    )
  )

disease_color <- c(
  "ALS only" = "#F6416C",
  "FTD only" = "#FFDE7D",
  "Both" = "#ff955a",
  "Neither" = "darkgrey"
)

p <- res %>%
  mutate(
    sig = factor(sig, levels = c("ALS only", "FTD only", "Both", "Neither")),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  sample_frac() %>%
  ggplot(aes(model_log2FC_ALS, model_log2FC_FTD))
p +
  geom_point_rast(
    aes(color = sig),
    raster.dpi = 300,
    scale = 0.5
  ) +
  facet_wrap(~ cell_type + region) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  "./DE_disease_effect_FC_in_ALS_and_FTD_scatterPlot.pdf",
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)


# correlation of FC
get_cor <- function(selected_cell_type, selected_region) {
  res_sub <- res %>%
    filter(cell_type == selected_cell_type, region == selected_region)
  test_res <- cor.test(
    res_sub$model_log2FC_ALS,
    res_sub$model_log2FC_FTD,
    method = "pearson",
    conf.level = 0.95
  )
  tibble(
    cell_type = selected_cell_type,
    region = selected_region,
    pearson_cor = test_res$estimate,
    cor_conf_low = test_res$conf.int[1],
    cor_conf_hi = test_res$conf.int[2],
    p = test_res$p.value
  )
}

params <- expand_grid(
  selected_cell_type = cell_types,
  selected_region = c("MCX", "mFCX")
)
cor_res <- pmap_dfr(params, get_cor)

region_color <- c("MCX" = "#432266", "mFCX" = "#FAA51B")

p <- cor_res %>%
  mutate(
    fdr = p.adjust(p, method = "fdr"),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  ggplot(aes(cell_type, pearson_cor))
p +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_errorbar(
    aes(ymin = cor_conf_low, ymax = cor_conf_hi, color = region),
    size = 0.75,
    width = 0.25
  ) +
  geom_point(aes(color = region), size = 1) +
  scale_color_manual(values = region_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  "./DE_disease_effect_FC_pearsonCor_in_ALS_and_FTD_using_all_expressed_genes.pdf",
  device = cairo_pdf(),
  width = 2.5,
  height = 2,
  useDingbats = FALSE
)


get_cor_onlySigBoth <- function(selected_cell_type, selected_region) {
  res_sub <- res %>%
    filter(sig == "Both") %>%
    filter(cell_type == selected_cell_type, region == selected_region)
  test_res <- cor.test(
    res_sub$model_log2FC_ALS,
    res_sub$model_log2FC_FTD,
    method = "pearson",
    conf.level = 0.95
  )
  tibble(
    cell_type = selected_cell_type,
    region = selected_region,
    pearson_cor = test_res$estimate,
    cor_conf_low = test_res$conf.int[1],
    cor_conf_hi = test_res$conf.int[2],
    p = test_res$p.value
  )
}

params2 <- expand_grid(
  selected_cell_type = c("Astro", "Oligo", "OPC"),
  selected_region = c("MCX", "mFCX")
)
cor_res2 <- pmap_dfr(params2, get_cor_onlySigBoth)

p <- cor_res2 %>%
  mutate(
    fdr = p.adjust(p, method = "fdr"),
    cell_type = factor(cell_type, levels = c("Astro", "Oligo", "OPC"))
  ) %>%
  ggplot(aes(cell_type, pearson_cor))
p +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_errorbar(
    aes(ymin = cor_conf_low, ymax = cor_conf_hi, color = region),
    size = 0.75,
    width = 0.25
  ) +
  geom_point(aes(color = region), size = 1) +
  scale_color_manual(values = region_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  "./DE_disease_effect_FC_pearsonCor_in_ALS_and_FTD_using_sigDE_genes_in_both.pdf",
  device = cairo_pdf(),
  width = 2.5,
  height = 2,
  useDingbats = FALSE
)

# compare with the original DE results
ori_DE_level_2 <- read_tsv("./MAST_res_level2_summary.txt")

ori_DE_level_2_sub <- ori_DE_level_2 %>%
  select(cond_1, cond_2, cell_type, region, gene, model_log2FC, FDR) %>%
  rename(
    diagnosis_1 = cond_1,
    diagnosis_2 = cond_2,
    model_log2FC_ori = model_log2FC,
    FDR_ori = FDR
  )

DE_level_2_sub <- DE_level_2 %>%
  select(diagnosis_1, diagnosis_2, cell_type, region, gene, model_log2FC, fdr)

merged_res <- DE_level_2_sub %>%
  left_join(ori_DE_level_2_sub)

p <- merged_res %>% ggplot(aes(model_log2FC, model_log2FC_ori))
p +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "darkgrey") +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "D", trans = "log10") +
  facet_wrap(~ cell_type + region + diagnosis_1) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  "./DE_disease_effect_FC_in_splitControls_vs_originalFullControls.pdf",
  device = cairo_pdf(),
  width = 8,
  height = 6,
  useDingbats = FALSE
)

# show the counts of DE genes in each group
DE_count <- res %>%
  count(cell_type, region, sig) %>%
  complete(cell_type, region, sig, fill = list(n = 0)) %>%
  mutate(
    sig = factor(sig, levels = c("ALS only", "FTD only", "Both", "Neither")),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  )

disease_color <- c(
  "ALS only" = "#F6416C",
  "FTD only" = "#FFDE7D",
  "Both" = "#ff955a",
  "Neither" = "darkgrey"
)

p <- DE_count %>% ggplot(aes(cell_type, n))

p +
  geom_bar(
    aes(fill = sig),
    stat = "identity",
    position = position_dodge2(preserve = "single"),
    width = 0.75
  ) +
  facet_wrap(~region) +
  xlab("Cell type") +
  ylab("Number of genes") +
  scale_fill_manual(name = "Sig. DE in", values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  "./num_sig_disease_DE_genes_in_splitControls_ALS_and_FTD.pdf",
  device = cairo_pdf(),
  width = 4,
  height = 2,
  useDingbats = FALSE
)


p <- DE_count %>%
  filter(sig != "Neither") %>%
  ggplot(aes(cell_type, n))

p +
  geom_bar(
    aes(fill = sig),
    stat = "identity",
    position = position_dodge2(preserve = "single"),
    width = 0.75
  ) +
  facet_wrap(~region) +
  xlab("Cell type") +
  ylab("Number of genes") +
  scale_fill_manual(name = "Sig. DE in", values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  "./num_sig_disease_DE_genes_in_splitControls_ALS_and_FTD_removeGroupNeither.pdf",
  device = cairo_pdf(),
  width = 4,
  height = 2,
  useDingbats = FALSE
)

# session info
sessionInfo()
