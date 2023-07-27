library(tidyverse)
library(ggrastr)
library(ggrepel)
library(shades)
library(plotly)
library(htmlwidgets)
library(writexl)

DE_level_2 <- read_tsv("MAST_res_level2_splitControls_summary.txt")
# write_xlsx(DE_level_2, "MAST_res_level2_splitControls_summary.xlsx")

full_DE_level_2 <- read_tsv("./MAST_res_level2_summary.txt")

FDR_threshold <- 0.05
FC_threshold <- 1.2

cell_types <- c(
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

full_DE_level_2_sub <- full_DE_level_2 %>%
  select(-num_cells_cond_1, -num_cells_cond_2) %>%
  filter(
    cell_type %in% cell_types
  ) %>%
  rename_with(function(x) paste0("full_model_", x), p_value:avg_logFC)

merged_DE <- DE_level_2 %>%
  rename_with(function(x) paste0("split_control_", x), p_value:avg_log2FC) %>%
  rename(
    cond_1 = diagnosis_1,
    cond_2 = diagnosis_2
  ) %>%
  left_join(full_DE_level_2_sub, by = c("gene", "cell_type", "region", "cond_1", "cond_2")) %>%
  mutate(
    sig_DE_full = if_else(
      (full_model_FDR < FDR_threshold &
        abs(full_model_model_log2FC) > log2(FC_threshold) &
        full_model_conv_C == TRUE &
        full_model_conv_D == TRUE &
        full_model_model_log2FC_ci_low * full_model_model_log2FC_ci_hi > 0 &
        abs(full_model_model_log2FC - full_model_avg_logFC) < 2
      ), "Sig", "nonSig"
    ),
    sig_DE_split_control = if_else(
      (split_control_fdr < FDR_threshold &
        abs(split_control_model_log2FC) > log2(FC_threshold) &
        split_control_conv_C == TRUE &
        split_control_conv_D == TRUE &
        split_control_ci.hi * split_control_ci.lo > 0 &
        abs(split_control_model_log2FC - split_control_avg_log2FC) < 2
      ), "Sig", "nonSig"
    )
  ) %>%
  select(cond_1, cond_2, region, everything())

# write_xlsx(merged_DE, "MAST_res_level2_merged_full_model_and_split_control_model_summary.xlsx")


merged_res <- merged_DE %>%
  select(
    cond_1, region, cell_type, gene,
    full_model_model_log2FC, split_control_model_log2FC,
    sig_DE_full, sig_DE_split_control
  ) %>%
  pivot_wider(
    names_from = "cond_1",
    values_from = c("full_model_model_log2FC", "split_control_model_log2FC", "sig_DE_full", "sig_DE_split_control")
  )
# write_xlsx(merged_res, "MAST_res_level2_merged_full_model_and_split_control_model_significance.xlsx")

merged_res_sub <- merged_res %>%
  filter(
    sig_DE_full_ALS == "Sig" & sig_DE_full_FTD == "Sig" &
      sig_DE_split_control_ALS == "Sig" & sig_DE_split_control_FTD == "Sig"
  ) %>%
  filter(
    sign(full_model_model_log2FC_ALS) != sign(full_model_model_log2FC_FTD),
    sign(split_control_model_log2FC_ALS) != sign(split_control_model_log2FC_FTD),
    sign(full_model_model_log2FC_ALS) == sign(split_control_model_log2FC_ALS),
    sign(full_model_model_log2FC_FTD) == sign(split_control_model_log2FC_FTD)
  )
# write_xlsx(merged_res_sub, "MAST_res_level2_merged_full_model_and_split_control_model_significance_opposite_DE_in_two_diseases.xlsx")


res <- DE_level_2 %>%
  mutate(
    sig_DE = if_else(
      (fdr < FDR_threshold &
        abs(model_log2FC) > log2(FC_threshold) &
        conv_C == TRUE &
        conv_D == TRUE &
        ci.hi * ci.lo > 0 &
        abs(model_log2FC - avg_log2FC) < 2
      # ci.hi * ci.lo > 0
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

write_tsv(res, "DE_disease_effect_FC_in_ALS_and_FTD.txt")
write_xlsx(res, "DE_disease_effect_FC_in_ALS_and_FTD.xlsx")

disease_color <- c(
  "ALS only" = "#F6416C",
  "FTD only" = "#FFDE7D",
  "Both" = "#ff955a",
  "Neither" = "darkgrey"
)

disease_color_dark <- lightness(disease_color, scalefac(0.9)) %>%
  as.vector()
names(disease_color_dark) <- names(disease_color)

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
dev.off()

# Astro
p <- res %>%
  mutate(
    sig = factor(sig, levels = c("ALS only", "FTD only", "Both", "Neither")),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  filter(sig != "Neither") %>%
  filter(cell_type == "Astro") %>%
  sample_frac() %>%
  ggplot(aes(model_log2FC_ALS, model_log2FC_FTD))
p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point_rast(
    aes(color = sig),
    raster.dpi = 300,
    scale = 0.4
  ) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_Astro_scatterPlot.pdf",
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
dev.off()

p2 <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(size = 1, aes(label = gene, color = sig)) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p2 <- ggplotly(p2)
saveWidget(
  p2,
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_Astro_scatterPlot.html",
  selfcontained = TRUE
)

selected_genes <- c("GFAP", "CD44", "ARHGEF4", "SLC7A11", "OSBPL3", "MIB1", "FLRT2", "RFX2")
df_sub <- res %>%
  mutate(
    sig = factor(sig, levels = c("ALS only", "FTD only", "Both", "Neither")),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  filter(sig != "Neither") %>%
  filter(cell_type == "Astro") %>%
  sample_frac()

p <- df_sub %>%
  ggplot(aes(model_log2FC_ALS, model_log2FC_FTD))
p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point_rast(
    aes(color = sig),
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    # aes(color = sig),
    color = "black",
    data = df_sub %>% filter(gene %in% selected_genes),
    size = 2,
    shape = 1
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% selected_genes),
    size = 5,
    nudge_x = 0.5,
    nudge_y = 0.5
  ) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# Oligo
p <- res %>%
  mutate(
    sig = factor(sig, levels = c("ALS only", "FTD only", "Both", "Neither")),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  filter(sig != "Neither") %>%
  filter(cell_type == "Oligo") %>%
  sample_frac() %>%
  ggplot(aes(model_log2FC_ALS, model_log2FC_FTD))
p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point_rast(
    aes(color = sig),
    raster.dpi = 300,
    scale = 0.4
  ) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_Oligo_scatterPlot.pdf",
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
dev.off()

p2 <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(size = 1, aes(label = gene, color = sig)) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p2 <- ggplotly(p2)
saveWidget(
  p2,
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_Oligo_scatterPlot.html",
  selfcontained = TRUE
)

selected_genes <- c(
  #   "STRN", "LDLRAD4", "ZKSCAN1",
  # "AGPAT4", "GAS7", "FGFR2", "ZFYVE16", "SLC4A8",
  "MBP", "MOG", "OPALIN"
)
df_sub <- res %>%
  mutate(
    sig = factor(sig, levels = c("ALS only", "FTD only", "Both", "Neither")),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  filter(sig != "Neither") %>%
  filter(cell_type == "Oligo") %>%
  sample_frac()

p <- df_sub %>%
  ggplot(aes(model_log2FC_ALS, model_log2FC_FTD))
p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point_rast(
    aes(color = sig),
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    # aes(color = sig),
    color = "black",
    data = df_sub %>% filter(gene %in% selected_genes),
    size = 2,
    shape = 1
  ) +
  geom_text_repel(
    aes(label = gene),
    data = df_sub %>% filter(gene %in% selected_genes),
    size = 5,
    nudge_x = 0.5,
    nudge_y = 0.5
  ) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# Micro
p <- res %>%
  mutate(
    sig = factor(sig, levels = c("ALS only", "FTD only", "Both", "Neither")),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  filter(sig != "Neither") %>%
  filter(cell_type == "Micro") %>%
  sample_frac() %>%
  ggplot(aes(model_log2FC_ALS, model_log2FC_FTD))
p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point_rast(
    aes(color = sig),
    raster.dpi = 300,
    scale = 0.4
  ) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_Micro_scatterPlot.pdf",
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
dev.off()

p2 <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(size = 1, aes(label = gene, color = sig)) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p2 <- ggplotly(p2)
saveWidget(
  p2,
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_Micro_scatterPlot.html",
  selfcontained = TRUE
)

# OPC
p <- res %>%
  mutate(
    sig = factor(sig, levels = c("ALS only", "FTD only", "Both", "Neither")),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  filter(sig != "Neither") %>%
  filter(cell_type == "OPC") %>%
  sample_frac() %>%
  ggplot(aes(model_log2FC_ALS, model_log2FC_FTD))
p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point_rast(
    aes(color = sig),
    raster.dpi = 300,
    scale = 0.4
  ) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_OPC_scatterPlot.pdf",
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
dev.off()

p2 <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(size = 1, aes(label = gene, color = sig)) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p2 <- ggplotly(p2)
saveWidget(
  p2,
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_OPC_scatterPlot.html",
  selfcontained = TRUE
)

# Endo
p <- res %>%
  mutate(
    sig = factor(sig, levels = c("ALS only", "FTD only", "Both", "Neither")),
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  filter(sig != "Neither") %>%
  filter(cell_type == "Endo") %>%
  sample_frac() %>%
  ggplot(aes(model_log2FC_ALS, model_log2FC_FTD))
p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point_rast(
    aes(color = sig),
    raster.dpi = 300,
    scale = 0.4
  ) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_Endo_scatterPlot.pdf",
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
dev.off()

p2 <- p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(size = 1, aes(label = gene, color = sig)) +
  facet_wrap(~ cell_type + region + sig) +
  scale_color_manual(values = disease_color_dark) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p2 <- ggplotly(p2)
saveWidget(
  p2,
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_Endo_scatterPlot.html",
  selfcontained = TRUE
)


# plot selected genes that are significant in both ALS and FTD;
# also significant in both the full model and the split control model
# in Astro, Oligo and OPC
cell_type_colors <- c(
  "Astro" = "#ef6075",
  "Oligo" = "#89520d",
  "OPC" = "#b98735"
)
res_sub <- res %>%
  mutate(
    cell_type = factor(cell_type, levels = c("Astro", "Micro", "Oligo", "OPC", "Endo", "VLMC"))
  ) %>%
  filter(sig == "Both") %>%
  filter(cell_type %in% c("Astro", "Oligo", "OPC"))

# astro_highlight_genes <- c(
#   "MAPK4", "KCNN3", "GABBR1", "SPON1",
#   "NEDD4L", "ARHGAP42", "LGI1", "ZFHX4", "EFNA5", "MYO18A",
#   "GABBR1", "GFAP"
# )
# oligo_highlight_genes <- c(
#   "EFHD1", "ENOX1", "TRA2A",
#   "PHLPP1", "NAV2", "GAS7", "JAM3", "STOX2", "PSD3", "ZNF280D",
#   "RASGRF2", "RASGRF1"
# )
# opc_highlight_genes <- c(
#   "MAML2", "PRRX1", "AFAP1L2",
#   "GRIA4", "KCND2", "KCNMA1", "RIN2", "GNPTAB"
# )

astro_highlight_genes <- c(
  "LGI1", "ZFHX4", "EFNA5", "MYO18A",
  "GABBR1", "GFAP", "TMEM132B"
)

oligo_highlight_genes <- c(
  "NAV2",
  "RASGRF2", "RASGRF1", "CNTNAP2", "LPAR1", "MAPK10", "PLCG2"
)
opc_highlight_genes <- c(
  "RIT2", "LINC00854",
  "GRIA4", "KCND2", "GNPTAB"
)


p <- res_sub %>%
  ggplot(aes(model_log2FC_ALS, model_log2FC_FTD))
p +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point_rast(
    aes(color = cell_type),
    raster.dpi = 300,
    scale = 0.4
  ) +
  geom_point(
    aes(fill = cell_type),
    color = "black",
    data = res_sub %>%
      filter(gene %in% astro_highlight_genes, cell_type == "Astro"),
    shape = 21,
    size = 2
  ) +
  geom_point(
    aes(fill = cell_type),
    color = "black",
    data = res_sub %>%
      filter(gene %in% oligo_highlight_genes, cell_type == "Oligo"),
    shape = 21,
    size = 2
  ) +
  geom_point(
    aes(fill = cell_type),
    color = "black",
    data = res_sub %>%
      filter(gene %in% opc_highlight_genes, cell_type == "OPC"),
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    data = res_sub %>%
      filter(gene %in% astro_highlight_genes, cell_type == "Astro"),
    aes(label = gene),
    size = 2
  ) +
  geom_text_repel(
    data = res_sub %>%
      filter(gene %in% oligo_highlight_genes, cell_type == "Oligo"),
    aes(label = gene),
    size = 2
  ) +
  geom_text_repel(
    data = res_sub %>%
      filter(gene %in% opc_highlight_genes, cell_type == "OPC"),
    aes(label = gene),
    size = 2
  ) +
  # facet_wrap(~ cell_type + region, ncol = 2) +
  facet_wrap(~ cell_type + region, ncol = 2, scale = "free") +
  xlab("log2FC in ALS vs. group 1 Control") +
  ylab("log2FC in FTD vs. group 2 Control") +
  scale_color_manual(values = cell_type_colors) +
  scale_fill_manual(values = cell_type_colors) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

ggsave(
  "./plot_splitControls/DE_disease_effect_FC_in_ALS_and_FTD_Astro_scatterPlot_highlightGenes.pdf",
  device = cairo_pdf(),
  width = 2.8,
  height = 5.3,
  useDingbats = FALSE
)
dev.off()


#####################################################################################

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
