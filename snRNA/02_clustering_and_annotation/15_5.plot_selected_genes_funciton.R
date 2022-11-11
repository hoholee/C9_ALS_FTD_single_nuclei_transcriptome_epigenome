# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(scales)
library(scico)
library(tictoc)
library(shades)

# read in annotated data
# data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_clean_SeuratV4_object.rds")
# meta_data <- data_obj_sub@meta.data

# load the gene x cell CPM matrix
cpm_mat <- readRDS("snRNA_geneByCell_dgCMatrix_RNA_CPM.rds")

# read metadata
meta_data <- read_tsv("metadata_all_cells_2nd_round_annotations.txt", col_types = "cciiccccddcddliicccccccc") %>%
  filter(rna_anno_2ndRound_level_3 != "Ambiguous")

# check order
all.equal(colnames(cpm_mat), meta_data$cell_id)

disease_color_panel <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")

# read DE list
DE_level2 <- read_tsv("./DE_MAST/second_round_level_2/MAST_res_level2_summary.txt")

# define function for split violin plot
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    # Original function by Jan Gleixner (@jan-glx)
    # Adjustments by Wouter van der Bijl (@Axeman)
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...)
  )
}

# define function to plot chosen genes in chosen cell types
plot_snRNA_exp_violin <- function(selected_region = "MCX",
                                  selected_disease_1 = "ALS",
                                  selected_disease_2 = "Control",
                                  selected_cell_types = c("Astro"),
                                  selected_genes = c("CD44"),
                                  FDR_threshold = 0.05,
                                  FC_threshold = 1.2,
                                  facet_n_col = 1,
                                  strip_position = "left") {
  selected_diseases <- c(selected_disease_1, selected_disease_2)
  cell_idx <- which(meta_data$region == selected_region &
    meta_data$disease %in% selected_diseases &
    meta_data$rna_anno_2ndRound_level_2 %in% selected_cell_types)

  meta_data_used <- meta_data %>%
    filter(
      region == selected_region,
      disease %in% selected_diseases,
      rna_anno_2ndRound_level_2 %in% selected_cell_types
    )

  cpm_mat_used <- cpm_mat[selected_genes, cell_idx, drop = FALSE]
  all.equal(colnames(cpm_mat_used), meta_data_used$cell_id)

  DE_used <- DE_level2 %>%
    filter(
      cond_1 == selected_disease_1,
      cond_2 == selected_disease_2,
      region == selected_region,
      cell_type %in% selected_cell_types,
      gene %in% selected_genes
    ) %>%
    mutate(
      sig = if_else(
        (FDR < FDR_threshold &
          abs(model_log2FC) > log2(FC_threshold) &
          conv_C == TRUE &
          conv_D == TRUE &
          model_log2FC_ci_hi * model_log2FC_ci_low > 0 &
          abs(model_log2FC - avg_logFC) < 2
        ),
        "DE",
        "nonDE"
      ),
      sig_sign = case_when(
        sig == "DE" & FDR < 0.05 & FDR >= 0.01 ~ "*",
        sig == "DE" & FDR < 0.01 & FDR >= 0.001 ~ "**",
        sig == "DE" & FDR < 0.001 & FDR >= 0.0001 ~ "***",
        sig == "DE" & FDR < 0.0001 ~ "****",
        TRUE ~ NA_character_
      ),
      DE_dir = if_else(model_log2FC > 0, "up-regulated", "down-regulated"),
      gene = factor(gene, levels = selected_genes)
    )

  df_genes_cpm <- cpm_mat_used %>%
    as.matrix() %>%
    t() %>%
    as_tibble() %>%
    rename_all(.funs = function(x) (paste0("gene_CPM_", x))) %>%
    bind_cols(meta_data_used)

  to_plot_df <- df_genes_cpm %>%
    select(
      starts_with("gene_CPM_"),
      disease,
      subject,
      rna_anno_2ndRound_level_2,
    ) %>%
    pivot_longer(
      cols = starts_with("gene_CPM_"),
      names_to = "gene",
      values_to = "CPM"
    ) %>%
    mutate(
      gene = str_replace(gene, "^gene_CPM_", ""),
      log10_CPM_1p = log10(CPM + 1)
    ) %>%
    left_join(
      DE_used,
      by = c("gene" = "gene", "rna_anno_2ndRound_level_2" = "cell_type")
    ) %>%
    mutate(
      gene = factor(gene, levels = selected_genes),
      rna_anno_2ndRound_level_2 = factor(
        rna_anno_2ndRound_level_2,
        levels = selected_cell_types
      ),
    ) %>%
    arrange(gene)

  p <- to_plot_df %>% ggplot(aes(rna_anno_2ndRound_level_2, log10_CPM_1p))

  p + geom_split_violin(
    aes(
      fill = disease,
      group = paste(rna_anno_2ndRound_level_2, disease)
    ),
    draw_quantiles = c(0.5),
    trim = TRUE,
    scale = "width",
    size = 0.2,
    width = 0.8
  ) +
    geom_text(
      aes(
        x = cell_type,
        y = 4,
        label = sig_sign,
        color = DE_dir
      ),
      data = DE_used,
      size = 2
    ) +
    facet_wrap(
      ~gene,
      ncol = facet_n_col,
      scales = "fixed",
      strip.position = strip_position
    ) +
    ggtitle(selected_region) +
    scale_fill_manual(
      name = "diagnosis",
      values = disease_color_panel,
      limits = force
    ) +
    scale_color_manual(
      name = "sig. DE",
      values = c(
        "up-regulated" = as.vector(disease_color_panel[selected_disease_1]),
        "down-regulated" = as.vector(disease_color_panel[selected_disease_2])
      )
    ) +
    xlab("Major cell types") +
    ylab("snRNA gene expression (log10(CPM+1))") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}

plot_snRNA_exp_violin(
  selected_region = "MCX",
  selected_disease_1 = "ALS",
  selected_disease_2 = "Control",
  selected_cell_types = c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep",
    "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST",
    "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
  ),
  selected_genes = c(
    "GFAP", "CD44", "CHI3L1", "RANBP3L",
    "MAOB", "TGFB2", "ATP2B4", "GOLGB1"
  ),
  FDR_threshold = 0.05,
  FC_threshold = 1.2,
  facet_n_col = 1,
  strip_position = "left"
)

ggsave(
  filename = "./plots/DE_ALS_vs_Control_MCX_example_genes_Astro_split.pdf",
  device = cairo_pdf(),
  width = 4,
  height = 4,
  useDingbats = FALSE
)

plot_snRNA_exp_violin(
  selected_region = "mFCX",
  selected_disease_1 = "ALS",
  selected_disease_2 = "Control",
  selected_cell_types = c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep",
    "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST",
    "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
  ),
  selected_genes = c(
    "GFAP", "CD44", "CHI3L1", "RANBP3L",
    "MAOB", "TGFB2", "ATP2B4", "GOLGB1"
  ),
  FDR_threshold = 0.05,
  FC_threshold = 1.2,
  facet_n_col = 1,
  strip_position = "left"
)

ggsave(
  filename = "./plots/DE_ALS_vs_Control_mFCX_example_genes_Astro_split.pdf",
  device = cairo_pdf(),
  width = 4,
  height = 4,
  useDingbats = FALSE
)

# Exc-neurons
# exc_genes <- c(
#   "NDUFS7", "NDUFC1", "NDUFB2", "NDUFAF6", "NDUFAF5",
#   "NDUFAB1", "NDUFA5", "NDUFA4", "NDUFA12", "NDUFA1"
# )

# exc_genes <- c(
#   "SDHC", "UQCRH", "UQCRFS1", "UQCRB",
#   "CYCS", "CYC1", "COX7B", "COX7A2", "COX4I1"
# )

# exc_genes <- c(
#   "ATP5PF", "ATP5PB", "ATP5MPL", "ATP5MD", "ATP5IF1",
#   "ATP5F1E", "ATP5F1C", "ATP5F1B", "ATP5F1A"
# )

# exc_genes <- c(
#   "SLC25A28", "SLC25A3", "SLC25A32", "SLC25A37", "SLC25A40",
#   "SLC25A4", "TIMMDC1", "TIMM23B", "SOD1", "PRDX1"
# )

# exc_genes <- c(
#   "HSPA9", "HSPA8", "HSPA4", "HSPB1", "HSPB11",
#   "HSP90AA1", "HSP90AB1", "CLU"
# )

# exc_genes <- c(
#   "DNAJC9", "DNAJC7", "DNAJC3-DT", "DNAJC11", "DNAJC1",
#   "DNAJB1", "DNAJA1", "DNAJA2", "DNAJA3"
# )

# exc_genes <- c(
#   "RPS29", "RPS23", "RPS28",
#   "MRPS31", "MRPS30", "MRPL50", "MRPL20"
# )

# exc_genes <- c(
#   "SRP14", "EIF2S3", "EIF2B1", "EIF1B"
# )

# exc_genes <- c(
#   "NUP107", "NUP205", "TPR", "IPO5", "IPO7",
#   "XPO1", "XPOT", "KPNA5", "RANGAP1", "RANBP2",
#   "RGPD2", "ENY2", "RPA1", "RPA3", "NEK1"
# )

# exc_genes <- c(
#   "NDUFA12", "NDUFA4", "UQCRH", "COX7A2", "ATP5PF",
#   "ATP5F1A", "SLC25A40", "SOD1", "PRDX1", "HSP90AA1",
#   "HSP90AB1", "CLU", "HSPB1", "EIF1B",
#   "DNAJB1", "DNAJA1", "DNAJA2", "DNAJA3", "RPS29", "KCND3",
# )

exc_genes <- c(
  "NDUFA12", "COX7A2", "ATP5PF", "SOD1", "PRDX1",
  "HSP90AA1", "CLU", "EIF1B", "RPS29", "KCND3"
)

plot_snRNA_exp_violin(
  selected_region = "MCX",
  selected_disease_1 = "ALS",
  selected_disease_2 = "Control",
  selected_cell_types = c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep"
    # "Exc_superficial", "Exc_intermediate", "Exc_deep",
    # "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST",
    # "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
  ),
  selected_genes = exc_genes,
  FDR_threshold = 0.05,
  FC_threshold = 1.2,
  facet_n_col = 1,
  strip_position = "left"
)

ggsave(
  filename = "./plots/DE_ALS_vs_Control_MCX_example_genes_ExcNeu_split.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 4,
  useDingbats = FALSE
)

plot_snRNA_exp_violin(
  selected_region = "mFCX",
  selected_disease_1 = "ALS",
  selected_disease_2 = "Control",
  selected_cell_types = c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep"
    # "Exc_superficial", "Exc_intermediate", "Exc_deep",
    # "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST",
    # "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
  ),
  selected_genes = exc_genes,
  FDR_threshold = 0.05,
  FC_threshold = 1.2,
  facet_n_col = 1,
  strip_position = "left"
)

ggsave(
  filename = "./plots/DE_ALS_vs_Control_mFCX_example_genes_ExcNeu_split.pdf",
  device = cairo_pdf(),
  width = 2,
  height = 4,
  useDingbats = FALSE
)

exc_genes <- c(
  "NDUFA12", "COX7A2", "ATP5PF", "SOD1", "PRDX1",
  "HSP90AA1", "CLU", "EIF1B", "RPS29", "KCND3", "SLIT3", "DNMT3A"
)

plot_snRNA_exp_violin(
  selected_region = "MCX",
  selected_disease_1 = "ALS",
  selected_disease_2 = "Control",
  selected_cell_types = c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep"
  ),
  selected_genes = exc_genes,
  FDR_threshold = 0.05,
  FC_threshold = 1.2,
  facet_n_col = length(exc_genes),
  strip_position = "top"
)

ggsave(
  filename = "./plots/DE_ALS_vs_Control_MCX_example_genes_ExcNeu_split_horizontal.pdf",
  device = cairo_pdf(),
  width = 6.5,
  height = 2,
  useDingbats = FALSE
)

plot_snRNA_exp_violin(
  selected_region = "mFCX",
  selected_disease_1 = "ALS",
  selected_disease_2 = "Control",
  selected_cell_types = c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep"
  ),
  selected_genes = exc_genes,
  FDR_threshold = 0.05,
  FC_threshold = 1.2,
  facet_n_col = length(exc_genes),
  strip_position = "top"
)

ggsave(
  filename = "./plots/DE_ALS_vs_Control_mFCX_example_genes_ExcNeu_split_horizontal.pdf",
  device = cairo_pdf(),
  width = 6.5,
  height = 2,
  useDingbats = FALSE
)

# plot C9orf72

# plot_snRNA_exp_violin(
#   selected_region = "MCX",
#   selected_disease_1 = "ALS",
#   selected_disease_2 = "Control",
#   selected_cell_types = c(
#     "Exc_superficial", "Exc_intermediate", "Exc_deep",
#     "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST",
#     "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
#   ),
#   selected_genes = c(
#     "C9orf72"
#   ),
#   FDR_threshold = 0.05,
#   FC_threshold = 1.2,
#   facet_n_col = 1
# )

# plot_snRNA_exp_violin(
#   selected_region = "mFCX",
#   selected_disease_1 = "ALS",
#   selected_disease_2 = "Control",
#   selected_cell_types = c(
#     "Exc_superficial", "Exc_intermediate", "Exc_deep",
#     "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST",
#     "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
#   ),
#   selected_genes = c(
#     "C9orf72"
#   ),
#   FDR_threshold = 0.05,
#   FC_threshold = 1.2,
#   facet_n_col = 1
# )

# plot_snRNA_exp_violin(
#   selected_region = "MCX",
#   selected_disease_1 = "FTD",
#   selected_disease_2 = "Control",
#   selected_cell_types = c(
#     "Exc_superficial", "Exc_intermediate", "Exc_deep",
#     "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST",
#     "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
#   ),
#   selected_genes = c(
#     "C9orf72"
#   ),
#   FDR_threshold = 0.05,
#   FC_threshold = 1.2,
#   facet_n_col = 1
# )

# plot_snRNA_exp_violin(
#   selected_region = "mFCX",
#   selected_disease_1 = "FTD",
#   selected_disease_2 = "Control",
#   selected_cell_types = c(
#     "Exc_superficial", "Exc_intermediate", "Exc_deep",
#     "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST",
#     "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
#   ),
#   selected_genes = c(
#     "C9orf72"
#   ),
#   FDR_threshold = 0.05,
#   FC_threshold = 1.2,
#   facet_n_col = 1
# )

selected_cell_types2 <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)
cell_idx2 <- which(meta_data$rna_anno_2ndRound_level_2 %in% selected_cell_types2)

meta_data_used2 <- meta_data %>%
  filter(
    rna_anno_2ndRound_level_2 %in% selected_cell_types2
  )

cpm_mat_C9orf72 <- cpm_mat["C9orf72", cell_idx2, drop = FALSE]
all.equal(colnames(cpm_mat_C9orf72), meta_data_used2$cell_id)

DE_C9orf72 <- DE_level2 %>%
  filter(
    gene == "C9orf72",
    cell_type %in% selected_cell_types2
  ) %>%
  mutate(
    sig = if_else(
      (FDR < 0.05 &
        abs(model_log2FC) > log2(1.2) &
        conv_C == TRUE &
        conv_D == TRUE &
        model_log2FC_ci_hi * model_log2FC_ci_low > 0 &
        abs(model_log2FC - avg_logFC) < 2
      ),
      "DE",
      "nonDE"
    ),
    sig_sign = case_when(
      sig == "DE" & FDR < 0.05 & FDR >= 0.01 ~ "*",
      sig == "DE" & FDR < 0.01 & FDR >= 0.001 ~ "**",
      sig == "DE" & FDR < 0.001 & FDR >= 0.0001 ~ "***",
      sig == "DE" & FDR < 0.0001 ~ "****",
      TRUE ~ NA_character_
    ),
    DE_dir = if_else(model_log2FC > 0, "up-regulated", "down-regulated")
  )

df_C9orf72_cpm <- cpm_mat_C9orf72 %>%
  as.matrix() %>%
  t() %>%
  as_tibble() %>%
  rename_all(.funs = function(x) (paste0("gene_CPM_", x))) %>%
  bind_cols(meta_data_used2)

to_plot_C9orf72 <- df_C9orf72_cpm %>%
  select(
    starts_with("gene_CPM_"),
    region,
    disease,
    subject,
    rna_anno_2ndRound_level_2,
  ) %>%
  pivot_longer(
    cols = starts_with("gene_CPM_"),
    names_to = "gene",
    values_to = "CPM"
  ) %>%
  mutate(
    gene = str_replace(gene, "^gene_CPM_", ""),
    log10_CPM_1p = log10(CPM + 1)
  ) %>%
  mutate(
    rna_anno_2ndRound_level_2 = factor(
      rna_anno_2ndRound_level_2,
      levels = selected_cell_types2
    ),
    disease = factor(disease, levels = c("ALS", "FTD", "Control"))
  ) %>%
  mutate(
    group = paste(rna_anno_2ndRound_level_2, disease),
    group = factor(
      group,
      levels = paste(
        rep(selected_cell_types2, each = 3),
        c("ALS", "FTD", "Control")
      )
    )
  ) %>%
  arrange(gene)

p <- to_plot_C9orf72 %>%
  ggplot(aes(rna_anno_2ndRound_level_2, log10_CPM_1p))

p +
  geom_violin(
    aes(
      fill = disease,
      group = group
    ),
    draw_quantiles = c(0.5),
    trim = TRUE,
    scale = "width",
    size = 0.2,
    width = 0.8
  ) +
  facet_wrap(
    ~region,
    ncol = 1,
    scales = "fixed",
    strip.position = "left"
  ) +
  scale_fill_manual(
    name = "diagnosis",
    values = disease_color_panel,
    limits = force
  ) +
  xlab("Major cell types") +
  ylab("snRNA gene expression (log10(CPM+1))") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none"
  )

ggsave(
  filename = "./plots/DE_violin_C9orf72.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 2.5,
  useDingbats = FALSE
)

# remove FTD

to_plot_C9orf72 <- df_C9orf72_cpm %>%
  select(
    starts_with("gene_CPM_"),
    region,
    disease,
    subject,
    rna_anno_2ndRound_level_2,
  ) %>%
  pivot_longer(
    cols = starts_with("gene_CPM_"),
    names_to = "gene",
    values_to = "CPM"
  ) %>%
  mutate(
    gene = str_replace(gene, "^gene_CPM_", ""),
    log10_CPM_1p = log10(CPM + 1)
  ) %>%
  filter(disease != "FTD") %>%
  mutate(
    rna_anno_2ndRound_level_2 = factor(
      rna_anno_2ndRound_level_2,
      levels = selected_cell_types2
    ),
    disease = factor(disease, levels = c("ALS", "Control"))
  ) %>%
  mutate(
    group = paste(rna_anno_2ndRound_level_2, disease),
    group = factor(
      group,
      levels = paste(
        rep(selected_cell_types2, each = 2),
        c("ALS", "Control")
      )
    )
  ) %>%
  arrange(gene)

p <- to_plot_C9orf72 %>%
  ggplot(aes(rna_anno_2ndRound_level_2, log10_CPM_1p))

p +
  geom_violin(
    aes(
      fill = disease,
      group = group
    ),
    draw_quantiles = c(0.5),
    trim = TRUE,
    scale = "width",
    size = 0.2,
    width = 0.8
  ) +
  facet_wrap(
    ~region,
    ncol = 1,
    scales = "fixed",
    strip.position = "left"
  ) +
  scale_fill_manual(
    name = "diagnosis",
    values = disease_color_panel,
    limits = force
  ) +
  xlab("Major cell types") +
  ylab("snRNA gene expression (log10(CPM+1))") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none"
  )

ggsave(
  filename = "./plots/DE_violin_C9orf72_rmFTD.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 2.5,
  useDingbats = FALSE
)


## log sessionInfo
sessionInfo()