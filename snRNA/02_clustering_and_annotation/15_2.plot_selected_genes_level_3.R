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
DE_level3 <- read_tsv("./DE_MAST/second_round_level_3/MAST_res_level3_summary.txt")

# todo: clean up and define a plotting function
selected_region <- "MCX"
# selected_region <- "mFCX"
selected_disease_1 <- "ALS"
selected_disease_2 <- "Control"

selected_cell_types <- c(
  "Astro_CD44",
  "Astro_HPSE2"
)

selected_genes <- c(
  "RANBP3L", "CD44", "GFAP", "CHI3L1"
)

FDR_threshold <- 0.05
FC_threshold <- 1.2

selected_diseases <- c(selected_disease_1, selected_disease_2)
cell_idx <- which(meta_data$region == selected_region &
  meta_data$disease %in% selected_diseases &
  meta_data$rna_anno_2ndRound_level_3 %in% selected_cell_types)

meta_data_used <- meta_data %>%
  filter(
    region == selected_region,
    disease %in% selected_diseases,
    rna_anno_2ndRound_level_3 %in% selected_cell_types
  )

cpm_mat_used <- cpm_mat[selected_genes, cell_idx, drop = FALSE]
all.equal(colnames(cpm_mat_used), meta_data_used$cell_id)

DE_used <- DE_level3 %>%
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
      sig == "DE" & FDR < 0.05 ~ "*",
      sig == "DE" & FDR < 0.01 ~ "**",
      sig == "DE" & FDR < 0.001 ~ "***",
      sig == "DE" & FDR < 0.0001 ~ "****",
      TRUE ~ "n.s."
    ),
    DE_dir = if_else(model_log2FC > 0, "up-regulated", "down-regulated")
  )

df_genes_cpm <- cpm_mat_used %>%
  as.matrix() %>%
  t() %>%
  as_tibble() %>%
  rename_all(.funs = function(x) (paste0("gene_CPM_", x))) %>%
  bind_cols(meta_data_used)

to_plot_df <- df_genes_cpm %>%
  select(starts_with("gene_CPM_"), disease, subject, rna_anno_2ndRound_level_2, rna_anno_2ndRound_level_3) %>%
  pivot_longer(cols = starts_with("gene_CPM_"), names_to = "gene", values_to = "CPM") %>%
  mutate(
    gene = str_replace(gene, "^gene_CPM_", ""),
    gene = factor(gene, levels = selected_genes),
    rna_anno_2ndRound_level_3 = factor(rna_anno_2ndRound_level_3, levels = selected_cell_types),
    log10_CPM_1p = log10(CPM + 1)
  ) %>%
  left_join(DE_used, by = c("gene" = "gene", "rna_anno_2ndRound_level_3" = "cell_type"))

p <- to_plot_df %>% ggplot(aes(rna_anno_2ndRound_level_3, log10_CPM_1p))
p + geom_violin(aes(
  fill = disease,
  # color = disease,
  group = paste(rna_anno_2ndRound_level_3, disease)
),
draw_quantiles = c(0.5),
trim = TRUE,
scale = "width",
size = 0.2,
width = 0.8,
position = position_dodge(width = 0.8)
) +
  geom_text(aes(x = cell_type, y = 3.2, label = sig_sign), data = DE_used, size = 2) +
  facet_wrap(~ DE_dir + gene, ncol = 2, scales = "free_y", strip.position = "left") +
  ggtitle(selected_region) +
  scale_fill_manual(values = disease_color_panel) +
  scale_color_manual(values = lightness(disease_color_panel, scalefac(0.5))) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# ggsave("./plots/DE_ALS_vs_Control_mFCX_Astro.pdf", device = cairo_pdf(), width = 4, height = 3, useDingbats = F)

# split violin plot
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

p + geom_split_violin(aes(
  fill = disease,
  # color = disease,
  group = paste(rna_anno_2ndRound_level_3, disease)
),
draw_quantiles = c(0.5),
trim = TRUE,
scale = "width",
size = 0.2,
width = 0.8
) +
  geom_text(aes(x = cell_type, y = 4, label = sig_sign), data = DE_used, size = 2) +
  facet_wrap(~ DE_dir + gene, nrow = 1, scales = "fixed", strip.position = "top") +
  ggtitle(selected_region) +
  scale_fill_manual(values = disease_color_panel) +
  scale_color_manual(values = lightness(disease_color_panel, scalefac(0.5))) +
  xlab("Level 3 subtypes") +
  ylab("log10(CPM+1)") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave("./plots/DE_ALS_vs_Control_MCX_example_genes_Astro_split_level_3.pdf",
  device = cairo_pdf(), width = 4.2, height = 2, useDingbats = F
)

## log sessionInfo
sessionInfo()