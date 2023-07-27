# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)
library(ggrastr)
library(shades)
library(ggbeeswarm)
library(ggpubr)
library(scico)

## read the sctransform-normalized seurat object
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBender_corrected_noQCFilters_scTransformed_clustered_SeuratV4_object.rds")

meta_data_full <- data_obj_sub@meta.data %>%
       rownames_to_column("cell_id") %>%
       as_tibble()

## read meta data from the filtered cells (keep "NK_cell", "Exc_unknown", "Ambiguous")
meta_data <- read_tsv("./metadata_all_cells_2nd_round_annotations.txt")

annotations <- meta_data %>%
       select(cell_id, rna_anno_2ndRound_level_1, rna_anno_2ndRound_level_2, rna_anno_2ndRound_level_3)

res <- meta_data_full %>%
       left_join(annotations, by = "cell_id") %>%
       mutate(
              rna_anno_2ndRound_level_1 = replace_na(rna_anno_2ndRound_level_1, "filtered"),
              rna_anno_2ndRound_level_2 = replace_na(rna_anno_2ndRound_level_2, "filtered"),
              rna_anno_2ndRound_level_3 = replace_na(rna_anno_2ndRound_level_3, "filtered")
       )

data_obj_sub$rna_anno_2ndRound_level_1 <- res$rna_anno_2ndRound_level_1
data_obj_sub$rna_anno_2ndRound_level_2 <- res$rna_anno_2ndRound_level_2
data_obj_sub$rna_anno_2ndRound_level_3 <- res$rna_anno_2ndRound_level_3

# annotate clusters with previous labels
res_annotated <- res %>%
       mutate(
              rna_anno_level1_with_filtered_cells = case_when(
                     seurat_clusters %in% c(2, 14, 33, 7, 28, 9, 15, 17, 27, 34, 8, 20) ~ "Exc_neuron",
                     seurat_clusters %in% c(25, 18, 19, 13, 23, 29) ~ "Inh_neuron",
                     TRUE ~ "Non_neuron"
              )
       )

data_obj_sub$rna_anno_level1_with_filtered_cells <- res_annotated$rna_anno_level1_with_filtered_cells

# load the gene x cell CPM matrix
cpm_mat <- readRDS("./snRNA_noQCfilters_geneByCell_dgCMatrix_RNA_CPM.rds")

# check order
all.equal(colnames(cpm_mat), res_annotated$cell_id)

disease_color_panel <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")

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

# use level 1 markers to check neuronal gene expression
# level_1_markers <- c("SNAP25", "MEF2C", "SLC17A7", "SATB2", "GAD1", "GAD2", "SOX10")
level_1_markers <- c("SNAP25", "MEF2C", "RBFOX3")
# level_1_markers <- c("SNAP25", "MEF2C")

cpm_mat_used <- cpm_mat[level_1_markers, , drop = FALSE]
all.equal(colnames(cpm_mat_used), res_annotated$cell_id)

df_genes_cpm <- cpm_mat_used %>%
       as.matrix() %>%
       t() %>%
       as_tibble() %>%
       rename_all(.funs = function(x) (paste0("gene_CPM_", x))) %>%
       bind_cols(res_annotated)

to_plot_df <- df_genes_cpm %>%
       select(
              starts_with("gene_CPM_"),
              cell_id,
              disease,
              region,
              subject,
              seurat_clusters,
              rna_anno_2ndRound_level_1,
              rna_anno_2ndRound_level_2,
              rna_anno_2ndRound_level_3,
              rna_anno_level1_with_filtered_cells
       ) %>%
       pivot_longer(
              cols = starts_with("gene_CPM_"),
              names_to = "gene",
              values_to = "CPM"
       ) %>%
       mutate(
              gene = str_replace(gene, "^gene_CPM_", ""),
              gene = factor(gene, levels = level_1_markers),
              log10_CPM_1p = log10(CPM + 1),
              is_filtered = if_else(rna_anno_2ndRound_level_1 == "filtered", "filtered", "kept")
       ) %>%
       arrange(gene)

p <- to_plot_df %>% ggplot(aes(disease, log10_CPM_1p))

p + geom_split_violin(
       aes(
              fill = is_filtered
       ),
       draw_quantiles = c(0.5),
       trim = TRUE,
       # scale = "width",
       size = 0.2,
       width = 0.8
) +
       facet_wrap(
              ~ gene + region,
              ncol = 2,
              scales = "fixed"
              # strip.position = strip_position
       ) +
       # ggtitle(selected_region) +
       # scale_fill_manual(
       #        name = "diagnosis",
       #        values = disease_color_panel,
       #        limits = force
       # ) +
       # xlab("Major cell types") +
       # ylab("snRNA gene expression (log10(CPM+1))") +
       theme_bw(base_size = 8, base_family = "Helvetica") +
       theme(
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              strip.background = element_blank(),
              strip.placement = "outside",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
       )


p <- to_plot_df %>% ggplot(aes(paste0(disease, "_", subject), log10_CPM_1p))

p + geom_split_violin(
       aes(
              fill = is_filtered
       ),
       draw_quantiles = c(0.5),
       trim = TRUE,
       scale = "width",
       size = 0.2,
       width = 0.8
) +
       facet_wrap(
              ~ gene + region,
              ncol = 2,
              scales = "fixed"
              # strip.position = strip_position
       ) +
       # ggtitle(selected_region) +
       # scale_fill_manual(
       #        name = "diagnosis",
       #        values = disease_color_panel,
       #        limits = force
       # ) +
       # xlab("Major cell types") +
       # ylab("snRNA gene expression (log10(CPM+1))") +
       theme_bw(base_size = 8, base_family = "Helvetica") +
       theme(
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              strip.background = element_blank(),
              strip.placement = "outside",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
       )



to_plot_df <- df_genes_cpm %>%
       select(
              starts_with("gene_CPM_"),
              cell_id,
              disease,
              region,
              subject,
              seurat_clusters,
              percent_mt,
              rna_anno_2ndRound_level_1,
              rna_anno_2ndRound_level_2,
              rna_anno_2ndRound_level_3,
              rna_anno_level1_with_filtered_cells
       ) %>%
       pivot_longer(
              cols = starts_with("gene_CPM_"),
              names_to = "gene",
              values_to = "CPM"
       ) %>%
       mutate(
              gene = str_replace(gene, "^gene_CPM_", ""),
              gene = factor(gene, levels = level_1_markers),
              log10_CPM_1p = log10(CPM + 1),
              is_filtered = if_else(rna_anno_2ndRound_level_1 == "filtered", "filtered", "kept")
       ) %>%
       arrange(gene)

p <- to_plot_df %>% ggplot(aes(log10_CPM_1p, percent_mt))

p +
       # geom_point_rast() +
       geom_hex(bins = 100) +
       facet_wrap(~ disease + region + gene + is_filtered) +
       scale_fill_viridis_c(trans = "log10") +
       theme_bw()

df_count <- df_genes_cpm %>%
       mutate(
              is_neuron_1 = if_else(gene_CPM_SNAP25 > 0, "neuron", "non-neuron"),
              is_neuron_2 = if_else(gene_CPM_MEF2C > 0, "neuron", "non-neuron"),
              is_neuron_3 = if_else(gene_CPM_RBFOX3 > 0, "neuron", "non-neuron"),
              is_neuron = if_else(gene_CPM_MEF2C > 0 | gene_CPM_SNAP25 > 0 | gene_CPM_RBFOX3 > 0, "neuron", "non-neuron"),
              is_neuron_cluster = if_else(rna_anno_level1_with_filtered_cells == "Non_neuron", "non-neuron", "neuron"),
              is_filtered = if_else(rna_anno_2ndRound_level_1 == "filtered", "filtered", "kept")
       ) %>%
       # count(region, disease, subject, is_neuron_cluster, is_filtered)
       count(region, disease, subject, is_neuron_3, is_filtered)
# count(region, disease, subject, is_neuron, is_filtered)

# p <- df_count %>% ggplot(aes(is_neuron, n))
# p <- df_count %>% ggplot(aes(is_neuron_cluster, n))
p <- df_count %>% ggplot(aes(is_neuron_3, n))
p + geom_bar(aes(fill = is_filtered), position = "fill", stat = "identity") +
       facet_wrap(~ region + disease + subject) +
       theme_bw()

p <- df_count %>% ggplot(aes(is_filtered, n))
# p + geom_bar(aes(fill = is_neuron), position = "fill", stat = "identity") +
# p + geom_bar(aes(fill = is_neuron_cluster), position = "fill", stat = "identity") +
p + geom_bar(aes(fill = is_neuron_3), position = "fill", stat = "identity") +
       facet_wrap(~ region + disease + subject) +
       theme_bw()

# plot the proportion of neurons in FANS-sorted samples
df_FANS <- read_tsv("./cell_proportion_FANS.txt")

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")
disease_color_dark <- lightness(disease_color, scalefac(0.8)) %>%
       as.vector()
names(disease_color_dark) <- names(disease_color)

p <- df_FANS %>%
       mutate(
              disease = factor(disease, levels = c("ALS", "FTD", "Control"))
       ) %>%
       ggplot(aes(disease, percent_neuron))

p +
       geom_boxplot(aes(color = disease), width = 0.5, outlier.shape = NA) +
       geom_beeswarm(aes(color = disease), size = 1, cex = 4, shape = 1) +
       stat_compare_means(
              # aes(label = ..p.signif..),
              method = "t.test",
              comparisons = list(
                     c("ALS", "Control"),
                     c("FTD", "Control"),
                     c("ALS", "FTD")
              )
       ) +
       xlab("Diagnosis") +
       ylab("Proportion of neurons (%)") +
       ggtitle("FANS-sorted proportion") +
       facet_wrap(~region, nrow = 1, strip.position = "bottom") +
       coord_cartesian(ylim = c(0, NA)) +
       scale_color_manual(values = disease_color_dark, name = "diagnosis") +
       theme_bw(base_size = 8, base_family = "Helvetica") +
       theme(
              strip.background = element_blank(),
              panel.grid.minor = element_blank(),
              strip.placement = "outside",
              legend.position = "none"
       )

ggsave("./plots/compare_neurons_proportion_between_disease_diagnosis_in_FANS.pdf",
       device = cairo_pdf(), width = 2, height = 2, useDingbats = FALSE
)
dev.off()

# normalize the proportion of neurons by the total number of cells
p <- df_FANS %>%
       mutate(
              disease = factor(disease, levels = c("ALS", "FTD", "Control")),
              sum_proportion = percent_neuron + percent_oligo + percent_glia,
              normalized_percent_neuron = percent_neuron / sum_proportion,
              normalized_percent_oligo = percent_oligo / sum_proportion,
              normalized_percent_glia = percent_glia / sum_proportion
       ) %>%
       ggplot(aes(disease, normalized_percent_neuron))

p +
       geom_boxplot(aes(color = disease), width = 0.5, outlier.shape = NA) +
       geom_beeswarm(aes(color = disease), size = 1, cex = 4, shape = 1) +
       stat_compare_means(
              aes(label = ..p.signif..),
              method = "t.test",
              comparisons = list(
                     c("ALS", "Control"),
                     c("FTD", "Control"),
                     c("ALS", "FTD")
              )
       ) +
       xlab("Diagnosis") +
       ylab("Proportion of neurons (%)") +
       ggtitle("FANS-sorted proportion") +
       facet_wrap(~region, nrow = 1, strip.position = "bottom") +
       coord_cartesian(ylim = c(0, NA)) +
       scale_color_manual(values = disease_color_dark, name = "diagnosis") +
       theme_bw(base_size = 8, base_family = "Helvetica") +
       theme(
              strip.background = element_blank(),
              panel.grid.minor = element_blank(),
              strip.placement = "outside",
              legend.position = "none"
       )
ggsave("./plots/compare_neurons_proportion_between_disease_diagnosis_in_FANS_normalized.pdf",
       device = cairo_pdf(), width = 2, height = 2, useDingbats = FALSE
)
dev.off()

# Try using centered log ratio transformation
clr_res <- df_FANS %>%
       mutate(percent_other = 100 - percent_neuron - percent_oligo - percent_glia) %>%
       select(starts_with("percent_")) %>%
       as.matrix() %>%
       apply(1, clr) %>%
       t() %>%
       as_tibble() %>%
       rename_all(~ str_replace(., "percent_", "clr_"))

df_FANS_clr <- bind_cols(df_FANS, clr_res)

p <- df_FANS_clr %>%
       mutate(
              disease = factor(disease, levels = c("ALS", "FTD", "Control"))
       ) %>%
       ggplot(aes(disease, clr_neuron))

p +
       geom_boxplot(aes(color = disease), width = 0.5, outlier.shape = NA) +
       geom_beeswarm(aes(color = disease), size = 1, cex = 4, shape = 1) +
       stat_compare_means(
              aes(label = ..p.signif..),
              method = "t.test",
              comparisons = list(
                     c("ALS", "Control"),
                     c("FTD", "Control"),
                     c("ALS", "FTD")
              )
       ) +
       xlab("Diagnosis") +
       ylab("CLR transformed proportion of neurons") +
       facet_wrap(~region, nrow = 1) +
       scale_color_manual(values = disease_color_dark, name = "diagnosis") +
       theme_bw(base_size = 8, base_family = "Helvetica") +
       theme(
              strip.background = element_blank(),
              panel.grid.minor = element_blank(),
              strip.placement = "outside",
              legend.position = "none"
       )
ggsave("./plots/compare_neurons_proportion_between_disease_diagnosis_in_FANS_CLR_transformed.pdf",
       device = cairo_pdf(), width = 2, height = 2, useDingbats = FALSE
)
dev.off()

## log sessionInfo
sessionInfo()
