# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(pheatmap)
library(viridis)
library(scico)
library(pdfCluster)
library(ggrastr)

# data_obj <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_clean_addAIBSprediction_SeuratV4_object.rds")
data_obj <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_cleanV2_addAIBSprediction_SeuratV4_object.rds")

meta_data <- data_obj@meta.data %>% as_tibble()

# compute adjusted rand index
adj.rand.index(meta_data$rna_anno_2ndRound_level_3, meta_data$predicted.id)

# plot cumulative distribution of maximum prediction score by cell type
color_palette_level_2 <- read_tsv("./color_palette_level_2.txt")
colors_level_2 <- color_palette_level_2$color
names(colors_level_2) <- color_palette_level_2$sub_cluster

p <- meta_data %>%
  mutate(
    region = factor(region, levels = c("MCX", "mFCX")),
    disease = factor(disease,
      levels = c("ALS", "FTD", "Control")
    ),
    rna_anno_2ndRound_level_2 = factor(rna_anno_2ndRound_level_2,
      levels = color_palette_level_2$sub_cluster
    )
  ) %>%
  ggplot(aes(prediction.score.max))

p +
  rasterize(stat_ecdf(aes(color = rna_anno_2ndRound_level_2), pad = FALSE), dpi = 600) +
  geom_vline(
    xintercept = 0.5,
    color = "darkgrey",
    linetype = 2
  ) +
  facet_grid(
    region ~ disease,
    scales = "fixed"
  ) +
  scale_color_manual(
    values = colors_level_2,
    name = "Cell type"
  ) +
  xlab("Maximum prediction score from label transfer") +
  ylab("Cumulative proportion") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "bottom",
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 6)
  )
ggsave(
  "./plots/cumulative_distribution_max_prediction_scores_in_labelTransfer_AIBS_by_level_2_cell_types.pdf",
  device = cairo_pdf(),
  width = 4, height = 4, useDingbats = FALSE
)

# plot confusion matrix
cluster_n_cells <- meta_data %>%
  count(rna_anno_2ndRound_level_3) %>%
  rename(total_n = n)
match_count_df <- meta_data %>%
  count(rna_anno_2ndRound_level_3, predicted.id) %>%
  left_join(cluster_n_cells) %>%
  mutate(n_percent = n / total_n)

uniq_row_name <- match_count_df$rna_anno_2ndRound_level_3 %>% unique()
row_order_id <- c(
  str_which(uniq_row_name, "Exc"),
  str_which(uniq_row_name, "Inh"),
  str_which(uniq_row_name, "Astro"),
  str_which(uniq_row_name, "Micro"),
  str_which(uniq_row_name, "Oligo"),
  str_which(uniq_row_name, "OPC"),
  str_which(uniq_row_name, "Endo"),
  str_which(uniq_row_name, "VLMC")
)
row_order <- uniq_row_name[row_order_id]

match_rank <- match_count_df %>%
  select(rna_anno_2ndRound_level_3, predicted.id, n_percent) %>%
  group_by(predicted.id) %>%
  top_n(1) %>%
  ungroup() %>%
  arrange(match(rna_anno_2ndRound_level_3, row_order))

match_count_ht <- match_count_df %>%
  select(rna_anno_2ndRound_level_3, predicted.id, n_percent) %>%
  arrange(match(rna_anno_2ndRound_level_3, row_order)) %>%
  pivot_wider(names_from = predicted.id, values_from = n_percent, values_fill = 0) %>%
  select(rna_anno_2ndRound_level_3, one_of(match_rank$predicted.id)) %>%
  as.data.frame() %>%
  column_to_rownames("rna_anno_2ndRound_level_3")

pheatmap(match_count_ht,
  cluster_rows = F,
  cluster_cols = F,
  color = colorRampPalette(scico(100, palette = "acton", direction = -1))(100),
  # color = colorRampPalette(viridis(100, option = "A", end = 0.97))(100),
  border_color = NA,
  cellwidth = 5.5,
  cellheight = 5.5,
  fontsize = 5.5,
  # filename = "./plots/confusion_matrix_cellBender_corrected_rna_anno_level_3_vs_AIBS_annotation.pdf",
  filename = "./plots/confusion_matrix_cellBender_corrected_cleanV2_rna_anno_level_3_vs_AIBS_annotation.pdf",
  width = 12, height = 6
)