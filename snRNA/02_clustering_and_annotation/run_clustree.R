# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(scales)
library(ggrastr)
library(scattermore)
library(clustree)
library(ggpubr)

# read metadata
meta_data <- read_tsv("./metadata_all_cells_2nd_round_annotations_clean.txt")
meta_Exc <- read_tsv("./meta_annotated_Exc_sub_clusters.txt") %>%
  select(cell_id, seurat_sub_clusters)
meta_Inh <- read_tsv("./meta_annotated_Inh_sub_clusters.txt") %>%
  select(cell_id, seurat_sub_clusters)
meta_NonNeu <- read_tsv("./meta_annotated_NonNeu_sub_clusters.txt") %>%
  select(cell_id, seurat_sub_clusters)

# read cluster assignments from multiple resolutions (excluding the selected resolutions in the manuscript)
df_full <- read_tsv("./snRNA_fullDataset_cellBender_corrected_leiden_clusters_res_merged.txt") %>%
  pivot_wider(
    names_from = leiden_res,
    values_from = seurat_clusters,
    names_prefix = "cluster_res_"
  )

df_Exc <- read_tsv("./snRNA_subclustering_Exc_cellBender_corrected_leiden_clusters_res_merged.txt") %>%
  pivot_wider(
    names_from = leiden_res,
    values_from = seurat_Exc_subclusters,
    names_prefix = "cluster_res_"
  )

df_Inh <- read_tsv("./snRNA_subclustering_Inh_cellBender_corrected_leiden_clusters_res_merged.txt") %>%
  pivot_wider(
    names_from = leiden_res,
    values_from = seurat_Inh_subclusters,
    names_prefix = "cluster_res_"
  )

df_NonNeu <- read_tsv("./snRNA_subclustering_NonNeu_cellBender_corrected_leiden_clusters_res_merged.txt") %>%
  pivot_wider(
    names_from = leiden_res,
    values_from = seurat_NonNeu_subclusters,
    names_prefix = "cluster_res_"
  )

# Add back selected resolutions in the manuscript
# Full dataset: r = 1; Exc: r = 1.5; Inh: r = 0.5; NonNeu: r = 1
# Also add a cluster r = 0 to show the root
# Remove some resolutions to keep the clustree simple
# Now keeping r values of 0, 0.05, 0.1, 0.5, 1, 1.5, 2, 3
df_full_merged <- meta_data %>%
  left_join(df_full, by = "cell_id") %>%
  mutate(
    cluster_res_1 = SCT_snn_res.1,
    cluster_res_0 = "1",
  ) %>%
  select(
    -cluster_res_0.3,
    -cluster_res_0.8,
    -cluster_res_1.8,
    -cluster_res_2.5,
    -cluster_res_3.5
  )

df_Exc_merged <- meta_data %>%
  filter(rna_anno_2ndRound_level_1 == "Exc_neuron") %>%
  left_join(meta_Exc, by = "cell_id") %>%
  left_join(df_Exc, by = "cell_id") %>%
  mutate(
    cluster_res_1.5 = seurat_sub_clusters,
    cluster_res_0 = "1",
  ) %>%
  select(
    -cluster_res_0.3,
    -cluster_res_0.8,
    -cluster_res_1.8,
    -cluster_res_2.5,
    -cluster_res_3.5
  )

df_Inh_merged <- meta_data %>%
  filter(rna_anno_2ndRound_level_1 == "Inh_neuron") %>%
  left_join(meta_Inh, by = "cell_id") %>%
  left_join(df_Inh, by = "cell_id") %>%
  mutate(
    cluster_res_0.5 = seurat_sub_clusters,
    cluster_res_0 = "1",
  ) %>%
  select(
    -cluster_res_0.3,
    -cluster_res_0.8,
    -cluster_res_1.8,
    -cluster_res_2.5,
    -cluster_res_3.5
  )

df_NonNeu_merged <- meta_data %>%
  filter(rna_anno_2ndRound_level_1 == "Non_neuron") %>%
  left_join(meta_NonNeu, by = "cell_id") %>%
  left_join(df_NonNeu, by = "cell_id") %>%
  mutate(
    cluster_res_1 = seurat_sub_clusters,
    cluster_res_0 = "1"
  ) %>%
  select(
    -cluster_res_0.3,
    -cluster_res_0.8,
    -cluster_res_1.8,
    -cluster_res_2.5,
    -cluster_res_3.5
  )

# define functions to find the label of the major cell type in each cluster

# label_position <- function(labels) {
#   if (length(unique(labels)) == 1) {
#     position <- as.character(unique(labels))
#   } else {
#     position <- "mixed"
#   }
#   return(position)
# }

label_count_percentage <- function(labels) {
  label_count <- table(labels)
  max_count <- max(label_count)
  max_label <- names(label_count[label_count == max_count])
  label_percentage <- signif(max_count / length(labels) * 100, digits = 2)
  label_output <- str_glue("{max_label} ({label_percentage}%)")
  return(label_output)
}

# plot clustree
selected_node_size_range <- c(8, 20)
selected_node_label_size <- 1
selected_SC3_stability_limits <- c(0, 1)
selected_SC3_stability_breaks <- c(0, 0.5, 1)
selected_in_prop_limits <- c(0, 1)
selected_in_prop_breaks <- c(0.25, 0.5, 0.75, 1)
selected_cluster_size_limits <- c(1, 110000)
selected_cluster_size_breaks <- c(10, 1000, 100000)
selected_num_cell_limits <- c(1, 110000)
selected_num_cell_breaks <- c(1, 10, 100, 1000, 10000, 100000)
selected_edge_width <- 0.5 

p1 <- clustree(
  df_full_merged,
  prefix = "cluster_res_",
  layout = "sugiyama",
  node_size_range = selected_node_size_range,
  node_colour = "sc3_stability",
  edge_width = selected_edge_width,
  # node_label = "rna_anno_2ndRound_level_3",
  # node_label_aggr = "label_count_percentage",
  # node_label_size = selected_node_label_size
) +
  scale_color_viridis_c(
    name = "SC3 stability",
    limits = selected_SC3_stability_limits,
    breaks = selected_SC3_stability_breaks,
    option = "plasma"
  ) +
  scale_fill_viridis_c(
    name = "SC3 stability",
    limits = selected_SC3_stability_limits,
    breaks = selected_SC3_stability_breaks,
    option = "plasma"
  ) +
  scale_edge_alpha(
    name = "In-proportion",
    limits = selected_in_prop_limits,
    breaks = selected_in_prop_breaks
  ) +
  scale_size(
    name = "Cluster size",
    limits = selected_cluster_size_limits,
    breaks = selected_cluster_size_breaks,
    # trans = "log10"
  ) +
  scale_edge_color_gradientn(
    colours = viridis::viridis(256, option = "viridis"),
    name = "# Nuclei",
    trans = "log10",
    limits = selected_num_cell_limits,
    breaks = selected_num_cell_breaks,
    oob = squish
  ) +
  theme(legend.position = "bottom")

p2 <- clustree(
  df_Exc_merged,
  prefix = "cluster_res_",
  layout = "sugiyama",
  node_size_range = selected_node_size_range,
  node_colour = "sc3_stability",
  edge_width = selected_edge_width,
  # node_label = "rna_anno_2ndRound_level_3",
  # node_label_aggr = "label_count_percentage",
  # node_label_size = selected_node_label_size
) +
  scale_color_viridis_c(
    name = "SC3 stability",
    limits = selected_SC3_stability_limits,
    breaks = selected_SC3_stability_breaks,
    option = "plasma"
  ) +
  scale_fill_viridis_c(
    name = "SC3 stability",
    limits = selected_SC3_stability_limits,
    breaks = selected_SC3_stability_breaks,
    option = "plasma"
  ) +
  scale_edge_alpha(
    name = "In-proportion",
    limits = selected_in_prop_limits,
    breaks = selected_in_prop_breaks
  ) +
  scale_size(
    name = "Cluster size",
    limits = selected_cluster_size_limits,
    breaks = selected_cluster_size_breaks,
    # trans = "log10"
  ) +
  scale_edge_color_gradientn(
    colours = viridis::viridis(256, option = "viridis"),
    name = "# Nuclei",
    trans = "log10",
    limits = selected_num_cell_limits,
    breaks = selected_num_cell_breaks,
    oob = squish
  ) +
  theme(legend.position = "bottom")

p3 <- clustree(
  df_Inh_merged,
  prefix = "cluster_res_",
  layout = "sugiyama",
  node_size_range = selected_node_size_range,
  node_colour = "sc3_stability",
  edge_width = selected_edge_width,
  # node_label = "rna_anno_2ndRound_level_3",
  # node_label_aggr = "label_count_percentage",
  # node_label_size = selected_node_label_size
) +
  scale_color_viridis_c(
    name = "SC3 stability",
    limits = selected_SC3_stability_limits,
    breaks = selected_SC3_stability_breaks,
    option = "plasma"
  ) +
  scale_fill_viridis_c(
    name = "SC3 stability",
    limits = selected_SC3_stability_limits,
    breaks = selected_SC3_stability_breaks,
    option = "plasma"
  ) +
  scale_edge_alpha(
    name = "In-proportion",
    limits = selected_in_prop_limits,
    breaks = selected_in_prop_breaks
  ) +
  scale_size(
    name = "Cluster size",
    limits = selected_cluster_size_limits,
    breaks = selected_cluster_size_breaks,
    # trans = "log10"
  ) +
  scale_edge_color_gradientn(
    colours = viridis::viridis(256, option = "viridis"),
    name = "# Nuclei",
    trans = "log10",
    limits = selected_num_cell_limits,
    breaks = selected_num_cell_breaks,
    oob = squish
  ) +
  theme(legend.position = "bottom")

p4 <- clustree(
  df_NonNeu_merged,
  prefix = "cluster_res_",
  layout = "sugiyama",
  node_size_range = selected_node_size_range,
  node_colour = "sc3_stability",
  edge_width = selected_edge_width,
  # node_label = "rna_anno_2ndRound_level_3",
  # node_label_aggr = "label_count_percentage",
  # node_label_size = selected_node_label_size
) +
  scale_color_viridis_c(
    name = "SC3 stability",
    limits = selected_SC3_stability_limits,
    breaks = selected_SC3_stability_breaks,
    option = "plasma"
  ) +
  scale_fill_viridis_c(
    name = "SC3 stability",
    limits = selected_SC3_stability_limits,
    breaks = selected_SC3_stability_breaks,
    option = "plasma"
  ) +
  scale_edge_alpha(
    name = "In-proportion",
    limits = selected_in_prop_limits,
    breaks = selected_in_prop_breaks
  ) +
  scale_size(
    name = "Cluster size",
    limits = selected_cluster_size_limits,
    breaks = selected_cluster_size_breaks,
    # trans = "log10"
  ) +
  scale_edge_color_gradientn(
    colours = viridis::viridis(256, option = "viridis"),
    name = "# Nuclei",
    trans = "log10",
    limits = selected_num_cell_limits,
    breaks = selected_num_cell_breaks,
    oob = squish
  ) +
  theme(legend.position = "bottom")

ggarrange(p1, p2, p3, p4, ncol = 1, common.legend = TRUE, legend = "bottom")

ggsave(
  "./plots/clustree.pdf",
  device = cairo_pdf(),
  width = 10,
  height = 20,
  useDingbats = F
)

## log sessionInfo
sessionInfo()
