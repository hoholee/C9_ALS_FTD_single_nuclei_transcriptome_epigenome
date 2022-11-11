# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(scales)
library(ggrastr)
library(scattermore)
library(cowplot)

# read meta data
meta_data <- read_tsv("./metadata_merged_addAnno.txt")

# shuffle the order of dots to alleviate overplotting
set.seed(seed = 666)
meta_data_shuffle <- meta_data %>%
  sample_frac() %>%
  mutate(sample = Sample) %>%
  separate(sample, c("region", "disease", "subject"), sep = "_") %>%
  mutate(group = paste(atac_anno_level_2, region, disease, sep = "_"))

# manually plot UMAP plots colored by meta data
p <- meta_data_shuffle %>% ggplot(aes(full_UMAP_1, full_UMAP_2))

# brain region
region_labels <- c("MCX" = "motor cortex", "mFCX" = "medial frontal cortex")
region_color <- c("MCX" = "#432266", "mFCX" = "#FAA51B")
relabel_region <- function(x) {
  region_labels[x]
}
p_region <- p +
  geom_scattermore(aes(color = region), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(labels = relabel_region(), values = region_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_region_panel <- p_region + theme(legend.position = "none")
p_region_legend <- get_legend(
  p_region +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(p_region_panel, p_region_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_full_colored_by_region.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

# disease
# disease_color <- c("ALS" = "#431853", "FTD" = "#27868E", "Control" = "#BAD535")
# disease_color <- c("ALS" = "#cc4b67", "FTD" = "#ccbd4b", "Control" = "#4b83cc")
# disease_color <- c("ALS" = "#bf4c8d", "FTD" = "#c0c252", "Control" = "#4cbfae")
disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")
p_disease <- p +
  geom_scattermore(aes(color = disease), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(values = disease_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_disease_panel <- p_disease + theme(legend.position = "none")
p_disease_legend <- get_legend(
  p_disease +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(p_disease_panel, p_disease_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_full_colored_by_disease.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

# sex
# sex_labels <- c("F" = "female", "M" = "male")
# sex_color <- c("F" = "#FC9E89", "M" = "#3D304C")
# relabel_sex <- function(x){sex_labels[x]}
# p_sex <- p +
#   geom_scattermore(aes(color = sex), pointsize = 4, pixels = c(2048, 2048)) +
#   scale_color_manual(labels = relabel_sex(), values = sex_color) +
#   theme_classic(base_size = 8, base_family = "Helvetica")
# p_sex_panel <- p_sex + theme(legend.position = "none")
# p_sex_legend <- get_legend(
#   p_sex +
#     guides(color = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom")
# )

# plot_grid(p_sex_panel, p_sex_legend, ncol = 1, rel_heights = c(1, 0.1))
# ggsave("./plots/umap_full_colored_by_sex.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)


# batch
# batch_labels <- c("batch_1" = "1", "batch_2" = "2", "batch_3" = "3", "batch_4" = "4", "batch_5" = "5", "batch_6" = "6")
# batch_color <- c("batch_1" = "#A8DBA8",
#                  "batch_2" = "#00A0B0",
#                  "batch_3" = "#6A4A3C",
#                  "batch_4" = "#CC333F",
#                  "batch_5" = "#EB6841",
#                  "batch_6" = "#EDC951"
#                  )
# relabel_batch <- function(x){batch_labels[x]}
# p_batch <- p +
#   geom_scattermore(aes(color = seq_batch), pointsize = 4, pixels = c(2048, 2048)) +
#   scale_color_manual(labels = relabel_batch(), values = batch_color) +
#   theme_classic(base_size = 8, base_family = "Helvetica")
# p_batch_panel <- p_batch + theme(legend.position = "none")
# p_batch_legend <- get_legend(
#   p_batch +
#     guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
#     theme(legend.position = "bottom")
# )

# plot_grid(p_batch_panel, p_batch_legend, ncol = 1, rel_heights = c(1, 0.1))
# ggsave("./plots/umap_full_colored_by_batch.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)


# individuals
individual_color <- c(
  "332" = "#92072a",
  "111" = "#c30938",
  "113" = "#f40b45",
  "52" = "#f63c6b",
  "388" = "#f86d90",
  "110" = "#fa9eb5",
  "36" = "#ffd866",
  "54" = "#ffcb33",
  "55" = "#ffbe00",
  "674" = "#cc9800",
  "61" = "#997200",
  "945" = "#00998d",
  "902" = "#00ccbb",
  "91" = "#00ffea",
  "1069" = "#33ffee",
  "904" = "#66fff3",
  "906" = "#99fff7"
)
p_individual <- p +
  geom_scattermore(aes(color = subject), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(values = individual_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_individual_panel <- p_individual + theme(legend.position = "none")
p_individual_legend <- get_legend(
  p_individual +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom")
)

plot_grid(p_individual_panel, p_individual_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_full_colored_by_individual.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

# cell class (level 1 annotation)
cell_class_labels <- c("Glia" = "non-neurons", "Exc" = "excitatory neurons", "Inh" = "inhibitory neurons")
cell_class_color <- c("Glia" = "#9C482B", "Exc" = "#96BB45", "Inh" = "#718DC7")
relabel_cell_class <- function(x) {
  cell_class_labels[x]
}
p_cell_class <- p +
  geom_scattermore(aes(color = atac_anno_level_1), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Cell class", labels = relabel_cell_class(), values = cell_class_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_cell_class_panel <- p_cell_class + theme(legend.position = "none")
p_cell_class_legend <- get_legend(
  p_cell_class +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(p_cell_class_panel, p_cell_class_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_full_colored_by_cell_class.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

# extract sub-clustering meta data
meta_data_exc <- meta_data %>%
  filter(atac_anno_level_1 == "Exc") %>%
  mutate(sample = Sample) %>%
  separate(sample, c("region", "disease", "subject"), sep = "_") %>%
  mutate(group = paste(atac_anno_level_2, region, disease, sep = "_"))
meta_data_inh <- meta_data %>%
  filter(atac_anno_level_1 == "Inh") %>%
  mutate(sample = Sample) %>%
  separate(sample, c("region", "disease", "subject"), sep = "_") %>%
  mutate(group = paste(atac_anno_level_2, region, disease, sep = "_"))
meta_data_glia <- meta_data %>%
  filter(atac_anno_level_1 == "Glia") %>%
  mutate(sample = Sample) %>%
  separate(sample, c("region", "disease", "subject"), sep = "_") %>%
  mutate(group = paste(atac_anno_level_2, region, disease, sep = "_"))

# shuffle the order of dots
set.seed(seed = 666)
meta_data_exc_shuffle <- meta_data_exc %>% sample_frac()
set.seed(seed = 666)
meta_data_inh_shuffle <- meta_data_inh %>% sample_frac()
set.seed(seed = 666)
meta_data_glia_shuffle <- meta_data_glia %>% sample_frac()

# manually plot UMAP plots, colored by subclusters (level 3 annotations)
# Exc
p <- meta_data_exc_shuffle %>% ggplot(aes(Exc_UMAP_1, Exc_UMAP_2))

exc_palette <- read_tsv("color_palette_exc_subclusters.txt")
exc_subcluster_color <- exc_palette$color
names(exc_subcluster_color) <- exc_palette$sub_cluster
relabel_exc_subcluster <- function(x) {
  str_replace_all(x, "_", " ")
}

p_exc_subcluster <- p +
  geom_scattermore(aes(color = atac_anno_level_3), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Exc subcluster", labels = relabel_exc_subcluster, values = exc_subcluster_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_exc_subcluster_panel <- p_exc_subcluster + theme(legend.position = "none")
p_exc_subcluster_legend <- get_legend(
  p_exc_subcluster +
    guides(color = guide_legend(ncol = 4)) +
    theme(legend.position = "bottom")
)

plot_grid(p_exc_subcluster_panel, p_exc_subcluster_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_exc_subclustering_colored_by_subtype_level3.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

exc_palette <- read_tsv("color_palette_exc_subclusters_level2.txt")
exc_subcluster_color <- exc_palette$color
names(exc_subcluster_color) <- exc_palette$sub_cluster
relabel_exc_subcluster <- function(x) {
  str_replace_all(x, "_", " ")
}

p_exc_subcluster <- p +
  geom_scattermore(aes(color = atac_anno_level_2), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Exc subcluster", labels = relabel_exc_subcluster, values = exc_subcluster_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_exc_subcluster_panel <- p_exc_subcluster + theme(legend.position = "none")
p_exc_subcluster_legend <- get_legend(
  p_exc_subcluster +
    guides(color = guide_legend(ncol = 4)) +
    theme(legend.position = "bottom")
)

plot_grid(p_exc_subcluster_panel, p_exc_subcluster_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_exc_subclustering_colored_by_subtype_level2.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

# Inh
p <- meta_data_inh_shuffle %>% ggplot(aes(Inh_UMAP_1, Inh_UMAP_2))

inh_palette <- read_tsv("color_palette_inh_subclusters.txt")
inh_subcluster_color <- inh_palette$color
names(inh_subcluster_color) <- inh_palette$sub_cluster
relabel_inh_subcluster <- function(x) {
  str_replace_all(x, "_", " ")
}

p_inh_subcluster <- p +
  geom_scattermore(aes(color = atac_anno_level_3), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Inh subcluster", labels = relabel_inh_subcluster, values = inh_subcluster_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_inh_subcluster_panel <- p_inh_subcluster + theme(legend.position = "none")
p_inh_subcluster_legend <- get_legend(
  p_inh_subcluster +
    guides(color = guide_legend(ncol = 4)) +
    theme(legend.position = "bottom")
)

plot_grid(p_inh_subcluster_panel, p_inh_subcluster_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_inh_subclustering_colored_by_subtype_level3.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

inh_palette <- read_tsv("color_palette_inh_subclusters_level2.txt")
inh_subcluster_color <- inh_palette$color
names(inh_subcluster_color) <- inh_palette$sub_cluster
relabel_inh_subcluster <- function(x) {
  str_replace_all(x, "_", " ")
}

p_inh_subcluster <- p +
  geom_scattermore(aes(color = atac_anno_level_2), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Inh subcluster", labels = relabel_inh_subcluster, values = inh_subcluster_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_inh_subcluster_panel <- p_inh_subcluster + theme(legend.position = "none")
p_inh_subcluster_legend <- get_legend(
  p_inh_subcluster +
    guides(color = guide_legend(ncol = 4)) +
    theme(legend.position = "bottom")
)

plot_grid(p_inh_subcluster_panel, p_inh_subcluster_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_inh_subclustering_colored_by_subtype_level2.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

# Glia
p <- meta_data_glia_shuffle %>% ggplot(aes(Glia_UMAP_1, Glia_UMAP_2))

glia_palette <- read_tsv("color_palette_glia_subclusters.txt")
glia_subcluster_color <- glia_palette$color
names(glia_subcluster_color) <- glia_palette$sub_cluster
relabel_glia_subcluster <- function(x) {
  str_replace_all(x, "_", " ")
}

p_glia_subcluster <- p +
  geom_scattermore(aes(color = atac_anno_level_3), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Non-neuron subcluster", labels = relabel_glia_subcluster, values = glia_subcluster_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_glia_subcluster_panel <- p_glia_subcluster + theme(legend.position = "none")
p_glia_subcluster_legend <- get_legend(
  p_glia_subcluster +
    guides(color = guide_legend(ncol = 4)) +
    theme(legend.position = "bottom")
)

plot_grid(p_glia_subcluster_panel, p_glia_subcluster_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_glia_subclustering_colored_by_subtype_level3.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

glia_palette <- read_tsv("color_palette_glia_subclusters_level2.txt")
glia_subcluster_color <- glia_palette$color
names(glia_subcluster_color) <- glia_palette$sub_cluster
relabel_glia_subcluster <- function(x) {
  str_replace_all(x, "_", " ")
}

p_glia_subcluster <- p +
  geom_scattermore(aes(color = atac_anno_level_2), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Non-neuron subcluster", labels = relabel_glia_subcluster, values = glia_subcluster_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_glia_subcluster_panel <- p_glia_subcluster + theme(legend.position = "none")
p_glia_subcluster_legend <- get_legend(
  p_glia_subcluster +
    guides(color = guide_legend(ncol = 4)) +
    theme(legend.position = "bottom")
)

plot_grid(p_glia_subcluster_panel, p_glia_subcluster_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_glia_subclustering_colored_by_subtype_level2.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

## log sessionInfo
sessionInfo()