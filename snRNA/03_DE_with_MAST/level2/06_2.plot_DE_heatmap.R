library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(scico)

DE_level_2 <- read_tsv("MAST_res_level2_summary.txt")

FDR_threshold <- 0.05
# FDR_threshold <- 0.2
# FC_threshold <- 1
# FC_threshold <- 1.2
FC_threshold <- 2

selected_cell_types <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

# selected_cell_types <- c(
#   "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
# )

DE_level_2_sub <- DE_level_2 %>%
  filter(
    cond_1 == "ALS",
    # cond_1 == "FTD",
    cond_2 == "Control",
    cell_type %in% selected_cell_types,
    FDR < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    abs(model_log2FC - avg_logFC) < 2
  )

mat <- DE_level_2_sub %>%
  mutate(group = paste(cell_type, region, sep = "__")) %>%
  select(gene, group, model_log2FC) %>%
  replace_na(list(model_log2FC = 0)) %>%
  pivot_wider(
    names_from = "group",
    values_from = "model_log2FC",
    values_fill = 0
  ) %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  as.matrix()

# define cell-type-specific scores using Tau
get_tau <- function(x) {
  n <- length(x)
  if (max(x) == 0) {
    tau <- 0
  } else {
    x_hat <- x / max(x)
    tau <- (n - sum(x_hat)) / (n - 1)
  }
  return(tau)
}

mat_abs <- abs(mat)
gene_tau <- apply(mat_abs, 1, get_tau)

color_fun <- colorRamp2(seq(-2, 2, length.out = 10), scico(10, palette = "vik"))
color_fun_tau <- colorRamp2(
  seq(0, 1, length.out = 10),
  scico(10, palette = "roma", direction = -1)
)
group_order <- paste(rep(selected_cell_types, each = 2),
  rep(c("MCX", "mFCX"), length(selected_cell_types)),
  sep = "__"
)
mat_col_order <- colnames(mat)[order(match(colnames(mat), group_order))]

region_color <- c("MCX" = "#432266", "mFCX" = "#FAA51B")
cell_type_color_df <- read_tsv("../../color_palette_level_2.txt") %>%
  filter(sub_cluster %in% selected_cell_types) %>%
  mutate(sub_cluster = factor(sub_cluster, levels = selected_cell_types))
cell_type_color <- cell_type_color_df$color
names(cell_type_color) <- cell_type_color_df$sub_cluster

cell_class_color <- c(
  "Non_neuron" = "#9C482B",
  "Exc_neuron" = "#96BB45",
  "Inh_neuron" = "#718DC7"
)

cell_count <- DE_level_2_sub %>%
  distinct(cond_1, cond_2, cell_type, region, num_cells_cond_1, num_cells_cond_2) %>%
  mutate(
    group = paste(cell_type, region, sep = "__"),
    group = factor(group, levels = colnames(mat)),
    num_cells_cond_1 = log10(num_cells_cond_1),
    num_cells_cond_2 = log10(num_cells_cond_2)
  ) %>%
  arrange(group)

col_meta <- tibble(
  group = colnames(mat),
  region = str_extract(colnames(mat), pattern = "(?<=__)\\w+"),
  cell_type = str_extract(colnames(mat), pattern = "\\w+(?=__)"),
  cell_class = case_when(
    grepl("Exc_", cell_type) ~ "Exc_neuron",
    grepl("Inh_", cell_type) ~ "Inh_neuron",
    TRUE ~ "Non_neuron"
  ),
  num_cell_control = cell_count$num_cells_cond_2,
  num_cell_disease = cell_count$num_cells_cond_1
)

mat_col_anno <- HeatmapAnnotation(
  region = col_meta$region,
  cell_class = col_meta$cell_class,
  cell_type = col_meta$cell_type,
  log10_num_cell_control = col_meta$num_cell_control,
  log10_num_cell_disease = col_meta$num_cell_disease,
  col = list(
    region = region_color,
    cell_class = cell_class_color,
    cell_type = cell_type_color
  )
)

selected_DE_genes <- read_tsv("./selected_DE_genes.txt", col_names = c("gene"))
selected_DE_genes_sub <- selected_DE_genes %>% sample_frac(0.1)

mat_row_anno <- rowAnnotation(
  tau = gene_tau,
  gene = anno_mark(
    at = match(selected_DE_genes_sub$gene, rownames(mat)),
    labels = selected_DE_genes_sub$gene
    # at = match(selected_DE_genes$gene, rownames(mat)),
    # labels = selected_DE_genes$gene
  ),
  col = list(tau = color_fun_tau)
)

pdf(("test_DE_heatmap10.pdf"), width = 8, height = 8)
ht <- Heatmap(
  matrix = mat,
  col = color_fun,
  name = "log2FC",
  na_col = "grey",
  # rect_gp = gpar(col = "white", lwd = 2),
  # border_gp = gpar(col = "black", lty = 2),
  color_space = "LAB",
  show_row_names = FALSE,
  show_row_dend = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_reorder = TRUE,
  cluster_columns = FALSE,
  column_order = mat_col_order,
  row_km = 7,
  row_km_repeats = 100,
  # row_split = 8,
  # cluster_row_slices = FALSE,
  top_annotation = mat_col_anno,
  right_annotation = mat_row_anno,
  use_raster = TRUE,
  raster_by_magick = TRUE,
  raster_magick_filter = "Lanczos",
  raster_device = "CairoPNG",
  width = unit(4, "inches"),
  height = unit(4, "inches"),
  # heatmap_width = unit(8, "cm"),
  # heatmap_height = unit(8, "cm"),
)
draw(ht)
dev.off()

# region_labels <- c("MCX" = "motor cortex", "mFCX" = "frontal cortex")
# DE_group_color <- c("up" = "#DE315D", "down" = "#14668E")

# show all FC for genes with FC > 2 in at least one group
selected_genes <- DE_level_2_sub %>%
  distinct(gene) %>%
  pull(gene)

DE_level_2_sub_2 <- DE_level_2 %>%
  filter(
    cond_1 == "ALS",
    cond_2 == "Control",
    cell_type %in% selected_cell_types,
    gene %in% selected_genes
  )

mat_2 <- DE_level_2_sub_2 %>%
  mutate(group = paste(cell_type, region, sep = "__")) %>%
  mutate(model_log2FC = if_else(FDR < 0.2, model_log2FC, 0)) %>%
  select(gene, group, model_log2FC) %>%
  replace_na(list(model_log2FC = 0)) %>%
  pivot_wider(
    names_from = "group",
    values_from = "model_log2FC",
    values_fill = 0
  ) %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  as.matrix()

mat_abs_2 <- abs(mat_2)
gene_tau_2 <- apply(mat_abs_2, 1, get_tau)

mat_col_order_2 <- colnames(mat_2)[order(match(colnames(mat_2), group_order))]

col_meta_2 <- tibble(
  group = colnames(mat_2),
  region = str_extract(colnames(mat_2), pattern = "(?<=__)\\w+"),
  cell_type = str_extract(colnames(mat_2), pattern = "\\w+(?=__)"),
  cell_class = case_when(
    grepl("Exc_", cell_type) ~ "Exc_neuron",
    grepl("Inh_", cell_type) ~ "Inh_neuron",
    TRUE ~ "Non_neuron"
  )
)

mat_col_anno_2 <- HeatmapAnnotation(
  region = col_meta_2$region,
  cell_class = col_meta_2$cell_class,
  cell_type = col_meta_2$cell_type,
  col = list(
    region = region_color,
    cell_class = cell_class_color,
    cell_type = cell_type_color
  )
)

mat_row_anno_2 <- rowAnnotation(
  tau = gene_tau_2,
  gene = anno_mark(
    # at = match(selected_DE_genes_sub$gene, rownames(mat)),
    # labels = selected_DE_genes_sub$gene
    at = match(selected_DE_genes_sub$gene, rownames(mat_2)),
    labels = selected_DE_genes_sub$gene
  ),
  col = list(tau = color_fun_tau)
)

pdf(("test_DE_heatmap11.pdf"), width = 8, height = 8)
ht <- Heatmap(
  matrix = mat_2,
  col = color_fun,
  name = "log2FC",
  na_col = "grey",
  # rect_gp = gpar(col = "white", lwd = 2),
  # border_gp = gpar(col = "black", lty = 2),
  color_space = "LAB",
  show_row_names = FALSE,
  show_row_dend = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_reorder = TRUE,
  cluster_columns = FALSE,
  column_order = mat_col_order_2,
  row_km = 12,
  row_km_repeats = 100,
  # row_split = 8,
  # cluster_row_slices = FALSE,
  top_annotation = mat_col_anno_2,
  right_annotation = mat_row_anno_2,
  use_raster = TRUE,
  raster_by_magick = TRUE,
  raster_magick_filter = "Lanczos",
  raster_device = "CairoPNG",
  width = unit(4, "inches"),
  height = unit(4, "inches"),
  # heatmap_width = unit(8, "cm"),
  # heatmap_height = unit(8, "cm"),
)
draw(ht)
dev.off()