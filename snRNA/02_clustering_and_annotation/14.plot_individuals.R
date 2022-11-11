# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(scales)
library(scico)
library(tictoc)
library(shades)
library(ComplexHeatmap)
library(circlize)
library(viridis)

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

# todo: clean up and define a plotting function
selected_region <- "MCX"
# selected_region <- "mFCX"
selected_disease_1 <- "ALS"
selected_disease_2 <- "Control"
selected_cell_types <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)
selected_genes <- c("RANBP3L", "CD44", "GFAP", "CHI3L1")
# selected_genes <- c("HSPA8")

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
    rna_anno_2ndRound_level_2 = factor(rna_anno_2ndRound_level_2, levels = selected_cell_types),
    log10_CPM_1p = log10(CPM + 1)
  )

p <- to_plot_df %>% ggplot(aes(rna_anno_2ndRound_level_2, log10_CPM_1p))
p + geom_violin(aes(
  fill = disease,
  # color = disease,
  group = paste(rna_anno_2ndRound_level_2, disease, subject)
),
draw_quantiles = c(0.5),
trim = TRUE,
scale = "width",
size = 0.2,
width = 0.8,
position = position_dodge(width = 0.8)
) +
  facet_wrap(~gene, ncol = 1, scales = "free_y") +
  ggtitle(selected_region) +
  scale_fill_manual(values = disease_color_panel) +
  scale_color_manual(values = lightness(disease_color_panel, scalefac(0.5))) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank()
  )
# ggsave("astro_genes_to_validate_MCX_ALS_vs_Control_level_2.pdf", device = cairo_pdf(), width = 3, height = 6, useDingbats = FALSE)
# ggsave("ATP2B2_ATP2B4_MCX_ALS_vs_Control_level_2_Astro.pdf", device = cairo_pdf(), width = 3, height = 6, useDingbats = FALSE)
# ggsave("ATP2B2_ATP2B4_mFCX_ALS_vs_Control_level_2_Astro.pdf", device = cairo_pdf(), width = 3, height = 6, useDingbats = FALSE)
# ggsave("ATP2B2_ATP2B4_MCX_ALS_vs_Control_level_2_Exc.pdf", device = cairo_pdf(), width = 8, height = 6, useDingbats = FALSE)
# ggsave("ATP2B2_ATP2B4_mFCX_ALS_vs_Control_level_2_Exc.pdf", device = cairo_pdf(), width = 8, height = 6, useDingbats = FALSE)

p <- to_plot_df %>% ggplot(aes(rna_anno_2ndRound_level_2, log10_CPM_1p))
p + geom_violin(aes(
  fill = disease,
  # color = disease,
  group = paste(rna_anno_2ndRound_level_2, disease)
),
draw_quantiles = c(0.5),
trim = TRUE,
scale = "width",
size = 0.2,
width = 0.8,
position = position_dodge(width = 0.8)
) +
  facet_wrap(~gene, ncol = 1, scales = "free_y", strip.position = "left") +
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
ggsave("./plots/DE_ALS_vs_Control_MCX_example_genes_Astro.pdf", device = cairo_pdf(), width = 4, height = 3, useDingbats = F)


p <- to_plot_df %>% ggplot(aes(rna_anno_2ndRound_level_3, log10_CPM_1p))
p + geom_violin(aes(fill = disease, color = disease, group = paste(rna_anno_2ndRound_level_3, disease, subject)),
  draw_quantiles = c(0.5),
  trim = TRUE,
  scale = "width",
  size = 0.2,
  width = 0.8,
  position = position_dodge(width = 0.8)
) +
  facet_wrap(~gene, ncol = 1, scales = "free_y") +
  ggtitle(selected_region) +
  scale_fill_manual(values = disease_color_panel) +
  scale_color_manual(values = lightness(disease_color_panel, scalefac(0.5))) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0)
  )
# ggsave("astro_genes_to_validate_MCX_ALS_vs_Control_level_3.pdf", device = cairo_pdf(), width = 8, height = 6, useDingbats = FALSE)
# ggsave("ATP2B2_ATP2B4_MCX_ALS_vs_Control_level_3_Astro.pdf", device = cairo_pdf(), width = 8, height = 6, useDingbats = FALSE)
# ggsave("ATP2B2_ATP2B4_mFCX_ALS_vs_Control_level_3_Astro.pdf", device = cairo_pdf(), width = 8, height = 6, useDingbats = FALSE)
# ggsave("ATP2B2_ATP2B4_mFCX_ALS_vs_Control_level_3_Exc.pdf", device = cairo_pdf(), width = 25, height = 6, useDingbats = FALSE)

# plot marker genes for clusters, show control samples only
# selected_region <- "MCX"
selected_region <- c("MCX", "mFCX")
selected_cell_types <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)
selected_genes <- c(
  "CUX2", "RORB", "THEMIS",
  "PVALB", "SST", "VIP", "LAMP5", "ADARB2",
  "AQP4", "VWF", "CSF1R", "MOBP", "PDGFRA"
)

selected_diseases <- "Control"
cell_idx <- which(meta_data$region %in% selected_region &
  meta_data$disease %in% selected_diseases &
  meta_data$rna_anno_2ndRound_level_2 %in% selected_cell_types)

meta_data_used <- meta_data %>%
  filter(
    region %in% selected_region,
    disease %in% selected_diseases,
    rna_anno_2ndRound_level_2 %in% selected_cell_types
  )

cpm_mat_used <- cpm_mat[selected_genes, cell_idx, drop = FALSE]
all.equal(colnames(cpm_mat_used), meta_data_used$cell_id)

df_genes_cpm <- cpm_mat_used %>%
  as.matrix() %>%
  t() %>%
  as_tibble() %>%
  rename_all(.funs = function(x) (paste0("gene_CPM_", x))) %>%
  bind_cols(meta_data_used)

to_plot_df <- df_genes_cpm %>%
  select(starts_with("gene_CPM_"), region, disease, subject, rna_anno_2ndRound_level_2, rna_anno_2ndRound_level_3) %>%
  pivot_longer(cols = starts_with("gene_CPM_"), names_to = "gene", values_to = "CPM") %>%
  mutate(
    gene = str_replace(gene, "^gene_CPM_", ""),
    gene = factor(gene, levels = selected_genes),
    log10_CPM_1p = log10(CPM + 1),
    rna_anno_2ndRound_level_2 = factor(rna_anno_2ndRound_level_2, levels = selected_cell_types)
  )

color_palette <- read_tsv("color_palette_level_2.txt")
colors <- color_palette$color
names(colors) <- color_palette$sub_cluster

p <- to_plot_df %>% ggplot(aes(rna_anno_2ndRound_level_2, log10_CPM_1p))
p + geom_violin(aes(
  fill = rna_anno_2ndRound_level_2,
  color = rna_anno_2ndRound_level_2,
  group = paste(rna_anno_2ndRound_level_2, disease, subject)
),
draw_quantiles = c(0.5),
trim = TRUE,
scale = "width",
size = 0.2,
width = 0.8,
position = position_dodge(width = 0.8)
) +
  facet_wrap(~ gene + region, ncol = 2, scales = "free_y", strip.position = "left") +
  # ggtitle(selected_region) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = lightness(colors, scalefac(0.5))) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
  )
ggsave("./plots/violin_markers_control_samples_level2_subclusters.pdf",
  device = cairo_pdf(),
  width = 10, height = 4, useDingbats = F
)

# plot marker genes for clusters as heatmap, show all samples by individual
selected_region <- c("MCX", "mFCX")
selected_cell_types <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)
selected_genes <- c(
  "CUX2", "RORB", "THEMIS",
  "PVALB", "SST", "VIP", "LAMP5", "ADARB2",
  "AQP4", "VWF", "CSF1R", "MOBP", "PDGFRA", "TBX18"
)

selected_diseases <- c("ALS", "FTD", "Control")
cell_idx <- which(meta_data$region %in% selected_region &
  meta_data$disease %in% selected_diseases &
  meta_data$rna_anno_2ndRound_level_2 %in% selected_cell_types)

meta_data_used <- meta_data %>%
  filter(
    region %in% selected_region,
    disease %in% selected_diseases,
    rna_anno_2ndRound_level_2 %in% selected_cell_types
  )

cpm_mat_used <- cpm_mat[selected_genes, cell_idx, drop = FALSE]
all.equal(colnames(cpm_mat_used), meta_data_used$cell_id)

df_genes_cpm <- cpm_mat_used %>%
  as.matrix() %>%
  t() %>%
  as_tibble() %>%
  rename_all(.funs = function(x) (paste0("gene_CPM_", x))) %>%
  bind_cols(meta_data_used)

to_plot_df <- df_genes_cpm %>%
  select(starts_with("gene_CPM_"), region, disease, subject, rna_anno_2ndRound_level_2, rna_anno_2ndRound_level_3) %>%
  pivot_longer(cols = starts_with("gene_CPM_"), names_to = "gene", values_to = "CPM") %>%
  mutate(
    gene = str_replace(gene, "^gene_CPM_", ""),
    gene = factor(gene, levels = selected_genes),
    log10_CPM_1p = log10(CPM + 1),
    rna_anno_2ndRound_level_2 = factor(rna_anno_2ndRound_level_2, levels = selected_cell_types)
  ) %>%
  group_by(region, disease, subject, rna_anno_2ndRound_level_2, gene) %>%
  summarize(median_log10_CPM_1p = median(log10_CPM_1p)) %>%
  ungroup() %>%
  mutate(group = paste(region, disease, subject, rna_anno_2ndRound_level_2, sep = "__"))

mat <- to_plot_df %>%
  select(gene, group, median_log10_CPM_1p) %>%
  pivot_wider(names_from = group, values_from = median_log10_CPM_1p, values_fill = NA) %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  as.matrix()

###
color_fun <- colorRamp2(seq(0, 4, length.out = 10), scico(10, palette = "lajolla", direction = -1))

region_color <- c("MCX" = "#432266", "mFCX" = "#FAA51B")
cell_type_color_df <- read_tsv("color_palette_level_2.txt") %>%
  mutate(sub_cluster = factor(sub_cluster, levels = selected_cell_types))
cell_type_color <- cell_type_color_df$color
names(cell_type_color) <- cell_type_color_df$sub_cluster

cell_class_color <- c(
  "Non_neuron" = "#9C482B",
  "Exc_neuron" = "#96BB45",
  "Inh_neuron" = "#718DC7"
)

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")

col_meta <- tibble(
  group = colnames(mat),
  group2 = colnames(mat)
) %>%
  separate(group2, c("region", "disease", "subject", "cell_type"), sep = "__") %>%
  mutate(
    cell_class = case_when(
      grepl("Exc_", cell_type) ~ "Exc_neuron",
      grepl("Inh_", cell_type) ~ "Inh_neuron",
      TRUE ~ "Non_neuron"
    ),
    cell_type = factor(cell_type, levels = selected_cell_types),
    disease = factor(disease, levels = c("ALS", "FTD", "Control")),
    region = factor(region, levels = c("MCX", "mFCX"))
  )

group_order <- col_meta %>%
  arrange(cell_type, disease, subject, region) %>%
  pull(group)


mat_col_anno <- HeatmapAnnotation(
  region = col_meta$region,
  disease = col_meta$disease,
  cell_class = col_meta$cell_class,
  cell_type = col_meta$cell_type,
  col = list(
    region = region_color,
    disease = disease_color,
    cell_class = cell_class_color,
    cell_type = cell_type_color
  )
)

# mat_row_anno <- rowAnnotation(
#   tau = gene_tau,
#   gene = anno_mark(
#     at = match(selected_DE_genes_sub$gene, rownames(mat)),
#     labels = selected_DE_genes_sub$gene
#     # at = match(selected_DE_genes$gene, rownames(mat)),
#     # labels = selected_DE_genes$gene
#   ),
#   col = list(tau = color_fun_tau)
# )

###

pdf(("./plots/marker_genes_by_individual_heatmap.pdf"), width = 8, height = 8)
ht <- Heatmap(
  matrix = mat,
  col = color_fun,
  name = "log10(CPM+1)",
  na_col = "grey",
  # rect_gp = gpar(col = "white", lwd = 2),
  # border_gp = gpar(col = "black", lty = 2),
  color_space = "LAB",
  show_row_names = TRUE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  # clustering_distance_rows = "euclidean",
  # clustering_method_rows = "ward.D2",
  # row_dend_reorder = TRUE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_order = group_order,
  # row_km = 7,
  # row_km_repeats = 100,
  # row_split = 8,
  # cluster_row_slices = FALSE,
  top_annotation = mat_col_anno,
  # right_annotation = mat_row_anno,
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

## log sessionInfo
sessionInfo()