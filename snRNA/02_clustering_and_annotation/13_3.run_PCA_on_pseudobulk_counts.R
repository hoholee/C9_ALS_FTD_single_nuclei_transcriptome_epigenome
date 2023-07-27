# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(DESeq2)
library(Hmisc)
library(ggpubr)
library(patchwork)

# read in pseudobulk grouping metadata
sample_info <- read_tsv("./pseudobulk_group_metadata.tsv")

# read the pseudobulk counts genereated from 13_2.get_pseudobulk_counts_by_individual.R
mat_summary_mm <- readRDS("snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_RNA_raw_count.rds")
mat_summary_mm_SCT <- readRDS("snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_SCT_raw_count.rds")

mat_cpm <- readRDS("snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_RNA_CPM.rds")
mat_cpm_SCT <- readRDS("snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_SCT_CPM.rds")

run_PCA_by_group <- function(selected_cell_type, selected_region, selected_diseases) {
  selected_samples <- sample_info %>%
    filter(
      cell_type == selected_cell_type,
      region == selected_region,
      disease %in% selected_diseases
    ) %>%
    pull(group)

  mat_sub <- mat_summary_mm[, selected_samples]

  sample_info_df <- sample_info %>%
    filter(group %in% selected_samples) %>%
    mutate(group = factor(group, levels = selected_samples)) %>%
    arrange(group) %>%
    as.data.frame() %>%
    column_to_rownames("group")

  dds <- DESeqDataSetFromMatrix(
    countData = mat_sub,
    colData = sample_info_df,
    design = ~ disease + sex
  )

  # Perform median ratio normalization
  dds_norm <- estimateSizeFactors(dds)
  norm_counts <- counts(dds_norm, normalized = TRUE)

  # Identify constant or zero-variance genes
  zero_variance_genes <- which(apply(norm_counts, 1, var) == 0)

  # Remove constant or zero-variance genes from the data
  norm_counts_filtered <- norm_counts[-zero_variance_genes, ]
  norm_counts_filtered_log <- log1p(norm_counts_filtered)

  # Perform PCA
  pca_res <- prcomp(t(norm_counts_filtered_log), scale = TRUE, center = TRUE)
  # pca_res <- prcomp(t(norm_counts_filtered_log), scale = FALSE, center = TRUE)

  # Convert PCA results and metadata into a data frame
  data.frame(pca_res$x, colData(dds_norm)) %>% as_tibble()
}

cell_type_list <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

cell_type_list_FTD <- c(
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")

# ALS vs. Control, motor cortex
res <- map_dfr(
  cell_type_list,
  run_PCA_by_group,
  selected_region = "MCX",
  selected_diseases = c("ALS", "Control")
) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_list))

# Create PCA plot with ellipses
res %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_MCX_PC1_vs_PC2.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_MCX_PC1_vs_PC2_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 5,
  useDingbats = FALSE
)
dev.off()

res %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_MCX_PC3_vs_PC4.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_MCX_PC3_vs_PC4_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 5,
  useDingbats = FALSE
)
dev.off()

res %>%
  ggplot(aes(x = PC5, y = PC6)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_MCX_PC5_vs_PC6.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_MCX_PC5_vs_PC6_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 5,
  useDingbats = FALSE
)
dev.off()

# only plot selected cell types, remove ellipses, plot PC 1-4
# selected_cell_types <- c("Astro", "Exc_superficial", "Exc_intermediate", "Exc_deep")
# selected_cell_type <- "Astro"
# selected_cell_type <- "Exc_superficial"
# selected_cell_type <- "Exc_intermediate"
# selected_cell_type <- "Exc_deep"

plot_PCA <- function(selected_cell_type) {
  res_sub <- res %>%
    filter(cell_type == selected_cell_type)
  p1 <- res_sub %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p2 <- res_sub %>%
    ggplot(aes(x = PC1, y = PC3)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p3 <- res_sub %>%
    ggplot(aes(x = PC1, y = PC4)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p4 <- res_sub %>%
    ggplot(aes(x = PC2, y = PC3)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p5 <- res_sub %>%
    ggplot(aes(x = PC2, y = PC4)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p6 <- res_sub %>%
    ggplot(aes(x = PC3, y = PC4)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )

  p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3, guides = "collect")
  ggsave(
    str_glue("./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_MCX_PC1to4_noEllipses_{selected_cell_type}.pdf"),
    device = cairo_pdf(),
    width = 4.5,
    height = 2,
    useDingbats = FALSE
  )
  dev.off()
}

walk(cell_type_list, plot_PCA)

# ALS vs. Control, frontal cortex
res <- map_dfr(
  cell_type_list,
  run_PCA_by_group,
  selected_region = "mFCX",
  selected_diseases = c("ALS", "Control")
) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_list))

# Create PCA plot with ellipses
res %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_mFCX_PC1_vs_PC2.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_mFCX_PC1_vs_PC2_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 5,
  useDingbats = FALSE
)
dev.off()

res %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_mFCX_PC3_vs_PC4.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_mFCX_PC3_vs_PC4_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 5,
  useDingbats = FALSE
)
dev.off()

res %>%
  ggplot(aes(x = PC5, y = PC6)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_mFCX_PC5_vs_PC6.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_mFCX_PC5_vs_PC6_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 5,
  useDingbats = FALSE
)
dev.off()

# only plot selected cell types, remove ellipses, plot PC 1-4
# selected_cell_types <- c("Astro", "Exc_superficial", "Exc_intermediate", "Exc_deep")
# selected_cell_type <- "Astro"
# selected_cell_type <- "Exc_superficial"
# selected_cell_type <- "Exc_intermediate"
# selected_cell_type <- "Exc_deep"

plot_PCA <- function(selected_cell_type) {
  res_sub <- res %>%
    filter(cell_type == selected_cell_type)
  p1 <- res_sub %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p2 <- res_sub %>%
    ggplot(aes(x = PC1, y = PC3)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p3 <- res_sub %>%
    ggplot(aes(x = PC1, y = PC4)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p4 <- res_sub %>%
    ggplot(aes(x = PC2, y = PC3)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p5 <- res_sub %>%
    ggplot(aes(x = PC2, y = PC4)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )
  p6 <- res_sub %>%
    ggplot(aes(x = PC3, y = PC4)) +
    geom_point(aes(color = disease, shape = sex), size = 1) +
    scale_color_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )

  p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3, guides = "collect")
  ggsave(
    str_glue("./plots/pseudobulk_PCA_by_cell_type_ALS_vs_Control_mFCX_PC1to4_noEllipses_{selected_cell_type}.pdf"),
    device = cairo_pdf(),
    width = 4.5,
    height = 2,
    useDingbats = FALSE
  )
  dev.off()
}

walk(cell_type_list, plot_PCA)


###############################################################################

# FTD vs. Control, motor cortex
res <- map_dfr(
  cell_type_list_FTD,
  run_PCA_by_group,
  selected_region = "MCX",
  selected_diseases = c("FTD", "Control")
) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_list_FTD))

# Create PCA plot with ellipses
res %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_MCX_PC1_vs_PC2.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_MCX_PC1_vs_PC2_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 2.5,
  useDingbats = FALSE
)
dev.off()

res %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_MCX_PC3_vs_PC4.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_MCX_PC3_vs_PC4_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 2.5,
  useDingbats = FALSE
)
dev.off()

res %>%
  ggplot(aes(x = PC5, y = PC6)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_MCX_PC5_vs_PC6.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_MCX_PC5_vs_PC6_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 2.5,
  useDingbats = FALSE
)
dev.off()

# only plot selected cell types, remove ellipses, plot PC 1-4
# selected_cell_type <- "Astro"
# selected_cell_type <- "Endo"
# selected_cell_type <- "Micro"
# selected_cell_type <- "Oligo"
# selected_cell_type <- "OPC"
selected_cell_type <- "VLMC"

res_sub <- res %>%
  filter(cell_type == selected_cell_type)
p1 <- res_sub %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p2 <- res_sub %>%
  ggplot(aes(x = PC1, y = PC3)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p3 <- res_sub %>%
  ggplot(aes(x = PC1, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p4 <- res_sub %>%
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p5 <- res_sub %>%
  ggplot(aes(x = PC2, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p6 <- res_sub %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3, guides = "collect")
ggsave(
  str_glue("./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_MCX_PC1to4_noEllipses_{selected_cell_type}.pdf"),
  device = cairo_pdf(),
  width = 4.5,
  height = 2,
  useDingbats = FALSE
)
dev.off()

# FTD vs. Control, frontal cortex
res <- map_dfr(
  cell_type_list_FTD,
  run_PCA_by_group,
  selected_region = "mFCX",
  selected_diseases = c("FTD", "Control")
) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_list_FTD))

# Create PCA plot with ellipses
res %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_mFCX_PC1_vs_PC2.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_mFCX_PC1_vs_PC2_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 2.5,
  useDingbats = FALSE
)
dev.off()

res %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_mFCX_PC3_vs_PC4.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_mFCX_PC3_vs_PC4_noEllipses.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 2.5,
  useDingbats = FALSE
)
dev.off()

res %>%
  ggplot(aes(x = PC5, y = PC6)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed", aes(color = disease)) +
  facet_wrap(~cell_type, ncol = 3) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
ggsave(
  # "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_mFCX_PC5_vs_PC6.pdf",
  "./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_mFCX_PC5_vs_PC6_noEllispes.pdf",
  device = cairo_pdf(),
  width = 4.5,
  height = 2.5,
  useDingbats = FALSE
)
dev.off()

# only plot selected cell types, remove ellipses, plot PC 1-4
# selected_cell_type <- "Astro"
# selected_cell_type <- "Endo"
# selected_cell_type <- "Micro"
# selected_cell_type <- "Oligo"
# selected_cell_type <- "OPC"
selected_cell_type <- "VLMC"

res_sub <- res %>%
  filter(cell_type == selected_cell_type)
p1 <- res_sub %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p2 <- res_sub %>%
  ggplot(aes(x = PC1, y = PC3)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p3 <- res_sub %>%
  ggplot(aes(x = PC1, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p4 <- res_sub %>%
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p5 <- res_sub %>%
  ggplot(aes(x = PC2, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )
p6 <- res_sub %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  scale_color_manual(values = disease_color) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "lines")
  )

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3, guides = "collect")
ggsave(
  str_glue("./plots/pseudobulk_PCA_by_cell_type_FTD_vs_Control_mFCX_PC1to4_noEllipses_{selected_cell_type}.pdf"),
  device = cairo_pdf(),
  width = 4.5,
  height = 2,
  useDingbats = FALSE
)
dev.off()

# compute pair-wise correlation

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  tibble(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

expression_cutoff <- 10
percent_sample <- 0.5

get_cor_by_group <- function(selected_cell_type, selected_region, selected_diseases) {
  selected_samples <- sample_info %>%
    filter(
      cell_type == selected_cell_type,
      region == selected_region,
      disease %in% selected_diseases
    ) %>%
    pull(group)

  mat_cpm_SCT_sub <- mat_cpm_SCT[, selected_samples]
  mat_cpm_SCT_sub_log <- log10(mat_cpm_SCT_sub + 1)

  # Identify constant or zero-variance genes
  zero_variance_genes <- which(apply(mat_cpm_SCT_sub_log, 1, var) == 0)
  mat_cpm_SCT_sub_log <- mat_cpm_SCT_sub_log[-zero_variance_genes, ]

  # Identify genes expressed in more than `percent_sample` samples,
  # with expression cutoff of `expression_cutoff`
  num_samples <- ceiling(ncol(mat_cpm_SCT_sub_log) * percent_sample)
  filtered_genes <- apply(
    mat_cpm_SCT_sub_log,
    1,
    function(x) sum(x > log10(expression_cutoff + 1)) >= num_samples
  )
  mat_cpm_SCT_sub_log <- mat_cpm_SCT_sub_log[filtered_genes, ]
  # message(dim(mat_cpm_SCT_sub_log)[1])

  # sample_info_df <- sample_info %>%
  #   filter(group %in% selected_samples) %>%
  #   mutate(group = factor(group, levels = selected_samples)) %>%
  #   arrange(group) %>%
  #   as.data.frame() %>%
  #   column_to_rownames("group")

  cor_mat <- rcorr(as.matrix(mat_cpm_SCT_sub_log), type = "pearson")

  flattenCorrMatrix(cor_mat$r, cor_mat$P) %>%
    mutate(
      cell_type = selected_cell_type,
      region = selected_region
    )
}

cor_res_sub <- get_cor_by_group("Astro", "MCX", c("ALS", "Control")) %>%
  mutate(
    disease_row = if_else(grepl("ALS", row), "ALS", "Control"),
    disease_column = if_else(grepl("ALS", column), "ALS", "Control"),
    group = if_else(disease_row == disease_column, "within", "between"),
    group2 = paste(disease_row, disease_column, sep = " vs. ")
  )

p <- cor_res_sub %>% ggplot(aes(group, cor))

p +
  geom_boxplot() +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("within", "between")),
    aes(label = ..p.signif..)
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank()
  )
ggsave(
  "./plots/pseudobulk_logCPM_pearson_correlation_Astro_MCX_within_disease_vs_between_disease.pdf",
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
dev.off()

p <- cor_res_sub %>% ggplot(aes(group2, cor))

p +
  geom_boxplot() +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("ALS vs. ALS", "Control vs. Control"),
      c("ALS vs. ALS", "ALS vs. Control"),
      c("Control vs. Control", "ALS vs. Control")
    ),
    aes(label = ..p.signif..)
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank()
  )
ggsave(
  "./plots/pseudobulk_logCPM_pearson_correlation_Astro_MCX_within_disease_vs_between_disease_split.pdf",
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
dev.off()


#####

get_cor_by_group <- function(selected_region, selected_diseases) {
  selected_samples <- sample_info %>%
    filter(
      region == selected_region,
      disease %in% selected_diseases
    ) %>%
    pull(group)

  mat_cpm_SCT_sub <- mat_cpm_SCT[, selected_samples]
  mat_cpm_SCT_sub_log <- log10(mat_cpm_SCT_sub + 1)

  # Identify constant or zero-variance genes
  zero_variance_genes <- which(apply(mat_cpm_SCT_sub_log, 1, var) == 0)
  mat_cpm_SCT_sub_log <- mat_cpm_SCT_sub_log[-zero_variance_genes, ]

  # Identify genes expressed in more than `percent_sample` samples,
  # with expression cutoff of `expression_cutoff`
  num_samples <- ceiling(ncol(mat_cpm_SCT_sub_log) * percent_sample)
  filtered_genes <- apply(
    mat_cpm_SCT_sub_log,
    1,
    function(x) sum(x > log10(expression_cutoff + 1)) >= num_samples
  )
  mat_cpm_SCT_sub_log <- mat_cpm_SCT_sub_log[filtered_genes, ]
  # message(dim(mat_cpm_SCT_sub_log)[1])

  # sample_info_df <- sample_info %>%
  #   filter(group %in% selected_samples) %>%
  #   mutate(group = factor(group, levels = selected_samples)) %>%
  #   arrange(group) %>%
  #   as.data.frame() %>%
  #   column_to_rownames("group")

  cor_mat <- rcorr(as.matrix(mat_cpm_SCT_sub_log), type = "pearson")

  flattenCorrMatrix(cor_mat$r, cor_mat$P) %>%
    mutate(
      region = selected_region
    )
}

cor_res_sub <- get_cor_by_group("MCX", c("ALS", "Control")) %>%
  mutate(
    row2 = row,
    column2 = column,
  ) %>%
  separate(
    row2,
    c("region_row", "disease_row", "subject_row", "cell_type_row"),
    sep = " "
  ) %>%
  separate(
    column2,
    c("region_column", "disease_column", "subject_column", "cell_type_column"),
    sep = " "
  ) %>%
  mutate(
    group = case_when(
      cell_type_row != cell_type_column ~ "between cell types",
      disease_row != disease_column ~ "between diseases",
      cell_type_row == cell_type_column & disease_row == disease_column ~ "between subjects"
    )
  )

p <- cor_res_sub %>% ggplot(aes(group, cor))

p +
  geom_boxplot() +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("between cell types", "between diseases"),
      c("between cell types", "between subjects"),
      c("between diseases", "between subjects")
    ),
    aes(label = ..p.signif..)
  ) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  "./plots/pseudobulk_logCPM_pearson_correlation_MCX_ALS_and_Control_effect_celltype_disease_subject.pdf",
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
dev.off()




#####
mat_cpm_SCT_log <- log10(mat_cpm_SCT + 1)

cor_mat <- rcorr(as.matrix(mat_cpm_SCT_log), type = "pearson")

cor_res <- flattenCorrMatrix(cor_mat$r, cor_mat$P)

cell_type_list <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

cell_type_list_FTD <- c(
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")

# ALS vs. Control, motor cortex
res <- map_dfr(
  cell_type_list,
  run_PCA_by_group,
  selected_region = "MCX",
  selected_diseases = c("ALS", "Control")
) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_list))
## log sessionInfo
sessionInfo()
