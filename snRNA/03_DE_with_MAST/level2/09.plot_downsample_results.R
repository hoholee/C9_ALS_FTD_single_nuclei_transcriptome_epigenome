# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ggrepel)
library(directlabels)

sem <- function(x) {
  sd(x, na.rm = T) / sqrt(length(na.omit(x)))
}
ci <- function(x, conf.interval = 0.95) {
  N <- sum(!is.na(x))
  sem <- sd(x, na.rm = T) / sqrt(N)
  qt(1 - (1 - conf.interval) / 2, N - 1) * sem
}

output_plots_path <- "/cndd2/junhao/ALS_FTD_singleCell/snRNA_postCellBender/plots/"

read_MAST_res <- function(selected_celltype, selected_region, selected_comparison, selected_seed) {
  read_tsv(paste0(
    "./results_downsample/MAST_res_level2_", selected_celltype, "_", selected_region,
    "_", selected_comparison, "_seed_", selected_seed, ".txt"
  )) %>%
    mutate(
      celltype = selected_celltype,
      region = selected_region,
      comparison = selected_comparison,
      seed = selected_seed
    )
}

read_MAST_res_obs <- function(selected_celltype, selected_region, selected_comparison) {
  read_tsv(paste0(
    "./results/MAST_res_level2_", selected_celltype, "_", selected_region,
    "_", selected_comparison, ".txt"
  )) %>%
    mutate(
      celltype = selected_celltype,
      region = selected_region,
      comparison = selected_comparison
    )
}

# ALS vs. Control
param <- expand_grid(
  selected_celltype = c(
    "Astro", "Endo", "Micro", "OPC", "Oligo", "VLMC",
    "Exc_superficial", "Exc_intermediate", "Exc_deep", "Exc_unknown",
    "Inh_LAMP5", "Inh_VIP", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
  ),
  selected_region = c("MCX", "mFCX"),
  selected_comparison = "ALS_vs_Control",
  selected_seed = c(4, 8, 15, 16, 23, 42, 78, 173, 666, 2021)
)

param_obs <- expand_grid(
  selected_celltype = c(
    "Astro", "Endo", "Micro", "OPC", "Oligo", "VLMC",
    "Exc_superficial", "Exc_intermediate", "Exc_deep", "Exc_unknown",
    "Inh_LAMP5", "Inh_VIP", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
  ),
  selected_region = c("MCX", "mFCX"),
  selected_comparison = "ALS_vs_Control"
)

merged_res <- pmap_dfr(param, read_MAST_res)
merged_res_obs <- pmap_dfr(param_obs, read_MAST_res_obs)

FDR_threshold <- 0.05
FC_threshold <- 1.2

DE_count <- merged_res %>%
  filter(
    fdr < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    ci.hi * ci.lo > 0,
    abs(model_log2FC - avg_log2FC) < 2
  ) %>%
  group_by(celltype, region, seed) %>%
  count() %>%
  ungroup()

DE_count_obs <- merged_res_obs %>%
  filter(
    fdr < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    ci.hi * ci.lo > 0,
    abs(model_log2FC - avg_log2FC) < 2
  ) %>%
  group_by(celltype, region) %>%
  count() %>%
  ungroup()

p <- DE_count %>% ggplot(aes(celltype, n))
p + geom_boxplot(outlier.shape = NA) +
  facet_wrap(~region) +
  ylab("# DE genes") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(output_plots_path, "DE_cellBender_corrected_ALS_vs_Control_downsample_FC_", FC_threshold, "_num_DE_genes_by_level2_anno.pdf"),
  device = cairo_pdf(), width = 6, height = 4, useDingbats = F
)

p <- DE_count %>% ggplot(aes(celltype, n))
p + geom_boxplot(outlier.shape = NA) +
  geom_point(data = DE_count_obs, shape = 2) +
  facet_wrap(~region) +
  ylab("# DE genes") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(output_plots_path, "DE_cellBender_corrected_ALS_vs_Control_downsample_FC_", FC_threshold, "_num_DE_genes_by_level2_anno_add_obs.pdf"),
  device = cairo_pdf(), width = 6, height = 4, useDingbats = F
)


#
DE_count_summary <- DE_count %>%
  group_by(celltype, region) %>%
  summarise(
    n_mean = mean(n),
    n_sem = sem(n),
    n_ci = ci(n)
  ) %>%
  ungroup()

DE_count_sem_summary <- DE_count_summary %>%
  select(celltype, region, n_sem) %>%
  pivot_wider(names_from = "region", values_from = "n_sem", names_prefix = "sem_")

DE_count_summary_to_plot <- DE_count_summary %>%
  select(celltype, region, n_mean) %>%
  pivot_wider(names_from = "region", values_from = "n_mean", names_prefix = "mean_") %>%
  left_join(DE_count_sem_summary)

colors <- read_tsv("../../color_palette_level_2.txt")
celltype_colors <- colors$color
names(celltype_colors) <- colors$sub_cluster

p <- DE_count_summary_to_plot %>%
  ggplot(aes(mean_MCX, mean_mFCX))
p +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(aes(color = celltype)) +
  geom_errorbar(aes(ymin = mean_mFCX - sem_mFCX, ymax = mean_mFCX + sem_mFCX, color = celltype), width = 3) +
  geom_errorbarh(aes(xmin = mean_MCX - sem_MCX, xmax = mean_MCX + sem_MCX, color = celltype), height = 3) +
  geom_text_repel(aes(label = celltype)) +
  xlab("# of DE genes in motor cortex (downsampled)") +
  ylab("# of DE genes in medial frontal cortex (downsampled)") +
  scale_color_manual(values = celltype_colors, name = "Major cell types") +
  coord_equal() +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())
ggsave(paste0(output_plots_path, "DE_cellBender_corrected_ALS_vs_Control_downsample_FC_", FC_threshold, "_num_DE_genes_MCX_vs_mFCX.pdf"),
  device = cairo_pdf(),
  width = 4, height = 6, useDingbats = F
)


#
get_DE_count_by_FC_threshold <- function(selected_FC_threshold, selected_FDR_threshold = 0.05) {
  merged_res %>%
    filter(
      fdr < selected_FDR_threshold,
      abs(model_log2FC) > log2(selected_FC_threshold),
      conv_C == TRUE,
      conv_D == TRUE,
      ci.hi * ci.lo > 0,
      abs(model_log2FC - avg_log2FC) < 2
    ) %>%
    group_by(celltype, region, seed) %>%
    count() %>%
    ungroup() %>%
    mutate(FC_threshold = selected_FC_threshold)
}

FC_threshold_list <- seq(from = 1, to = 3.5, by = 0.5)

DE_count_by_FC_threshold <- map_dfr(FC_threshold_list, get_DE_count_by_FC_threshold, selected_FDR_threshold = 0.05)

DE_count_by_FC_threshold_to_plot <- DE_count_by_FC_threshold %>%
  group_by(celltype, region, FC_threshold) %>%
  summarise(mean_n = mean(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = "region", values_from = "mean_n", names_prefix = "mean_") %>%
  mutate(num_DE_diff = mean_MCX - mean_mFCX)

p <- DE_count_by_FC_threshold_to_plot %>% ggplot(aes(FC_threshold, num_DE_diff))
p <- p +
  geom_vline(xintercept = FC_threshold, linetype = 2, color = "darkgrey") +
  geom_line(aes(color = celltype)) +
  geom_point(aes(color = celltype)) +
  scale_color_manual(values = celltype_colors) +
  ylab("Differences of number DE genes detected in the downsample test\n(motor cortex - medial frontol cortex)") +
  xlab("Fold-change threshold\n(ALS vs. Control)") +
  coord_cartesian(xlim = c(0.5, 3.5)) +
  theme_bw(base_size = 8, base_family = "Helvetica")

direct.label.ggplot(p, method = "first.points")
ggsave(paste0(output_plots_path, "DE_cellBender_corrected_ALS_vs_Control_downsample_diff_num_MCX_vs_mFCX_by_FC_threshold.pdf"),
  device = cairo_pdf(),
  width = 6, height = 4, useDingbats = F
)


# FTD vs. Control
param <- expand_grid(
  selected_celltype = c(
    "Astro", "Endo", "Micro", "OPC", "Oligo", "VLMC"
    # "Exc_superficial", "Exc_intermediate", "Exc_deep", "Exc_unknown",
    # "Inh_LAMP5", "Inh_VIP", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
  ),
  selected_region = c("MCX", "mFCX"),
  selected_comparison = "FTD_vs_Control",
  selected_seed = c(4, 8, 15, 16, 23, 42, 78, 173, 666, 2021)
)

param_obs <- expand_grid(
  selected_celltype = c(
    "Astro", "Endo", "Micro", "OPC", "Oligo", "VLMC"
    # "Exc_superficial", "Exc_intermediate", "Exc_deep", "Exc_unknown",
    # "Inh_LAMP5", "Inh_VIP", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
  ),
  selected_region = c("MCX", "mFCX"),
  selected_comparison = "FTD_vs_Control"
)

merged_res <- pmap_dfr(param, read_MAST_res)
merged_res_obs <- pmap_dfr(param_obs, read_MAST_res_obs)

FDR_threshold <- 0.05
FC_threshold <- 1.2

DE_count <- merged_res %>%
  filter(
    fdr < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    ci.hi * ci.lo > 0,
    abs(model_log2FC - avg_log2FC) < 2
  ) %>%
  group_by(celltype, region, seed) %>%
  count() %>%
  ungroup()

DE_count_obs <- merged_res_obs %>%
  filter(
    fdr < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    ci.hi * ci.lo > 0,
    abs(model_log2FC - avg_log2FC) < 2
  ) %>%
  group_by(celltype, region) %>%
  count() %>%
  ungroup()

p <- DE_count %>% ggplot(aes(celltype, n))
p + geom_boxplot(outlier.shape = NA) +
  facet_wrap(~region) +
  ylab("# DE genes") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(output_plots_path, "DE_cellBender_corrected_FTD_vs_Control_downsample_FC_", FC_threshold, "_num_DE_genes_by_level2_anno.pdf"),
  device = cairo_pdf(), width = 6, height = 4, useDingbats = F
)

p <- DE_count %>% ggplot(aes(celltype, n))
p + geom_boxplot(outlier.shape = NA) +
  geom_point(data = DE_count_obs, shape = 2) +
  facet_wrap(~region) +
  ylab("# DE genes") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(output_plots_path, "DE_cellBender_corrected_FTD_vs_Control_downsample_FC_", FC_threshold, "_num_DE_genes_by_level2_anno_add_obs.pdf"),
  device = cairo_pdf(), width = 6, height = 4, useDingbats = F
)

#
DE_count_summary <- DE_count %>%
  group_by(celltype, region) %>%
  summarise(
    n_mean = mean(n),
    n_sem = sem(n),
    n_ci = ci(n)
  ) %>%
  ungroup()

DE_count_sem_summary <- DE_count_summary %>%
  select(celltype, region, n_sem) %>%
  pivot_wider(names_from = "region", values_from = "n_sem", names_prefix = "sem_")

DE_count_summary_to_plot <- DE_count_summary %>%
  select(celltype, region, n_mean) %>%
  pivot_wider(names_from = "region", values_from = "n_mean", names_prefix = "mean_") %>%
  left_join(DE_count_sem_summary)

colors <- read_tsv("../../color_palette_level_2.txt")
celltype_colors <- colors$color
names(celltype_colors) <- colors$sub_cluster

p <- DE_count_summary_to_plot %>%
  ggplot(aes(mean_MCX, mean_mFCX))
p +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(aes(color = celltype), size = 1) +
  geom_errorbar(
    aes(
      ymin = mean_mFCX - sem_mFCX,
      ymax = mean_mFCX + sem_mFCX,
      color = celltype
    ),
    # width = 12
    width = 8
  ) +
  geom_errorbarh(
    aes(
      xmin = mean_MCX - sem_MCX,
      xmax = mean_MCX + sem_MCX,
      color = celltype
    ),
    # height = 12
    height = 8
  ) +
  geom_text_repel(aes(label = celltype)) +
  xlab("# of DE genes in motor cortex (downsampled)") +
  ylab("# of DE genes in medial frontal cortex (downsampled)") +
  scale_color_manual(values = celltype_colors, name = "Major cell types") +
  coord_equal() +
  # xlim(0, 650) +
  # ylim(0, 650) +
  xlim(0, 300) +
  ylim(0, 300) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank())
# ggsave(paste0(output_plots_path, "DE_cellBender_corrected_FTD_vs_Control_downsample_FC_", FC_threshold, "_num_DE_genes_MCX_vs_mFCX.pdf"),
#   device = cairo_pdf(),
#   width = 4, height = 6, useDingbats = F
# )
ggsave(paste0(output_plots_path, "DE_cellBender_corrected_FTD_vs_Control_downsample_FC_", FC_threshold, "_num_DE_genes_MCX_vs_mFCX_rmNeu.pdf"),
  device = cairo_pdf(),
  width = 4, height = 6, useDingbats = F
)


#
get_DE_count_by_FC_threshold <- function(selected_FC_threshold, selected_FDR_threshold = 0.05) {
  merged_res %>%
    filter(
      fdr < selected_FDR_threshold,
      abs(model_log2FC) > log2(selected_FC_threshold),
      conv_C == TRUE,
      conv_D == TRUE,
      ci.hi * ci.lo > 0,
      abs(model_log2FC - avg_log2FC) < 2
    ) %>%
    group_by(celltype, region, seed) %>%
    count() %>%
    ungroup() %>%
    mutate(FC_threshold = selected_FC_threshold)
}

FC_threshold_list <- seq(from = 1, to = 3.5, by = 0.5)

DE_count_by_FC_threshold <- map_dfr(FC_threshold_list, get_DE_count_by_FC_threshold, selected_FDR_threshold = 0.05)

DE_count_by_FC_threshold_to_plot <- DE_count_by_FC_threshold %>%
  group_by(celltype, region, FC_threshold) %>%
  summarise(mean_n = mean(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = "region", values_from = "mean_n", names_prefix = "mean_") %>%
  mutate(num_DE_diff = mean_MCX - mean_mFCX)

p <- DE_count_by_FC_threshold_to_plot %>% ggplot(aes(FC_threshold, num_DE_diff))
p <- p +
  geom_vline(xintercept = FC_threshold, linetype = 2, color = "darkgrey") +
  geom_line(aes(color = celltype)) +
  geom_point(aes(color = celltype)) +
  scale_color_manual(values = celltype_colors) +
  ylab("Differences of number DE genes detected in the downsample test\n(motor cortex - medial frontol cortex)") +
  xlab("Fold-change threshold\n(FTD vs. Control)") +
  coord_cartesian(xlim = c(0.5, 3.5)) +
  theme_bw(base_size = 8, base_family = "Helvetica")

direct.label.ggplot(p, method = "first.points")
ggsave(paste0(output_plots_path, "DE_cellBender_corrected_FTD_vs_Control_downsample_diff_num_MCX_vs_mFCX_by_FC_threshold.pdf"),
  device = cairo_pdf(),
  width = 6, height = 4, useDingbats = F
)