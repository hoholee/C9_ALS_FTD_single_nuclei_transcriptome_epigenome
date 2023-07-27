# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(scales)
library(future)
library(ggbeeswarm)
library(ggpubr)
library(shades)
library(ggrastr)

# enable parallelization
plan("multisession", workers = 8)
options(future.globals.maxSize = 2000 * 1024^2)

## read the raw seurat object
data_obj <- readRDS("./seurat_objects/snRNA_cellBender_corrected_SeuratV4_object.rds")

## quick check on the metadata
head(data_obj@meta.data)
meta_data <- data_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

## add scrublet doublet score
doublet_score <- read_tsv("/cndd2/junhao/ALS_FTD_singleCell/run_cellBender_on_raw_snRNA/run_scrublet/doublet_score_summary_postCellBender_callThreshold_0.2.txt")

meta_data_add_doubletScores <- meta_data %>%
  mutate(barcode = str_replace(cell_id, "-1$", "")) %>%
  left_join(doublet_score, by = c("barcode"))

data_obj$doublet_scores <- meta_data_add_doubletScores$doublet_scores
data_obj$predicted_doublets <- meta_data_add_doubletScores$predicted_doublets

## filter cells
# DO NOT set upper bound for number of genes detected here, since some neurons (eg. Betz cells) will have lots of genes expressed (>10k)
# nFeature_upper <- 8000
# nFeature_upper <- quantile(data_obj@meta.data$nFeature_RNA, 0.95)
# nFeature_upper
# the BICCN paper use these standard for snRNA datasets:
# neurons with fewer than 1000 detected genes and non-neuronal cells with fewer than 500 detected genes were removed
# here I set it to 500 for all cells before we get the annotations
nFeature_lower <- 500
# set an upper bound for the %reads mapped to the mitochondrial genome to remove low-quality / dying cells
percent_mt_upper <- 1
# percent_mt_upper <- 5
# percent_mt_upper <- 10
## manually set the threshold of `doublet_scores` to 0.2
doublet_scores_upper <- 0.2

data_obj_sub <- subset(data_obj,
  subset = nFeature_RNA > nFeature_lower & percent_mt < percent_mt_upper & doublet_scores <= doublet_scores_upper
)

meta_data_sub <- data_obj_sub@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

raw_cell_count <- data_obj@meta.data %>%
  group_by(orig.ident) %>%
  count() %>%
  ungroup() %>%
  rename(total_cell_count = n)
filtered_cell_count <- data_obj_sub@meta.data %>%
  group_by(orig.ident) %>%
  count() %>%
  ungroup() %>%
  rename(pass_QC_cell_count = n)
cell_count_merge <- raw_cell_count %>%
  left_join(filtered_cell_count) %>%
  mutate(
    pass_QC_cell_count = if_else(is.na(pass_QC_cell_count), 0L, pass_QC_cell_count),
    yield_rate = pass_QC_cell_count / total_cell_count
  )

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")
disease_color_dark <- lightness(disease_color, scalefac(0.4)) %>%
  saturation(scalefac(0.8)) %>%
  as.vector()
names(disease_color_dark) <- names(disease_color)
fill_colors <- c(
  disease_color_dark["ALS"],
  disease_color["ALS"],
  disease_color_dark["FTD"],
  disease_color["FTD"],
  disease_color_dark["Control"],
  disease_color["Control"]
)
names(fill_colors) <- c(
  "ALS fail", "ALS pass",
  "FTD fail", "FTD pass",
  "Control fail", "Control pass"
)

p <- cell_count_merge %>%
  mutate(fail_QC_cell_count = total_cell_count - pass_QC_cell_count) %>%
  select(-yield_rate, -total_cell_count) %>%
  gather(QC_type, cell_count, -orig.ident) %>%
  mutate(group = orig.ident) %>%
  separate(orig.ident, c("region", "disease", "subject"), sep = "_") %>%
  mutate(
    disease = factor(disease, levels = c("ALS", "FTD", "Control")),
    region = factor(region, levels = c("MCX", "mFCX")),
    fill_color = paste(disease, str_replace(QC_type, "_QC_cell_count", ""), sep = " "),
    fill_color = factor(fill_color,
      levels = c(
        "ALS fail",
        "ALS pass",
        "FTD fail",
        "FTD pass",
        "Control fail",
        "Control pass"
      )
    )
  ) %>%
  ggplot(aes(subject, cell_count))

p + geom_bar(aes(fill = fill_color),
  stat = "identity",
  position = position_stack(),
  width = 0.8
) +
  scale_fill_manual(values = fill_colors, name = "QC") +
  facet_grid(~ region + disease,
    scales = "free_x",
    drop = TRUE,
    space = "free",
    switch = "x"
  ) +
  ylab("Number of nuclei") +
  xlab("Sample") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    # panel.spacing = unit(0.5, "lines"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  )
ggsave(
  paste0(
    "./plots/QC_all_samples_cellBender_corrected_postFilter_nFeature_", nFeature_lower,
    "_percentMT_", percent_mt_upper, "_doubletScore_", doublet_scores_upper, "_yieldRate.pdf"
  ),
  device = cairo_pdf(), width = 4, height = 3, useDingbats = F
)

disease_color_dark2 <- lightness(disease_color, scalefac(0.8)) %>%
  # saturation(scalefac(0.8)) %>%
  as.vector()
names(disease_color_dark2) <- names(disease_color)

p <- cell_count_merge %>%
  mutate(sample = orig.ident) %>%
  separate(orig.ident, c("region", "disease", "subject"), sep = "_") %>%
  mutate(
    region = factor(region, levels = c("MCX", "mFCX")),
    disease = factor(disease, levels = c("ALS", "FTD", "Control"))
  ) %>%
  ggplot(aes(disease, 100 * yield_rate))

p +
  geom_boxplot(aes(color = disease), width = 0.5, outlier.shape = NA) +
  geom_beeswarm(aes(color = disease), size = 1, cex = 4) +
  scale_color_manual(values = disease_color_dark2) +
  facet_wrap(~region, strip.position = "bottom") +
  stat_compare_means(aes(label = ..p.signif..),
    method = "t.test",
    comparisons = list(c("ALS", "Control"), c("FTD", "Control"), c("ALS", "FTD"))
  ) +
  xlab("Disease") +
  ylab("Nuclei passing QC (%)") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "none"
  )
ggsave(
  paste0(
    "./plots/QC_all_samples_cellBender_corrected_postFilter_nFeature_", nFeature_lower,
    "_percentMT_", percent_mt_upper, "_doubletScore_", doublet_scores_upper, "_yieldRateComparison.pdf"
  ),
  device = cairo_pdf(), width = 2, height = 2, useDingbats = F
)

# plot cumulative distribution of the three criteria used to filter cells
# number of detected genes, percent mitochrondia reads and doublet scores
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

to_plot_df <- meta_data_add_doubletScores %>%
  mutate(
    region = factor(region, levels = c("MCX", "mFCX")),
    disease = factor(disease, levels = c("ALS", "FTD", "Control"))
  )

# p <- to_plot_df %>% ggplot(aes(log10(nFeature_RNA)))
# p +
#   rasterize(stat_ecdf(aes(color = subject), pad = FALSE), dpi = 600) +
#   geom_vline(
#     xintercept = log10(nFeature_lower),
#     color = "darkgrey",
#     linetype = 2
#   ) +
#   annotate(
#     geom = "rect",
#     xmin = -Inf,
#     xmax = log10(nFeature_lower),
#     ymin = -Inf,
#     ymax = Inf,
#     alpha = 0.3,
#     fill = "darkgrey"
#   ) +
#   annotation_logticks(
#     base = 10,
#     sides = "b",
#     short = unit(0.05, "cm"),
#     mid = unit(0.1, "cm"),
#     long = unit(0.15, "cm"),
#     size = 0.3
#   ) +
#   facet_grid(
#     region ~ disease,
#     scales = "fixed"
#   ) +
#   scale_color_manual(
#     values = individual_color,
#     name = "Donor ID",
#     guide = guide_legend(ncol = 3)
#   ) +
#   xlab("Number of detected genes") +
#   ylab("Cumulative proportion") +
#   scale_x_continuous(
#     breaks = c(3, 4),
#     labels = scales::math_format(10^.x)
#   ) +
#   theme_bw(base_size = 8, base_family = "Helvetica") +
#   theme(
#     panel.grid.minor = element_blank(),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     legend.position = "bottom",
#     legend.key.size = unit(0.5, "lines"),
#     legend.text = element_text(size = 6)
#   )
# ggsave(
#   "./plots/preQC_cumulative_distribution_num_detected_genes_by_individual.pdf",
#   device = cairo_pdf(),
#   width = 3, height = 3, useDingbats = FALSE
# )

p <- to_plot_df %>% ggplot(aes(percent_mt))
p +
  rasterize(stat_ecdf(aes(color = subject), pad = FALSE), dpi = 600) +
  geom_vline(xintercept = percent_mt_upper, color = "darkgrey", linetype = 2) +
  annotate(
    geom = "rect",
    xmin = percent_mt_upper,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.3,
    fill = "darkgrey"
  ) +
  facet_grid(
    region ~ disease,
    scales = "fixed"
  ) +
  scale_color_manual(
    values = individual_color,
    name = "Donor ID",
    guide = guide_legend(ncol = 3)
  ) +
  # coord_cartesian(xlim = c(0, 5)) +
  coord_cartesian(xlim = c(0, 20)) +
  xlab("% mitochrondia reads") +
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
  "./plots/preQC_cumulative_distribution_percent_mt_by_individual_xExtend20.pdf",
  device = cairo_pdf(),
  width = 3, height = 3, useDingbats = FALSE
)

# p <- to_plot_df %>% ggplot(aes(doublet_scores))
# p +
#   rasterize(stat_ecdf(aes(color = subject), pad = FALSE), dpi = 600) +
#   geom_vline(xintercept = doublet_scores_upper, color = "darkgrey", linetype = 2) +
#   annotate(
#     geom = "rect",
#     xmin = doublet_scores_upper,
#     xmax = Inf,
#     ymin = -Inf,
#     ymax = Inf,
#     alpha = 0.3,
#     fill = "darkgrey"
#   ) +
#   facet_grid(
#     region ~ disease,
#     scales = "fixed"
#   ) +
#   scale_color_manual(
#     values = individual_color,
#     name = "Donor ID",
#     guide = guide_legend(ncol = 3)
#   ) +
#   coord_cartesian(xlim = c(0, 0.5)) +
#   xlab("Doublet score") +
#   ylab("Cumulative proportion") +
#   theme_bw(base_size = 8, base_family = "Helvetica") +
#   theme(
#     panel.grid.minor = element_blank(),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     legend.position = "bottom",
#     legend.key.size = unit(0.5, "lines"),
#     legend.text = element_text(size = 6)
#   )
# ggsave(
#   "./plots/preQC_cumulative_distribution_doublet_score_by_individual.pdf",
#   device = cairo_pdf(),
#   width = 3, height = 3, useDingbats = FALSE
# )

# save removed cells with chosen QC cutoffs in a Seurat object
nFeature_lower <- 500
percent_mt_upper <- 1
doublet_scores_upper <- 0.2

data_obj_sub <- subset(data_obj,
  subset = nFeature_RNA > nFeature_lower & percent_mt < percent_mt_upper & doublet_scores <= doublet_scores_upper,
  invert = TRUE
)

meta_data_sub <- data_obj_sub@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

write_tsv(meta_data_sub, "./meta_QC_filtered_cells_post_cellBender.txt")
saveRDS(data_obj_sub, file = "./seurat_objects/snRNA_cellBender_corrected_QC_filteredCells_SeuratV4_object.rds")

# save all cells without any QC filtering in a Seurat object
saveRDS(data_obj, file = "./seurat_objects/snRNA_cellBender_corrected_noQCFilters_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
