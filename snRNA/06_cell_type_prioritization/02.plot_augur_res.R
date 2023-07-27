library(tidyverse)
library(ggrepel)
library(patchwork)

auc_df <- bind_rows(
  "ALS_vs_Control_MCX" = read_tsv("augur_out_obj_ALS_vs_Control_MCX_level2_nSubsamples_50_subsampleSize_200_folds_3.tsv"),
  "ALS_vs_Control_mFCX" = read_tsv("augur_out_obj_ALS_vs_Control_mFCX_level2_nSubsamples_50_subsampleSize_200_folds_3.tsv"),
  "FTD_vs_Control_MCX" = read_tsv("augur_out_obj_FTD_vs_Control_MCX_level2_nSubsamples_50_subsampleSize_200_folds_3.tsv"),
  "FTD_vs_Control_mFCX" = read_tsv("augur_out_obj_FTD_vs_Control_mFCX_level2_nSubsamples_50_subsampleSize_200_folds_3.tsv"),
  .id = "group"
)

auc_df_mod <- auc_df %>%
  pivot_wider(names_from = group, values_from = auc)

auc_df_ordered <- auc_df %>%
  mutate(
    disease = if_else(grepl("ALS_vs_", group), "ALS", "FTD"),
    region = if_else(grepl("_MCX", group), "Motor cortex", "Frontal Cortex"),
    region = factor(region, levels = c("Motor cortex", "Frontal Cortex"))
  ) %>%
  mutate(
    cell_type = factor(
      cell_type,
      levels = rev(c(
        "Astro", "Micro", "Oligo", "OPC", "VLMC",
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5"
      ))
    )
  ) %>%
  filter(
    !(disease == "FTD" & grepl("Exc_|Inh_", cell_type))
  )

# p <- auc_df_ordered  %>%
#   ggplot(aes(cell_type, auc))
# p + geom_point(stat = "identity") +
#   facet_wrap(~group) +
#   xlab("Cell type level 2") +
#   ylab("Augur AUC") +
#   coord_flip() +
#   theme_bw(base_size = 10, base_family = "Helvetica") +
#   theme(
#     panel.grid.minor = element_blank(),
#     strip.background = element_blank()
#   )
# ggsave("Augur_res_level2_nSubsamples_50_subsampleSize_200_folds_3.pdf",
#   device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE
# )

p + geom_point(stat = "identity") +
  facet_grid(disease ~ region, scales = "free", space = "free") +
  xlab("Cell types") +
  ylab("Augur AUC") +
  coord_flip() +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
ggsave("Augur_res_level2_nSubsamples_50_subsampleSize_200_folds_3_mod.pdf",
  device = cairo_pdf(), width = 4, height = 2.5, useDingbats = FALSE
)
dev.off()


# add shuffle results
shuffle_files <- list.files("./Augur_shuffle/", pattern = ".*ALS_vs_Control_MCX.*tsv", full.names = TRUE)
seed_list <- str_extract(shuffle_files, "seed_(\\d+).tsv", group = 1)

read_shuffle_files <- function(file) {
  read_tsv(file) %>%
    mutate(seed = str_extract(file, "seed_(\\d+).tsv", group = 1))
}

shuffle_auc_df <- map_dfr(shuffle_files, read_shuffle_files)

# define summary functions
sem <- function(x) {
  N <- sum(!is.na(x))
  sd(x, na.rm = T) / sqrt(N)
}

ci <- function(x, conf.interval = 0.95) {
  N <- sum(!is.na(x))
  sem <- sd(x, na.rm = T) / sqrt(N)
  qt(1 - (1 - conf.interval) / 2, N - 1) * sem
}

shuffle_auc_summary <- shuffle_auc_df %>%
  group_by(cell_type) %>%
  summarize(
    mean_auc = mean(auc),
    sd_auc = sd(auc),
    meidan_auc = median(auc),
    sem_auc = sem(auc),
    ci_auc = ci(auc)
  ) %>%
  rename(auc = mean_auc) %>%
  mutate(
    cell_type = factor(
      cell_type,
      levels = rev(c(
        "Astro", "Micro", "Oligo", "OPC", "VLMC",
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5"
      ))
    )
  )


p <- auc_df_ordered %>%
  filter(group == "ALS_vs_Control_MCX") %>%
  ggplot(aes(cell_type, auc))
p1 <- p +
  geom_point(aes(cell_type, auc), data = shuffle_auc_summary, color = "darkgrey") +
  # geom_errorbar(aes(cell_type, ymin = auc - sd_auc, ymax = auc + sd_auc),
  #   data = shuffle_auc_summary, width = 0.2, color = "darkgrey"
  # ) +
  # geom_errorbar(aes(cell_type, ymin = auc - 2 * sd_auc, ymax = auc + 2 * sd_auc),
  #   data = shuffle_auc_summary, width = 0.2, color = "darkgrey"
  # ) +
  geom_errorbar(aes(cell_type, ymin = auc - ci_auc, ymax = auc + ci_auc),
    data = shuffle_auc_summary, width = 0.2, color = "darkgrey"
  ) +
  geom_point(stat = "identity", color = "red") +
  xlab("Cell type") +
  ylab("Augur AUC") +
  coord_flip() +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
# ggsave("Augur_res_level2_nSubsamples_50_subsampleSize_200_folds_3_ALS_vs_Control_MCX_add_shuffles_mod.pdf",
#   device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE
# )
# ggsave("Augur_res_level2_nSubsamples_50_subsampleSize_200_folds_3_ALS_vs_Control_MCX_add_shuffles.pdf",
#   device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE
# )
# ggsave("Augur_res_level2_nSubsamples_50_subsampleSize_200_folds_3_ALS_vs_Control_MCX_add_shuffles_2SD.pdf",
#   device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE
# )
# dev.off()

# mFCX
shuffle_files <- list.files("./Augur_shuffle/", pattern = ".*ALS_vs_Control_mFCX.*tsv", full.names = TRUE)
seed_list <- str_extract(shuffle_files, "seed_(\\d+).tsv", group = 1)

shuffle_auc_df <- map_dfr(shuffle_files, read_shuffle_files)

shuffle_auc_summary <- shuffle_auc_df %>%
  filter(!cell_type %in% c("Endo", "Inh_ADARB2_Other")) %>%
  group_by(cell_type) %>%
  summarize(
    mean_auc = mean(auc),
    sd_auc = sd(auc),
    meidan_auc = median(auc),
    sem_auc = sem(auc),
    ci_auc = ci(auc)
  ) %>%
  rename(auc = mean_auc) %>%
  mutate(
    cell_type = factor(
      cell_type,
      levels = rev(c(
        "Astro", "Micro", "Oligo", "OPC", "VLMC",
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5"
      ))
    )
  )

p <- auc_df %>%
  filter(group == "ALS_vs_Control_mFCX") %>%
  filter(!cell_type %in% c("Endo", "Inh_ADARB2_Other", "VLMC")) %>%
  ggplot(aes(cell_type, auc))
p2 <- p +
  geom_point(aes(cell_type, auc), data = shuffle_auc_summary, color = "darkgrey") +
  # geom_errorbar(aes(cell_type, ymin = auc - sd_auc, ymax = auc + sd_auc),
  #   data = shuffle_auc_summary, width = 0.2, color = "darkgrey"
  # ) +
  # geom_errorbar(aes(cell_type, ymin = auc - 2 * sd_auc, ymax = auc + 2 * sd_auc),
  #   data = shuffle_auc_summary, width = 0.2, color = "darkgrey"
  # ) +
  geom_errorbar(aes(cell_type, ymin = auc - ci_auc, ymax = auc + ci_auc),
    data = shuffle_auc_summary, width = 0.2, color = "darkgrey"
  ) +
  geom_point(stat = "identity", color = "red") +
  xlab("Cell type") +
  ylab("Augur AUC") +
  coord_flip() +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
# ggsave("Augur_res_level2_nSubsamples_50_subsampleSize_200_folds_3_ALS_vs_Control_mFCX_add_shuffles.pdf",
#   device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE
# )
# ggsave("Augur_res_level2_nSubsamples_50_subsampleSize_200_folds_3_ALS_vs_Control_mFCX_add_shuffles_2SD.pdf",
#   device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE
# )
# dev.off()

# merge plots
p1 + p2
ggsave("Augur_res_level2_nSubsamples_50_subsampleSize_200_folds_3_ALS_vs_Control_add_shuffles_mod.pdf",
  device = cairo_pdf(), width = 4.8, height = 2.5, useDingbats = FALSE
)



# auc_df <- bind_rows(
#   "ALS_vs_Control_MCX" = read_tsv("augur_auc_ALS_vs_Control_MCX_level3.txt"),
#   "ALS_vs_Control_mFCX" = read_tsv("augur_auc_ALS_vs_Control_mFCX_level3.txt"),
#   "FTD_vs_Control_MCX" = read_tsv("augur_auc_FTD_vs_Control_MCX_level3.txt"),
#   "FTD_vs_Control_mFCX" = read_tsv("augur_auc_FTD_vs_Control_mFCX_level3.txt"),
#   .id = "group"
# )

# auc_df_mod <- auc_df %>%
#   pivot_wider(names_from = group, values_from = auc)

# p <- auc_df %>% ggplot(aes(cell_type, auc))
# p + geom_point(stat = "identity") +
#   facet_wrap(~group) +
#   xlab("Cell type level 3") +
#   ylab("Augur AUC") +
#   coord_flip() +
#   theme_bw(base_size = 8, base_family = "Helvetica") +
#   theme(panel.grid.minor = element_blank(),
#         strip.background = element_blank())
# ggsave("Augur_res_level3.pdf", device = cairo_pdf(), width = 8, height = 8, useDingbats = F)
#
# p <- auc_df_mod %>% ggplot(aes(ALS_MCX, ALS_mFCX))
# p +
#   geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
#   geom_point() +
#   geom_text_repel(aes(label = cell_type)) +
#   xlab("Augur AUC, ALS vs. Control, MCX") +
#   ylab("Augur AUC, ALS vs. Control, mFCX") +
#   xlim(0.49, 0.72) +
#   ylim(0.49, 0.72) +
#   coord_equal() +
#   theme_bw(base_size = 10, base_family = "Helvetica")
# ggsave("auguc_auc_ALS_Control_MCX_vs_mFCX.pdf", device = cairo_pdf(),
#        width = 4, height = 4, useDingbats = FALSE)
#
# p <- auc_df_mod %>% ggplot(aes(FTD_MCX, FTD_mFCX))
# p +
#   geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
#   geom_point() +
#   geom_text_repel(aes(label = cell_type)) +
#   xlab("Augur AUC, FTD vs. Control, MCX") +
#   ylab("Augur AUC, FTD vs. Control, mFCX") +
#   xlim(0.49, 0.72) +
#   ylim(0.49, 0.72) +
#   coord_equal() +
#   theme_bw(base_size = 10, base_family = "Helvetica")
# ggsave("auguc_auc_FTD_Control_MCX_vs_mFCX.pdf", device = cairo_pdf(),
#        width = 4, height = 4, useDingbats = FALSE)
#
# p <- auc_df_mod %>% ggplot(aes(ALS_MCX, FTD_MCX))
# p +
#   geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
#   geom_point() +
#   geom_text_repel(aes(label = cell_type)) +
#   xlab("Augur AUC, ALS vs. Control, MCX") +
#   ylab("Augur AUC, FTD vs. Control, MCX") +
#   xlim(0.49, 0.72) +
#   ylim(0.49, 0.72) +
#   coord_equal() +
#   theme_bw(base_size = 10, base_family = "Helvetica")
# ggsave("auguc_auc_MCX_ALS_Control_vs_FTD_Control.pdf", device = cairo_pdf(),
#        width = 4, height = 4, useDingbats = FALSE)
#
# p <- auc_df_mod %>% ggplot(aes(ALS_mFCX, FTD_mFCX))
# p +
#   geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
#   geom_point() +
#   geom_text_repel(aes(label = cell_type)) +
#   xlab("Augur AUC, ALS vs. Control, mFCX") +
#   ylab("Augur AUC, FTD vs. Control, mFCX") +
#   xlim(0.49, 0.72) +
#   ylim(0.49, 0.72) +
#   coord_equal() +
#   theme_bw(base_size = 10, base_family = "Helvetica")
# ggsave("auguc_auc_mFCX_ALS_Control_vs_FTD_Control.pdf", device = cairo_pdf(),
#        width = 4, height = 4, useDingbats = FALSE)
#
# ## read cell counts
# meta_data <- read_tsv("../all_snRNA_dataset_merged_Seurat_object_subset_nFeature_500to8000_percentMT_1_scTransform_normalized_clustered_annotated_cell_metadata_annotatedCluster.txt",
#                       col_types = "cciicccdiiiicdic")
#
# atac_metadata <- read_tsv("/cndd2/junhao/ALS_FTD_singleCell/integration/data/organized/atac_metadata.tsv")
#
# rna_metadata_sub <- meta_data %>%
#   select(cell, subject, disease, region, seurat_clusters) %>%
#   mutate(subject = as.character(subject)) %>%
#   rename(cluster = seurat_clusters,
#          cell_id = cell)
#
# atac_metadata_sub <- atac_metadata %>%
#   mutate(sample_cell_id = paste(barcode, sample, sep = "_")) %>%
#   mutate(subject = str_split_fixed(sample, "_", 4)[,2]) %>%
#   select(sample_cell_id, subject, disease, brain_region, cluster) %>%
#   rename(cell_id = sample_cell_id,
#          region = brain_region)
#
# metadata_merge <- bind_rows(rna_metadata_sub, atac_metadata_sub) %>%
#   mutate(disease = str_to_upper(disease),
#          region = if_else(region == "FCX", "mFCX", region))
#
# integrate_clusters <- read_tsv("/cndd2/junhao/ALS_FTD_singleCell/integration/results/intg_summary_alsftd_intg_cleanRNA.tsv") %>%
#   rename(cell_id = sample) %>%
#   left_join(metadata_merge)
#
# rna_metadata_plus_joint_clusters <- meta_data %>%
#   left_join(integrate_clusters, by = c("cell" = "cell_id"))
#
# rna_metadata_plus_joint_clusters_mod <- rna_metadata_plus_joint_clusters %>%
#   mutate(cluster_joint_high = case_when(
#     cluster_joint_r1 %in% c(9, 14, 20, 22) ~ "Astrocyte",
#     cluster_joint_r1 %in% c(12, 23) ~ "Microglia",
#     cluster_joint_r1 %in% c(1, 2, 3, 7, 10, 16, 15) ~ "Oligodendrocyte",
#     cluster_joint_r1 %in% c(6) ~ "OPC",
#     cluster_joint_r1 == 21 & cluster_joint_r2 == 45 ~ "Endothelial",
#     cluster_joint_r1 == 21 & cluster_joint_r2 == 35 ~ "Endothelial",
#     cluster_joint_r1 %in% c(5, 8, 13, 17, 24) ~ "Exc-neurons",
#     cluster_joint_r1 %in% c(4, 11, 19) ~ "Inh-neurons",
#     TRUE ~ "Mixed"
#   ))
#
# meta_data_mod <- rna_metadata_plus_joint_clusters_mod %>%
#   select(-subject.y, -disease.y, -region.y) %>%
#   rename(subject = subject.x, disease = disease.x, region = region.x)
#
# cell_counts <- meta_data_mod %>%
#   count(region, disease, cluster_joint_high)
#
# cell_counts_1 <- cell_counts %>%
#   filter(cluster_joint_high != "Mixed",
#          disease != "FTD") %>%
#   pivot_wider(names_from = "disease", values_from = "n") %>%
#   mutate(n_sum = ALS + Control,
#          n_mean = (ALS + Control)/2,
#          disease = "ALS",
#          group = paste(disease, region, sep = "_")) %>%
#   rename(cell_type = cluster_joint_high) %>%
#   select(group, cell_type, n_sum, n_mean)
#
# cell_counts_2 <- cell_counts %>%
#   filter(cluster_joint_high != "Mixed",
#          disease != "ALS") %>%
#   pivot_wider(names_from = "disease", values_from = "n") %>%
#   mutate(n_sum = FTD + Control,
#          n_mean = (FTD + Control)/2,
#          disease = "FTD",
#          group = paste(disease, region, sep = "_")) %>%
#   rename(cell_type = cluster_joint_high) %>%
#   select(group, cell_type, n_sum, n_mean)
#
# cell_counts_mod <- bind_rows(cell_counts_1, cell_counts_2)
#
# ## read DE genes count
# rna_DE <- read_tsv("/cndd2/junhao/ALS_FTD_singleCell/snRNA/MAST_DE_genes_disease_by_cell_type_by_region_integrated_clusters_highLevel_sig.txt") %>%
#   mutate(cell_type = str_replace(cell_type, "Neu", "-neu"))
#
# rna_DE_count <- rna_DE %>%
#   filter(cell_type != "Mixed", grepl("vs_Control", comparison)) %>%
#   count(comparison, region, cell_type) %>%
#   mutate(group = paste(str_replace(comparison, "_vs_Control", ""), region, sep = "_")) %>%
#   rename(n_DE = n) %>%
#   select(group, cell_type, n_DE)
#
# merged_df <- cell_counts_mod %>%
#   left_join(auc_df) %>%
#   left_join(rna_DE_count)
#
# p <- merged_df %>% ggplot(aes(n_mean, n_DE))
# p + geom_point() +
#   geom_text_repel(aes(label = cell_type), size = 3) +
#   facet_wrap(~group) +
#   scale_x_log10() +
#   scale_y_log10() +
#   annotation_logticks() +
#   xlab("Number of cells\n(average across disease and control)") +
#   ylab("Number of DE genes\n(compared to control)") +
#   theme_bw(base_size = 8, base_family = "Helvetica") +
#   theme(strip.background = element_blank())
# ggsave("num_cells_vs_num_DE_genes.pdf", device = cairo_pdf(),
#        width = 6, height = 6, useDingbats = F)
#
# p <- merged_df %>% ggplot(aes(n_mean, auc))
# p + geom_point() +
#   geom_text_repel(aes(label = cell_type), size = 3) +
#   facet_wrap(~group) +
#   scale_x_log10() +
#   annotation_logticks(sides = "b") +
#   xlab("Number of cells\n(average across disease and control)") +
#   ylab("Augur AUC") +
#   theme_bw(base_size = 8, base_family = "Helvetica") +
#   theme(strip.background = element_blank())
# ggsave("num_cells_vs_augur_auc.pdf", device = cairo_pdf(),
#        width = 6, height = 6, useDingbats = F)
#
# p <- merged_df %>% ggplot(aes(n_DE, auc))
# p + geom_point() +
#   geom_text_repel(aes(label = cell_type), size = 3) +
#   facet_wrap(~group) +
#   scale_x_log10() +
#   annotation_logticks(sides = "b") +
#   xlab("Number of DE genes\n(compared to control)") +
#   ylab("Augur AUC") +
#   theme_bw(base_size = 8, base_family = "Helvetica") +
#   theme(strip.background = element_blank())
# ggsave("num_DE_genes_vs_augur_auc.pdf", device = cairo_pdf(),
#        width = 6, height = 6, useDingbats = F)
#
# library(ggpubr)
# # sp <- ggscatter(merged_df, x = "n_mean", y = "n_DE",
# # sp <- ggscatter(merged_df, x = "n_mean", y = "auc",
# sp <- ggscatter(merged_df, x = "n_DE", y = "auc",
#                 color = "group", palette = "jco",
#                 add = "reg.line", conf.int = TRUE)
# sp + stat_cor(aes(color = group), label.x = 3)  + facet_wrap(~group)
