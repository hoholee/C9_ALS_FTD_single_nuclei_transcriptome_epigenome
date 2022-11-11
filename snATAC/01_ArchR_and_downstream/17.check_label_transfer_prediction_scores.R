# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)

metadata <- read_tsv("metadata_merged.txt")

# full dataset
p <- metadata %>% ggplot(aes(predictedScore_from_snRNA_full_level2))
p + geom_histogram() +
    facet_wrap(~Clusters_full, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_full_dataset_clusters_vs_full_level_2_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~Clusters_full) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_full_dataset_clusters_vs_full_level_2_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata %>% ggplot(aes(predictedScore_from_snRNA_full_level3))
p + geom_histogram() +
    facet_wrap(~Clusters_full, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_full_dataset_clusters_vs_full_level_3_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~Clusters_full) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_full_dataset_clusters_vs_full_level_3_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata %>% ggplot(aes(predictedScore_from_snRNA_full_level3))
p + geom_histogram() +
    facet_wrap(~predictedGroup_from_snRNA_full_level3, scales = "free_y") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_full_dataset_predicted_labels_vs_full_level_3_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~predictedGroup_from_snRNA_full_level3) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_full_dataset_predicted_labels_vs_full_level_3_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata %>% ggplot(aes(predictedScore_from_snRNA_full_level2))
p + geom_histogram() +
    facet_wrap(~predictedGroup_from_snRNA_full_level2, scales = "free_y") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_full_dataset_predicted_labels_vs_full_level_2_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~predictedGroup_from_snRNA_full_level2) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_full_dataset_predicted_labels_vs_full_level_2_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)


# sub-clustering
metadata_Exc <- metadata %>% filter(!is.na(Clusters_sub_Exc))
metadata_Inh <- metadata %>% filter(!is.na(Clusters_sub_Inh))
metadata_Glia <- metadata %>% filter(!is.na(Clusters_sub_Glia))

# Exc
p <- metadata_Exc %>% ggplot(aes(predictedScore_from_snRNA_Exc_level2))
p + geom_histogram() +
    facet_wrap(~Clusters_sub_Exc, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Exc_subset_clusters_vs_full_level_2_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Exc %>% ggplot(aes(predictedScore_from_snRNA_Exc_level2))
p + stat_ecdf(pad = FALSE) +
    facet_wrap(~Clusters_sub_Exc, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Exc_subset_clusters_vs_full_level_2_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Exc %>% ggplot(aes(predictedScore_from_snRNA_Exc_level3))
p + geom_histogram() +
    facet_wrap(~Clusters_sub_Exc, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Exc_subset_clusters_vs_full_level_3_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Exc %>% ggplot(aes(predictedScore_from_snRNA_Exc_level3))
p + stat_ecdf(pad = FALSE) +
    facet_wrap(~Clusters_sub_Exc, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Exc_subset_clusters_vs_full_level_3_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Exc %>% ggplot(aes(predictedScore_from_snRNA_Exc_level2))
p + geom_histogram() +
    facet_wrap(~predictedGroup_from_snRNA_Exc_level2, scales = "free_y") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Exc_subset_predicted_labels_vs_full_level_2_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~predictedGroup_from_snRNA_Exc_level2) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Exc_subset_predicted_labels_vs_full_level_2_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Exc %>% ggplot(aes(predictedScore_from_snRNA_Exc_level3))
p + geom_histogram() +
    facet_wrap(~predictedGroup_from_snRNA_Exc_level3, scales = "free_y") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Exc_subset_predicted_labels_vs_full_level_3_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~predictedGroup_from_snRNA_Exc_level3) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Exc_subset_predicted_labels_vs_full_level_3_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

# p <- metadata_Exc %>%
#     filter(predictedScore_from_snRNA_Exc_level2 >= 0.8) %>%
#     ggplot(aes(Exc_UMAP_1, Exc_UMAP_2))
# p + geom_point(aes(color = predictedScore_from_snRNA_Exc_level2))

# p <- metadata_Exc %>%
#     # filter(predictedScore_from_snRNA_Exc_level3 >= 0.5) %>%
#     ggplot(aes(Exc_UMAP_1, Exc_UMAP_2))
# p + geom_point(aes(color = predictedScore_from_snRNA_Exc_level3))

# Inh
p <- metadata_Inh %>% ggplot(aes(predictedScore_from_snRNA_Inh_level2))
p + geom_histogram() +
    facet_wrap(~Clusters_sub_Inh, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Inh_subset_clusters_vs_full_level_2_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Inh %>% ggplot(aes(predictedScore_from_snRNA_Inh_level2))
p + stat_ecdf(pad = FALSE) +
    facet_wrap(~Clusters_sub_Inh, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Inh_subset_clusters_vs_full_level_2_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Inh %>% ggplot(aes(predictedScore_from_snRNA_Inh_level3))
p + geom_histogram() +
    facet_wrap(~Clusters_sub_Inh, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Inh_subset_clusters_vs_full_level_3_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Inh %>% ggplot(aes(predictedScore_from_snRNA_Inh_level3))
p + stat_ecdf(pad = FALSE) +
    facet_wrap(~Clusters_sub_Inh, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Inh_subset_clusters_vs_full_level_3_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Inh %>% ggplot(aes(predictedScore_from_snRNA_Inh_level2))
p + geom_histogram() +
    facet_wrap(~predictedGroup_from_snRNA_Inh_level2, scales = "free_y") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Inh_subset_predicted_labels_vs_full_level_2_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~predictedGroup_from_snRNA_Inh_level2) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Inh_subset_predicted_labels_vs_full_level_2_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Inh %>% ggplot(aes(predictedScore_from_snRNA_Inh_level3))
p + geom_histogram() +
    facet_wrap(~predictedGroup_from_snRNA_Inh_level3, scales = "free_y") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Inh_subset_predicted_labels_vs_full_level_3_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~predictedGroup_from_snRNA_Inh_level3) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Inh_subset_predicted_labels_vs_full_level_3_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)
# p <- metadata_Inh %>% ggplot(aes(predictedScore_from_snRNA_Inh_level3))
# p + geom_histogram() +
#     facet_wrap(~predictedGroup_from_snRNA_Inh_level3, scales = "free_y")

# p <- metadata_Inh %>%
#     # filter(predictedScore_from_snRNA_Inh_level3 >= 0.5) %>%
#     ggplot(aes(Inh_UMAP_1, Inh_UMAP_2))
# p + geom_point(aes(color = predictedScore_from_snRNA_Inh_level3))


### Glia
p <- metadata_Glia %>% ggplot(aes(predictedScore_from_snRNA_Glia_level2))
p + geom_histogram() +
    facet_wrap(~Clusters_sub_Glia, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Glia_subset_clusters_vs_full_level_2_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Glia %>% ggplot(aes(predictedScore_from_snRNA_Glia_level2))
p + stat_ecdf(pad = FALSE) +
    facet_wrap(~Clusters_sub_Glia, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Glia_subset_clusters_vs_full_level_2_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Glia %>% ggplot(aes(predictedScore_from_snRNA_Glia_level3))
p + geom_histogram() +
    facet_wrap(~Clusters_sub_Glia, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Glia_subset_clusters_vs_full_level_3_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Glia %>% ggplot(aes(predictedScore_from_snRNA_Glia_level3))
p + stat_ecdf(pad = FALSE) +
    facet_wrap(~Clusters_sub_Glia, scales = "free_y") +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Glia_subset_clusters_vs_full_level_3_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Glia %>% ggplot(aes(predictedScore_from_snRNA_Glia_level2))
p + geom_histogram() +
    facet_wrap(~predictedGroup_from_snRNA_Glia_level2, scales = "free_y") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Glia_subset_predicted_labels_vs_full_level_2_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~predictedGroup_from_snRNA_Glia_level2) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Glia_subset_predicted_labels_vs_full_level_2_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p <- metadata_Glia %>% ggplot(aes(predictedScore_from_snRNA_Glia_level3))
p + geom_histogram() +
    facet_wrap(~predictedGroup_from_snRNA_Glia_level3, scales = "free_y") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Glia_subset_predicted_labels_vs_full_level_3_hist.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)

p + stat_ecdf(pad = FALSE) +
    facet_wrap(~predictedGroup_from_snRNA_Glia_level3) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/snRNA_label_transfer_prediction_score_distribution_Glia_subset_predicted_labels_vs_full_level_3_ecdf.pdf",
    device = cairo_pdf(),
    width = 8, height = 6, useDingbats = FALSE
)
# p <- metadata_Glia %>%
#     # filter(predictedScore_from_snRNA_Glia_level3 >= 0.5) %>%
#     ggplot(aes(Glia_UMAP_1, Glia_UMAP_2))
# p + geom_point(aes(color = predictedScore_from_snRNA_Glia_level3))


# define final annotations based on confusion matrix
# first only keep cells with prediciton scores > `prediciton_score_cutoff`
# then get confusion matrix of snATAC clusters vs. snRNA transferred labels (level 2 and 3)
# take major votes for each cluster
# require votes percentage > `percent_cutoff`
# otherwise label the cluster as "Mixed", appended with the major vote and the original cluster ID
predition_score_cutoff <- 0.7
percent_cutoff <- 0.7

# Glia
metadata_Glia_level2_high_confidence_sub <- metadata_Glia %>%
    filter(predictedScore_from_snRNA_Glia_level2 > predition_score_cutoff)
metadata_Glia_level3_high_confidence_sub <- metadata_Glia %>%
    filter(predictedScore_from_snRNA_Glia_level3 > predition_score_cutoff)

mat_Glia_level2 <- confusionMatrix(
    metadata_Glia_level2_high_confidence_sub$Clusters_sub_Glia,
    metadata_Glia_level2_high_confidence_sub$predictedGroup_from_snRNA_Glia_level2
)

mat_Glia_level3 <- confusionMatrix(
    metadata_Glia_level3_high_confidence_sub$Clusters_sub_Glia,
    metadata_Glia_level3_high_confidence_sub$predictedGroup_from_snRNA_Glia_level3
)

percent_votes_Glia_level2 <- rowMaxs(mat_Glia_level2) / rowSums(mat_Glia_level2)
percent_votes_Glia_level3 <- rowMaxs(mat_Glia_level3) / rowSums(mat_Glia_level3)

Glia_level2_label_old <- rownames(mat_Glia_level2)
Glia_level2_label_new <- colnames(mat_Glia_level2)[apply(mat_Glia_level2, 1, which.max)]
idx <- which(!percent_votes_Glia_level2 > percent_cutoff)
Glia_level2_label_new[idx] <- paste("Glia_level2_Mixed", Glia_level2_label_new[idx], Glia_level2_label_old[idx], sep = "__")
names(Glia_level2_label_new) <- Glia_level2_label_old

Glia_level3_label_old <- rownames(mat_Glia_level3)
Glia_level3_label_new <- colnames(mat_Glia_level3)[apply(mat_Glia_level3, 1, which.max)]
idx <- which(!percent_votes_Glia_level3 > percent_cutoff)
Glia_level3_label_new[idx] <- paste("Glia_level3_Mixed", Glia_level3_label_new[idx], Glia_level3_label_old[idx], sep = "__")
names(Glia_level3_label_new) <- Glia_level3_label_old

Glia_annotation <- tibble(
    Clusters_sub_Glia = Glia_level2_label_old,
    atac_anno_level_1 = "Glia",
    atac_anno_level_2 = Glia_level2_label_new[Glia_level2_label_old],
    atac_anno_level_3 = Glia_level3_label_new[Glia_level2_label_old]
)

# p <- metadata_Glia_level2_high_confidence_sub %>% ggplot(aes(Glia_UMAP_1, Glia_UMAP_2))
# p + geom_point(aes(color = predictedGroup_from_snRNA_Glia_level2))

# p <- metadata_Glia_level3_high_confidence_sub %>% ggplot(aes(Glia_UMAP_1, Glia_UMAP_2))
# p + geom_point(aes(color = predictedGroup_from_snRNA_Glia_level3))

# Exc
metadata_Exc_level2_high_confidence_sub <- metadata_Exc %>%
    filter(predictedScore_from_snRNA_Exc_level2 > predition_score_cutoff)
metadata_Exc_level3_high_confidence_sub <- metadata_Exc %>%
    filter(predictedScore_from_snRNA_Exc_level3 > predition_score_cutoff)

mat_Exc_level2 <- confusionMatrix(
    metadata_Exc_level2_high_confidence_sub$Clusters_sub_Exc,
    metadata_Exc_level2_high_confidence_sub$predictedGroup_from_snRNA_Exc_level2
)

mat_Exc_level3 <- confusionMatrix(
    metadata_Exc_level3_high_confidence_sub$Clusters_sub_Exc,
    metadata_Exc_level3_high_confidence_sub$predictedGroup_from_snRNA_Exc_level3
)

percent_votes_Exc_level2 <- rowMaxs(mat_Exc_level2) / rowSums(mat_Exc_level2)
percent_votes_Exc_level3 <- rowMaxs(mat_Exc_level3) / rowSums(mat_Exc_level3)

Exc_level2_label_old <- rownames(mat_Exc_level2)
Exc_level2_label_new <- colnames(mat_Exc_level2)[apply(mat_Exc_level2, 1, which.max)]
idx <- which(!percent_votes_Exc_level2 > percent_cutoff)
Exc_level2_label_new[idx] <- paste("Exc_level2_Mixed", Exc_level2_label_new[idx], Exc_level2_label_old[idx], sep = "__")
names(Exc_level2_label_new) <- Exc_level2_label_old

Exc_level3_label_old <- rownames(mat_Exc_level3)
Exc_level3_label_new <- colnames(mat_Exc_level3)[apply(mat_Exc_level3, 1, which.max)]
idx <- which(!percent_votes_Exc_level3 > percent_cutoff)
Exc_level3_label_new[idx] <- paste("Exc_level3_Mixed", Exc_level3_label_new[idx], Exc_level3_label_old[idx], sep = "__")
names(Exc_level3_label_new) <- Exc_level3_label_old

Exc_annotation <- tibble(
    Clusters_sub_Exc = Exc_level2_label_old,
    atac_anno_level_1 = "Exc",
    atac_anno_level_2 = Exc_level2_label_new[Exc_level2_label_old],
    atac_anno_level_3 = Exc_level3_label_new[Exc_level2_label_old]
)

# Inh
metadata_Inh_level2_high_confidence_sub <- metadata_Inh %>%
    filter(predictedScore_from_snRNA_Inh_level2 > predition_score_cutoff)
metadata_Inh_level3_high_confidence_sub <- metadata_Inh %>%
    filter(predictedScore_from_snRNA_Inh_level3 > predition_score_cutoff)

mat_Inh_level2 <- confusionMatrix(
    metadata_Inh_level2_high_confidence_sub$Clusters_sub_Inh,
    metadata_Inh_level2_high_confidence_sub$predictedGroup_from_snRNA_Inh_level2
)

mat_Inh_level3 <- confusionMatrix(
    metadata_Inh_level3_high_confidence_sub$Clusters_sub_Inh,
    metadata_Inh_level3_high_confidence_sub$predictedGroup_from_snRNA_Inh_level3
)

percent_votes_Inh_level2 <- rowMaxs(mat_Inh_level2) / rowSums(mat_Inh_level2)
percent_votes_Inh_level3 <- rowMaxs(mat_Inh_level3) / rowSums(mat_Inh_level3)

Inh_level2_label_old <- rownames(mat_Inh_level2)
Inh_level2_label_new <- colnames(mat_Inh_level2)[apply(mat_Inh_level2, 1, which.max)]
idx <- which(!percent_votes_Inh_level2 > percent_cutoff)
Inh_level2_label_new[idx] <- paste("Inh_level2_Mixed", Inh_level2_label_new[idx], Inh_level2_label_old[idx], sep = "__")
names(Inh_level2_label_new) <- Inh_level2_label_old

Inh_level3_label_old <- rownames(mat_Inh_level3)
Inh_level3_label_new <- colnames(mat_Inh_level3)[apply(mat_Inh_level3, 1, which.max)]
idx <- which(!percent_votes_Inh_level3 > percent_cutoff)
Inh_level3_label_new[idx] <- paste("Inh_level3_Mixed", Inh_level3_label_new[idx], Inh_level3_label_old[idx], sep = "__")
names(Inh_level3_label_new) <- Inh_level3_label_old

Inh_annotation <- tibble(
    Clusters_sub_Inh = Inh_level2_label_old,
    atac_anno_level_1 = "Inh",
    atac_anno_level_2 = Inh_level2_label_new[Inh_level2_label_old],
    atac_anno_level_3 = Inh_level3_label_new[Inh_level2_label_old]
)

final_metadata_Glia <- metadata_Glia %>%
    left_join(Glia_annotation, by = "Clusters_sub_Glia") %>%
    select(cell_id, atac_anno_level_1, atac_anno_level_2, atac_anno_level_3)

final_metadata_Exc <- metadata_Exc %>%
    left_join(Exc_annotation, by = "Clusters_sub_Exc") %>%
    select(cell_id, atac_anno_level_1, atac_anno_level_2, atac_anno_level_3)

final_metadata_Inh <- metadata_Inh %>%
    left_join(Inh_annotation, by = "Clusters_sub_Inh") %>%
    select(cell_id, atac_anno_level_1, atac_anno_level_2, atac_anno_level_3)

final_anno <- bind_rows(final_metadata_Glia, final_metadata_Exc, final_metadata_Inh)

final_metadata <- metadata %>% left_join(final_anno, by = "cell_id")
write_tsv(final_metadata, "./metadata_merged_addAnno.txt")