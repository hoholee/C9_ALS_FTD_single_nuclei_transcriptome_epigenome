# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)
library(harmony)
library(ggrastr)
library(scCustomize)
library(shades)
library(ggbeeswarm)
library(ggpubr)

## read the sctransform-normalized seurat object
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBender_corrected_noQCFilters_scTransformed_clustered_SeuratV4_object.rds")

meta_data_full <- data_obj_sub@meta.data %>%
       rownames_to_column("cell_id") %>%
       as_tibble()

## read meta data from the filtered cells (keep "NK_cell", "Exc_unknown", "Ambiguous")
meta_data <- read_tsv("./metadata_all_cells_2nd_round_annotations.txt")

annotations <- meta_data %>%
       select(cell_id, rna_anno_2ndRound_level_1, rna_anno_2ndRound_level_2, rna_anno_2ndRound_level_3)

res <- meta_data_full %>%
       left_join(annotations, by = "cell_id") %>%
       mutate(
              rna_anno_2ndRound_level_1 = replace_na(rna_anno_2ndRound_level_1, "filtered"),
              rna_anno_2ndRound_level_2 = replace_na(rna_anno_2ndRound_level_2, "filtered"),
              rna_anno_2ndRound_level_3 = replace_na(rna_anno_2ndRound_level_3, "filtered")
       )

data_obj_sub$rna_anno_2ndRound_level_1 <- res$rna_anno_2ndRound_level_1
data_obj_sub$rna_anno_2ndRound_level_2 <- res$rna_anno_2ndRound_level_2
data_obj_sub$rna_anno_2ndRound_level_3 <- res$rna_anno_2ndRound_level_3

DimPlot_scCustom(
       seurat_object = data_obj_sub,
       # split.by = "rna_anno_2ndRound_level_2",
       split.by = "disease",
       # group.by = "rna_anno_2ndRound_level_1",
       raster = TRUE
)

FeaturePlot_scCustom(
       seurat_object = data_obj_sub,
       features = "percent_mt",
       split.by = "disease",
       raster = TRUE
)

FeaturePlot_scCustom(
       seurat_object = data_obj_sub,
       features = "percent_mt",
       split.by = "rna_anno_2ndRound_level_1",
       raster = TRUE
)

level_1_markers <- c("SNAP25", "MEF2C", "SLC17A7", "SATB2", "GAD1", "GAD2")

level_2_glia_markers <- c(
       "AQP4", "APOE", "GFAP", "HPSE2", "CABLES1", "CD44", "WDR49", # Astro
       "PECAM1", "VWF", "CLDN5", "NOSTRIN", # Endo
       "C1QB", "CSF1R", # Micro
       "MOG", "MBP", "MOBP", "CNP", "ENPP6", "OPALIN", # Oligo
       "PDGFRA", "VCAN", # OPC
       "TBX18", "COLEC12", "CEMIP", "CYP1B1", "P2RY14" # VLMC
)

level_2_inh_markers <- c(
       "GAD1", "GAD2",
       "LHX6", "SOX6",
       "PVALB", "CUX2", "MYBPC1", "PIEZO2", "SCUBE3",
       "SST", "NPY", "GPC5", "EDNRA", "CDH12", "KLHL14",
       "ADARB2",
       "LAMP5", "SV2C", "CHST9", "CPLX3", "KIT", "CXCL14", "NDNF", "RELN", "PAX6", "SEMA3C", "LINC01470", "SCML4",
       "VIP", "ABI3BP", "CLSTN2", "DACH2", "FLT1", "ZBTB20"
)

level_2_exc_markers <- c(
       "SLC17A7", "SATB2",
       "CUX2",
       "ACVR1C",
       "LAMP5", "SERPINE2", "LINC00507",
       "PDGFD", "CCBE1",
       "LRRC2", "GLIS3",
       "SV2C", "SYT2",
       "RORB",
       "COBLL1", "PLCH1", "PRSS12", "CCDC68",
       "OTOGL", "COL22A1", "ALDH1A1", "GABRG1", "NPFFR2",
       "LRRK1", "TRABD2A", "GRIN3A", "ADAMTSL1", "RPRM",
       "THEMIS", "SMYD1",
       "LINC00299", "CFLAR",
       "TLE4", "PCSK5", "MDFIC", "KCNK2",
       "SULF1", "SEMA5A", "FEZF2", "HTR2C", "KCNIP1",
       "CRYM", "ADRA1A", "MYO16", "NEFH", "POU3F1"
)

DotPlot_scCustom(
       seurat_object = data_obj_sub,
       # features = level_1_markers,
       # features = level_2_glia_markers,
       # features = level_2_inh_markers,
       features = level_2_exc_markers,
       flip_axes = TRUE
       # split.by = "rna_anno_2ndRound_level_1",
)

# annotate clusters with previous labels
res_annotated <- res %>%
       mutate(
              rna_anno_level1_with_filtered_cells = case_when(
                     seurat_clusters %in% c(2, 14, 33, 7, 28, 9, 15, 17, 27, 34, 8, 20) ~ "Exc_neuron",
                     seurat_clusters %in% c(25, 18, 19, 13, 23, 29) ~ "Inh_neuron",
                     TRUE ~ "Non_neuron"
              )
       )

# group by cell class, separated by individuals
# filter by QC metrics cutoff
num_genes_cutoff <- 500
doublet_score_cutoff <- 0.2
# percent_mt_cutoff <- 1
# percent_mt_cutoff <- 5
percent_mt_cutoff <- 10
# no filters
# num_genes_cutoff <- 0
# doublet_score_cutoff <- 1.1
# percent_mt_cutoff <- 101

cell_count_class <- res_annotated %>%
       filter(
              nFeature_RNA > num_genes_cutoff,
              doublet_scores < doublet_score_cutoff,
              percent_mt < percent_mt_cutoff
       ) %>%
       count(orig.ident, region, disease, subject, rna_anno_level1_with_filtered_cells) %>%
       mutate(rna_anno_level1_with_filtered_cells = factor(rna_anno_level1_with_filtered_cells,
              levels = c("Exc_neuron", "Inh_neuron", "Non_neuron")
       ))

cell_class_color <- c(
       "Exc_neuron" = "#96BB45",
       "Inh_neuron" = "#718DC7",
       "Non_neuron" = "#9C482B"
)

#
# p <- cell_count_class %>% ggplot(aes(subject, n))

# p + geom_bar(aes(fill = rna_anno_level1_with_filtered_cells),
#        position = position_fill(),
#        stat = "identity"
# ) +
#        facet_wrap(~ disease + region,
#               dir = "v",
#               drop = TRUE,
#               scales = "free_x",
#               nrow = 2
#        ) +
#        scale_fill_manual(values = cell_class_color, name = "Cell class") +
#        xlab("Individual IDs") +
#        ylab("Proportion of types (%)") +
#        theme_bw(base_size = 8, base_family = "Helvetica") +
#        theme(
#               strip.background = element_blank(),
#               panel.grid.minor = element_blank()
#        )

# statistical test of the neuronal cells proportion differences in FTD using individuals
cell_proportion <- cell_count_class %>%
       select(-orig.ident) %>%
       pivot_wider(names_from = rna_anno_level1_with_filtered_cells, values_from = n) %>%
       mutate(
              total_count = Exc_neuron + Inh_neuron + Non_neuron,
              neuron_proportion = (Exc_neuron + Inh_neuron) / total_count,
              exc_proportion = Exc_neuron / (Exc_neuron + Inh_neuron),
              disease = factor(disease, levels = c("ALS", "FTD", "Control")),
              region = factor(region, levels = c("MCX", "mFCX"))
       )

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")
disease_color_dark <- lightness(disease_color, scalefac(0.8)) %>%
       as.vector()
names(disease_color_dark) <- names(disease_color)

p <- cell_proportion %>% ggplot(aes(disease, 100 * neuron_proportion))
p +
       geom_boxplot(aes(color = disease), width = 0.5, outlier.shape = NA) +
       geom_beeswarm(aes(color = disease), size = 1, cex = 4) +
       facet_wrap(~region, nrow = 1, strip.position = "bottom") +
       stat_compare_means(
              aes(label = ..p.signif..),
              method = "t.test",
              comparisons = list(
                     c("ALS", "Control"),
                     c("FTD", "Control"),
                     c("ALS", "FTD")
              )
       ) +
       ggtitle(str_glue("Num.genes > {num_genes_cutoff}, Doublet score < {doublet_score_cutoff}, %Mt < {percent_mt_cutoff}")) +
       xlab("Diagnosis") +
       ylab("Proportion of neurons (%)") +
       coord_cartesian(ylim = c(0, 100)) +
       scale_color_manual(values = disease_color_dark, name = "diagnosis") +
       theme_bw(base_size = 8, base_family = "Helvetica") +
       theme(
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              strip.placement = "outside",
              legend.position = "none"
       )

ggsave(str_glue("./plots/compare_neurons_proportion_between_disease_diagnosis_fullDatasetClustering_NumGenes_{num_genes_cutoff}_DoubletScore_{doublet_score_cutoff}_PercentMt_{percent_mt_cutoff}.pdf"),
       device = cairo_pdf(), width = 2, height = 2, useDingbats = FALSE
)

p <- cell_proportion %>% ggplot(aes(disease, 100 * exc_proportion))
p +
       geom_boxplot(aes(color = disease), width = 0.5, outlier.shape = NA) +
       geom_beeswarm(aes(color = disease), size = 1, cex = 4) +
       facet_wrap(~region, nrow = 1, strip.position = "bottom") +
       stat_compare_means(
              aes(label = ..p.signif..),
              method = "t.test",
              comparisons = list(
                     c("ALS", "Control"),
                     c("FTD", "Control"),
                     c("ALS", "FTD")
              )
       ) +
       ggtitle(str_glue("Num.genes > {num_genes_cutoff}, Doublet score < {doublet_score_cutoff}, %Mt < {percent_mt_cutoff}")) +
       xlab("Diagnosis") +
       ylab("Proportion of excitatory neurons in all neurons (%)") +
       coord_cartesian(ylim = c(0, 100)) +
       scale_color_manual(values = disease_color_dark, name = "diagnosis") +
       theme_bw(base_size = 8, base_family = "Helvetica") +
       theme(
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              strip.placement = "outside",
              legend.position = "none"
       )

ggsave(str_glue("./plots/compare_ExcNeurons_proportion_between_disease_diagnosis_fullDatasetClustering_NumGenes_{num_genes_cutoff}_DoubletScore_{doublet_score_cutoff}_PercentMt_{percent_mt_cutoff}.pdf"),
       device = cairo_pdf(), width = 2, height = 2, useDingbats = FALSE
)

## log sessionInfo
sessionInfo()
