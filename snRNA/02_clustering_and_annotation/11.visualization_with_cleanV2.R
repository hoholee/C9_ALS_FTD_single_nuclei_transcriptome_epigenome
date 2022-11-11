# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)
library(ggrastr)
library(scattermore)
library(cowplot)

# read in clean dataset with 2nd round annotation (V2, remove NK and Exc_TUBB2A)
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_cleanV2_SeuratV4_object.rds")

# extract umap coordinates
umap_coordinates <- Embeddings(data_obj_sub, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

# add to meta data
meta_data <- data_obj_sub@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble() %>% 
  left_join(umap_coordinates, by = "cell_id")

# shuffle the order of dots to alleviate overplotting
set.seed(seed = 666)
# meta_data_shuffle <- meta_data[sample(1:nrow(meta_data)), ]
meta_data_shuffle <- meta_data %>% sample_frac()

# manually plot UMAP plots colored by meta data
p <- meta_data_shuffle %>% ggplot(aes(UMAP_1, UMAP_2))

# brain region
region_labels <- c("MCX" = "motor cortex", "mFCX" = "medial frontal cortex")
region_color <- c("MCX" = "#432266", "mFCX" = "#FAA51B")
relabel_region <- function(x){region_labels[x]}
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
ggsave("./plots/umap_cleanV2_colored_by_region.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

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
ggsave("./plots/umap_cleanV2_colored_by_disease.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)


# sex
sex_labels <- c("F" = "female", "M" = "male")
sex_color <- c("F" = "#FC9E89", "M" = "#3D304C")
relabel_sex <- function(x){sex_labels[x]}
p_sex <- p +
  geom_scattermore(aes(color = sex), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(labels = relabel_sex(), values = sex_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_sex_panel <- p_sex + theme(legend.position = "none")
p_sex_legend <- get_legend(
  p_sex +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(p_sex_panel, p_sex_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_cleanV2_colored_by_sex.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)


# batch
batch_labels <- c("batch_1" = "1", "batch_2" = "2", "batch_3" = "3", "batch_4" = "4", "batch_5" = "5", "batch_6" = "6")
batch_color <- c("batch_1" = "#A8DBA8",
                 "batch_2" = "#00A0B0",
                 "batch_3" = "#6A4A3C",
                 "batch_4" = "#CC333F",
                 "batch_5" = "#EB6841",
                 "batch_6" = "#EDC951"
                 )
relabel_batch <- function(x){batch_labels[x]}
p_batch <- p +
  geom_scattermore(aes(color = seq_batch), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(labels = relabel_batch(), values = batch_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_batch_panel <- p_batch + theme(legend.position = "none")
p_batch_legend <- get_legend(
  p_batch +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom")
)

plot_grid(p_batch_panel, p_batch_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_cleanV2_colored_by_batch.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)


# individuals
individual_color <- c("332" = "#92072a",
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
ggsave("./plots/umap_cleanV2_colored_by_individual.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

# cell class (level 1 annotation)
cell_class_labels <- c("Non_neuron" = "non-neurons", "Exc_neuron" = "excitatory neurons", "Inh_neuron" = "inhibitory neurons")
cell_class_color <- c("Non_neuron" = "#9C482B", "Exc_neuron" = "#96BB45", "Inh_neuron" = "#718DC7")
relabel_cell_class <- function(x){cell_class_labels[x]}
p_cell_class <- p +
  geom_scattermore(aes(color = rna_anno_2ndRound_level_1), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Cell class", labels = relabel_cell_class(), values = cell_class_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_cell_class_panel <- p_cell_class + theme(legend.position = "none")
p_cell_class_legend <- get_legend(
  p_cell_class +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(p_cell_class_panel, p_cell_class_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_cleanV2_colored_by_cell_class.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

# extract sub-clustering umap coordinates
data_exc <- readRDS("./seurat_objects/snRNA_cellBender_corrected_dataset_SeuratV4_object_subset_Exc_scTransform_renormalized_annotated.rds")
data_inh <- readRDS("./seurat_objects/snRNA_cellBender_corrected_dataset_SeuratV4_object_subset_Inh_scTransform_renormalized_annotated.rds")
data_glia <- readRDS("./seurat_objects/snRNA_cellBender_corrected_dataset_SeuratV4_object_subset_NonNeu_scTransform_renormalized_annotated.rds")

umap_coordinates_exc <- Embeddings(data_exc, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble()
umap_coordinates_inh <- Embeddings(data_inh, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble()
umap_coordinates_glia <- Embeddings(data_glia, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

# add to meta data 
meta_data_exc <- data_obj_sub@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble() %>% 
  filter(rna_anno_2ndRound_level_1 == "Exc_neuron") %>% 
  left_join(umap_coordinates_exc, by = "cell_id")
meta_data_inh <- data_obj_sub@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble() %>% 
  filter(rna_anno_2ndRound_level_1 == "Inh_neuron") %>% 
  left_join(umap_coordinates_inh, by = "cell_id")
meta_data_glia <- data_obj_sub@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble() %>% 
  filter(rna_anno_2ndRound_level_1 == "Non_neuron") %>% 
  left_join(umap_coordinates_glia, by = "cell_id")


# shuffle the order of dots
set.seed(seed = 666)
meta_data_exc_shuffle <- meta_data_exc %>% sample_frac()
set.seed(seed = 666)
meta_data_inh_shuffle <- meta_data_inh %>% sample_frac()
set.seed(seed = 666)
meta_data_glia_shuffle <- meta_data_glia %>% sample_frac()

# manually plot UMAP plots, colored by subclusters (level 3 annotations)
p <- meta_data_exc_shuffle %>% ggplot(aes(UMAP_1, UMAP_2))

exc_palette <- read_tsv("color_palette_exc_subclusters.txt")
exc_subcluster_color <- exc_palette$color
names(exc_subcluster_color) <- exc_palette$sub_cluster
relabel_exc_subcluster <- function(x){str_replace_all(x, "_", " ")}

p_exc_subcluster <- p +
  geom_scattermore(aes(color = rna_anno_2ndRound_level_3), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Exc subcluster", labels = relabel_exc_subcluster, values = exc_subcluster_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_exc_subcluster_panel <- p_exc_subcluster + theme(legend.position = "none")
p_exc_subcluster_legend <- get_legend(
  p_exc_subcluster +
    guides(color = guide_legend(ncol = 4)) +
    theme(legend.position = "bottom")
)

plot_grid(p_exc_subcluster_panel, p_exc_subcluster_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_cleanV2_exc_subclustering_colored_by_subtype.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)


p <- meta_data_inh_shuffle %>% ggplot(aes(UMAP_1, UMAP_2))

inh_palette <- read_tsv("color_palette_inh_subclusters.txt")
inh_subcluster_color <- inh_palette$color
names(inh_subcluster_color) <- inh_palette$sub_cluster
relabel_inh_subcluster <- function(x){str_replace_all(x, "_", " ")}

p_inh_subcluster <- p +
  geom_scattermore(aes(color = rna_anno_2ndRound_level_3), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Inh subcluster", labels = relabel_inh_subcluster, values = inh_subcluster_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_inh_subcluster_panel <- p_inh_subcluster + theme(legend.position = "none")
p_inh_subcluster_legend <- get_legend(
  p_inh_subcluster +
    guides(color = guide_legend(ncol = 4)) +
    theme(legend.position = "bottom")
)

plot_grid(p_inh_subcluster_panel, p_inh_subcluster_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_cleanV2_inh_subclustering_colored_by_subtype.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

p <- meta_data_glia_shuffle %>% ggplot(aes(UMAP_1, UMAP_2))

glia_palette <- read_tsv("color_palette_glia_subclusters.txt")
glia_subcluster_color <- glia_palette$color
names(glia_subcluster_color) <- glia_palette$sub_cluster
relabel_glia_subcluster <- function(x){str_replace_all(x, "_", " ")}

p_glia_subcluster <- p +
  geom_scattermore(aes(color = rna_anno_2ndRound_level_3), pointsize = 4, pixels = c(2048, 2048)) +
  scale_color_manual(name = "Non-neuron subcluster", labels = relabel_glia_subcluster, values = glia_subcluster_color) +
  theme_classic(base_size = 8, base_family = "Helvetica")
p_glia_subcluster_panel <- p_glia_subcluster + theme(legend.position = "none")
p_glia_subcluster_legend <- get_legend(
  p_glia_subcluster +
    guides(color = guide_legend(ncol = 4)) +
    theme(legend.position = "bottom")
)

plot_grid(p_glia_subcluster_panel, p_glia_subcluster_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave("./plots/umap_cleanV2_glia_subclustering_colored_by_subtype.pdf", device = cairo_pdf(), dpi = 300, width = 2, height = 2.2)

# umap_nCount <- FeaturePlot(data_obj_sub, "nCount_RNA")
# umap_nFeature <- FeaturePlot(data_obj_sub, "nFeature_RNA")
# umap_percent_mt <- FeaturePlot(data_obj_sub, "percent_mt")
# umap_age <- FeaturePlot(data_obj_sub, "age")
# umap_PMI <- FeaturePlot(data_obj_sub, "PMI")
# umap_doubletScore <- FeaturePlot(data_obj_sub, "doublet_scores")
# 
# CombinePlots(plots = list(umap_nCount, umap_nFeature, umap_percent_mt, umap_doubletScore, umap_age, umap_PMI), ncol = 3)
# ggsave(paste0("umap_colored_by_QC_PCA_", pca_chosen, ".pdf"),
#        device = cairo_pdf(), width = 12, height = 8, dpi = 300) 

# build a tree using the `average` cell by each cluster, using PCs as input
data_obj_sub <- BuildClusterTree(data_obj_sub, dims = c(1:50))
Tool(object = data_obj_sub, slot = 'BuildClusterTree')
PlotClusterTree(data_obj_sub, direction = "rightwards")
Tool(object = data_obj_sub, slot = 'BuildClusterTree') -> tmp
library(dendextend)
tmp2 <- click_rotate(as.dendrogram(tmp))
tmp3 <- click_rotate(as.dendrogram(tmp2))
tmp4 <- click_rotate(as.dendrogram(tmp3))
tmp5 <- click_rotate(as.dendrogram(tmp4))
tmp6 <- click_rotate(as.dendrogram(tmp5))
tmp7 <- click_rotate(as.dendrogram(tmp6))
tmp8 <- click_rotate(as.dendrogram(tmp7))
tmp9 <- click_rotate(as.dendrogram(tmp8))
saveRDS(tmp9, "dendrogram_cleanV2_rotated.rds")
cluster_order <- labels(tmp9)
ggplot(as.ggdend(tmp9), horiz = TRUE) +theme_bw(base_size = 8, base_family = "Helvetica")
ggsave("./plots/cleanV2_dendrogram.pdf", device = cairo_pdf(), width = 7, height = 7)

# plot meta data side by side next to dendrogram
individual_order <- meta_data %>% distinct(disease, subject) %>% arrange(disease)
meta_data_sort <- meta_data %>% 
  mutate(cluster = factor(rna_anno_2ndRound_level_3, levels = cluster_order),
         subject = factor(subject, levels = individual_order$subject))
cluster_color <- c(exc_subcluster_color, inh_subcluster_color, glia_subcluster_color)

p <- meta_data_sort %>% ggplot(aes(cluster))

p_cellcount <- p + geom_bar(stat = "count", aes(fill = cluster)) +
  scale_fill_manual(values = cluster_color) +
  scale_y_log10()+
  theme_bw(base_size = 10, base_family = "Helvetica") +
  coord_flip() +
  annotation_logticks(sides = "b") +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

p_region <- p + geom_bar(aes(fill = region), stat = "count", position = position_fill()) +
  scale_fill_manual(name = "brain region", labels = relabel_region,  values = region_color) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  coord_flip() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

p_disease <- p + geom_bar(aes(fill = disease), stat = "count", position = position_fill()) +
  scale_fill_manual(values = disease_color) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  coord_flip() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

p_sex <- p + geom_bar(aes(fill = sex), stat = "count", position = position_fill()) +
  scale_fill_manual(labels = relabel_sex(), values = sex_color) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  coord_flip() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

p_batch <- p + geom_bar(aes(fill = seq_batch), stat = "count", position = position_fill()) +
  scale_fill_manual(labels = relabel_batch(), values = batch_color) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  coord_flip() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

p_individual <- p + geom_bar(aes(fill = subject), stat = "count", position = position_fill()) +
  scale_fill_manual(values = individual_color) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  coord_flip() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

p_num_count <-  p + geom_violin(aes(fill = cluster, y = nCount_RNA), scale = "width") +
  scale_fill_manual(values = cluster_color) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  coord_flip() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

p_num_gene <-  p + geom_violin(aes(fill = cluster, y = nFeature_RNA), scale = "width") +
  scale_fill_manual(values = cluster_color) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  coord_flip() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

plot_grid(p_cellcount, p_region, p_disease, p_sex, p_individual, p_batch, p_num_count, p_num_gene,
          align = "h", nrow = 1)
ggsave("./plots/cleanV2_meta.pdf", device = cairo_pdf(), width = 7, height = 7)

## plot marker genes
# level 1
data_obj_sub$rna_anno_2ndRound_level_1 <- factor(data_obj_sub$rna_anno_2ndRound_level_1,
                                                 levels = rev(c("Exc_neuron", "Inh_neuron", "Non_neuron")))

DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
        features = c("SNAP25", "MEF2C", "SLC17A7", "SATB2", "GAD1", "GAD2"
                     # "SOX10", "RBFOX3", "DCX", "MAP2", "TUBB3", "ENO2"
                     ),
        group.by = "rna_anno_2ndRound_level_1") +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")
ggsave("./plots/marker_genes_cell_class.pdf", device = cairo_pdf(), width = 3, height = 2, useDingbats = F)

# non-neurons 
data_obj_glia <- subset(data_obj_sub, subset = rna_anno_2ndRound_level_1 == "Non_neuron")

data_obj_glia$rna_anno_2ndRound_level_3 <- factor(data_obj_glia$rna_anno_2ndRound_level_3,
                                                 levels = rev(c("Astro_HPSE2", "Astro_CD44",
                                                                "Endo", "Micro",
                                                                "Oligo_ENPP6", "Oligo_OPALIN",
                                                                "OPC", "VLMC_CEMIP", "VLMC_P2RY14")))
DotPlot(data_obj_glia,
        cols = c("lightgrey", "#C73D4C"),
        features = c(
                     "AQP4", "APOE", "GFAP", "HPSE2", "CABLES1", "CD44", "WDR49", # Astro
                     "PECAM1", "VWF", "CLDN5", "NOSTRIN",  # Endo
                     "C1QB", "CSF1R", # Micro
                     "MOG", "MBP", "MOBP", "CNP", "ENPP6", "OPALIN",  # Oligo
                     "PDGFRA", "VCAN", # OPC
                     "TBX18", "COLEC12", "CEMIP", "CYP1B1", "P2RY14" # VLMC
                     ),
        group.by = "rna_anno_2ndRound_level_3") +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")

ggsave("./plots/marker_genes_glia_subclusters.pdf", device = cairo_pdf(), width = 6, height = 3.2, useDingbats = F)


# Inh neurons 
data_obj_inh <- subset(data_obj_sub, subset = rna_anno_2ndRound_level_1 == "Inh_neuron")

data_obj_inh$rna_anno_2ndRound_level_3 <- factor(data_obj_inh$rna_anno_2ndRound_level_3,
                                                 levels = rev(c("Inh_PVALB_CUX2", "Inh_PVALB_MYBPC1", "Inh_PVALB_PIEZO2", "Inh_PVALB_SCUBE3",
                                                                "Inh_SST_NPY", "Inh_SST_GPC5", "Inh_SST_EDNRA", "Inh_SST_CDH12", "Inh_SST_KLHL14",
                                                                "Inh_LAMP5_CHST9", "Inh_LAMP5_CPLX3_Rosehip", "Inh_LAMP5_NDNF",
                                                                "Inh_ADARB2_PAX6", "Inh_ADARB2_SEMA3C", "Inh_ADARB2_LINC01470", "Inh_ADARB2_SCML4",
                                                                "Inh_VIP_ABI3BP", "Inh_VIP_CLSTN2","Inh_VIP_DACH2","Inh_VIP_FLT1","Inh_VIP_ZBTB20")))
DotPlot(data_obj_inh,
        cols = c("lightgrey", "#C73D4C"),
        features = c(
          "GAD1", "GAD2",
          "LHX6", "SOX6",
          "PVALB", "CUX2", "MYBPC1", "PIEZO2", "SCUBE3",
          "SST", "NPY", "GPC5", "EDNRA", "CDH12", "KLHL14",
          "ADARB2", 
          "LAMP5", "SV2C", "CHST9", "CPLX3", "KIT", "CXCL14", "NDNF", "RELN", "PAX6", "SEMA3C", "LINC01470", "SCML4", 
          "VIP", "ABI3BP", "CLSTN2", "DACH2", "FLT1",  "ZBTB20"
        ),
        group.by = "rna_anno_2ndRound_level_3") +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")

ggsave("./plots/marker_genes_inh_subclusters.pdf", device = cairo_pdf(), width = 7.5, height = 5, useDingbats = F)


# Exc neurons 
data_obj_exc <- subset(data_obj_sub, subset = rna_anno_2ndRound_level_1 == "Exc_neuron")

data_obj_exc$rna_anno_2ndRound_level_3 <- factor(data_obj_exc$rna_anno_2ndRound_level_3,
                                                 levels = rev(c(
                                                                "Exc_L2_IT_CUX2_PDGFD",
                                                                "Exc_L2_IT_CUX2_CCBE1",
                                                                "Exc_L2_IT_CUX2_LRRC2", 
                                                                "Exc_L2_IT_CUX2_SV2C",
                                                                "Exc_L2-3_IT_RORB_PRSS12",
                                                                "Exc_L3_IT_RORB_OTOGL",
                                                                "Exc_L3-5_IT_RORB_GABRG1",
                                                                "Exc_L5_IT_RORB_NPFFR2",
                                                                "Exc_L3-5_IT_RORB_GRIN3A",
                                                                "Exc_L3-5_IT_RORB_ADAMTSL1",
                                                                "Exc_L3-5_IT_RORB_RPRM",
                                                                "Exc_L5-6_IT_THEMIS_SMYD1",
                                                                "Exc_L6_IT_THEMIS_LINC00299", 
                                                                "Exc_L6_IT_THEMIS_CFLAR",
                                                                "Exc_L6b_TLE4_MDFIC",
                                                                "Exc_L6b_TLE4_KCNK2",
                                                                "Exc_L6_CT_TLE4_SEMA5A",
                                                                "Exc_L5-6_NP_FOXP2_HTR2C",
                                                                "Exc_L5_ET_FEZF2_ADRA1A"
                                                 )))

DotPlot(data_obj_exc,
        cols = c("lightgrey", "#C73D4C"),
        features = c("SLC17A7", "SATB2",
                     "CUX2",
                     "ACVR1C",
                     "LAMP5", "SERPINE2", "LINC00507",
                     "PDGFD", "CCBE1",
                     "LRRC2","GLIS3",
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
        ),
        group.by = "rna_anno_2ndRound_level_3") +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")

ggsave("./plots/marker_genes_exc_subclusters.pdf", device = cairo_pdf(), width = 9.5, height = 4.9, useDingbats = F)

## log sessionInfo
sessionInfo()
