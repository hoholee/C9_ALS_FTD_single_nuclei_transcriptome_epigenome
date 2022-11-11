# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)
library(harmony)
library(ggrastr)
# library(future)

# AugmentPlot_mod <- AugmentPlot
# fix(AugmentPlot_mod) # delete noLegend()

# enable parallelization
# plan("multiprocess", workers = 8)
# options(future.globals.maxSize = 2000 * 1024^2)

## read the sctransform-normalized seurat object
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBender_postQC_scTransformed_SeuratV4_object.rds")

## run PCA
data_obj_sub <- RunPCA(data_obj_sub, features = VariableFeatures(object = data_obj_sub), npcs = 50, seed.use = 666)
# quick check on PCs
print(data_obj_sub[["pca"]], dims = 1:5, nfeature = 5)
VizDimLoadings(data_obj_sub, dims = 1:4, reduction = "pca", nfeatures = 10)
ggsave("./plots/PCA_all_cellBender_corrected_QCed_scTransformed_cells_top4_PCs_loadings.pdf", device = cairo_pdf(),
       width = 8, height = 4, useDingbats = F)

DimPlot(data_obj_sub, reduction = "pca",  ncol = 6, raster = TRUE)
ggsave("./plots/PCA_all_cellBender_corrected_QCed_scTransformed_cells_PCs_dotplot.pdf", device = cairo_pdf(),
       width = 8, height = 4, useDingbats = F)
DimPlot(data_obj_sub, reduction = "pca", split.by = "orig.ident", ncol = 6)
ggsave("./plots/PCA_all_cellBender_corrected_QCed_scTransformed_cells_PCs_dotplot_splitBySample.pdf", device = cairo_pdf(),
       width = 13, height = 8, useDingbats = F)

pdf(file = "./plots/PCA_all_cellBender_corrected_QCed_scTransformed_cells_PCs_heatmap.pdf", width = 13, height = 8)
DimHeatmap(data_obj_sub, dims = 1:6, cells = 500, balanced = TRUE, nfeatures = 100)
dev.off()

## determine the 'dimensionality' of the dataset (select how many PCs to include in the downstream analysis)
# jackStarw method is not appropriate for scTransfromed data; use elbow plot instead
# data_obj_sub <- JackStraw(data_obj_sub, num.replicate = 100, dims = 50)
# data_obj_sub <- ScoreJackStraw(data_obj_sub, dims = 1:50)
# JackStrawPlot(data_obj_sub, dims = 1:50)

ElbowPlot(data_obj_sub, ndims = 50)
ggsave("./plots/PCA_all_cellBender_corrected_QCed_scTransformed_cells_elbowPlot.pdf", device = cairo_pdf(),
       width = 6, height = 4, useDingbats = F)
pca_chosen <- 30

# run harmony
data_obj_sub <- RunHarmony(data_obj_sub,
                           group.by.vars = c("orig.ident", "subject", "seq_batch"),
                           dims.use = 1:pca_chosen,
                           assay.use = "SCT",
                           kmeans_init_nstart = 20,
                           kmeans_init_iter_max = 5000)

## clustering
data_obj_sub <- FindNeighbors(data_obj_sub,
                              reduction = "harmony",
                              dims = 1:pca_chosen,
                              k.param = 30)
# use the Leiden algorithm
# switch back to single worker to avoid serialize error
# plan("multiprocess", workers = 1)
# took ~1.5h to finish one resolution
data_obj_sub <- FindClusters(data_obj_sub,
                             # resolution = c(0.5, 0.8, 1, 1.2, 1.4, 2, 3),
                             resolution = 1,
                             # algorithm = 1, # 1 = original Louvain
                             algorithm = 4, # 4 = Leiden
                             method = "igraph", # modify from "matrix" to "igraph" for large dataset
                             random.seed = 666,
                             verbose = TRUE)
# quick check on the clusters
Idents(data_obj_sub) %>% head()

## run UMAP/tSNE
data_obj_sub <- RunUMAP(data_obj_sub,
                        umap.method = "uwot",
                        metric = "cosine",
                        n.components = 2,
                        n.neighbors = 30,
                        dims = 1:pca_chosen,
                        reduction = "harmony",
                        # min.dist = 0.7,
                        a = 2,
                        b = 1.5,
                        seed.use = 666)


## building a tree relating the ‘average’ cell from each cluster
#
# data_obj_sub <- BuildClusterTree(data_obj_sub,
#                                  dims = c(1:pca_chosen),
#                                  reorder = FALSE)
# pdf(file = paste0("clustering_tree_PCA_", pca_chosen, ".pdf"), width = 8, height = 8)
# PlotClusterTree(data_obj_sub)
# dev.off()

## check batch effect on UMAP
meta_data <- data_obj_sub@meta.data %>% rownames_to_column("cell_id") %>% as_tibble()

umap_coordinates <- data_obj_sub@reductions$umap@cell.embeddings %>%
        as.data.frame() %>%
        rownames_to_column("cell_id") %>%
        as_tibble()
meta_data_mod <- meta_data %>% left_join(umap_coordinates, by = c("cell_id"))

p <- ggplot(meta_data_mod, aes(UMAP_1, UMAP_2))
p + geom_point_rast(size = 0.1, aes(color = seurat_clusters), raster.dpi = 300, scale = 0.3) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  # scale_color_manual(values = color_used) +
  theme(panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 6)))
ggsave("./plots/umap_cellBender_corrected_colored_by_cluster_r1.pdf", device = cairo_pdf(), width = 6, height = 4, dpi = 300)

p <- ggplot(meta_data_mod, aes(UMAP_1, UMAP_2))
p + geom_point_rast(size = 0.1, aes(color = region), raster.dpi = 300, scale = 0.3) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  scale_color_manual(values = c("MCX" = "#6A2D70", "mFCX" = "#F08A5C")) +
  theme(panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 6)))
ggsave("./plots/umap_cellBender_corrected_colored_by_region.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

p <- ggplot(meta_data_mod, aes(UMAP_1, UMAP_2))
p + geom_point_rast(size = 0.1, aes(color = disease), raster.dpi = 300, scale = 0.3) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  # scale_color_manual(values = c("ALS" = "#5E70B6", "FTD" = "#DE9D4D", "Control" = "#9A9393")) +
  scale_color_manual(values = c("ALS" = "#5E70B6", "Control" = "#9A9393", "FTD" = "#DE9D4D")) +
  theme(panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 6)))
ggsave("./plots/umap_cellBender_corrected_colored_by_disease.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)


umap_nCount <- FeaturePlot(data_obj_sub, "nCount_RNA")
umap_nFeature <- FeaturePlot(data_obj_sub, "nFeature_RNA")
umap_percent_mt <- FeaturePlot(data_obj_sub, "percent_mt")
umap_age <- FeaturePlot(data_obj_sub, "age")
umap_PMI <- FeaturePlot(data_obj_sub, "PMI")
umap_doubletScore <- FeaturePlot(data_obj_sub, "doublet_scores")

CombinePlots(plots = list(umap_nCount, umap_nFeature, umap_percent_mt, umap_doubletScore, umap_age, umap_PMI), ncol = 3)
ggsave(paste0("./plots/umap_cellBender_corrected_colored_by_QC_PCA_", pca_chosen, ".pdf"),
       device = cairo_pdf(), width = 12, height = 8, dpi = 300)


# save Seurat object
saveRDS(data_obj_sub, file = "./seurat_objects/snRNA_cellBender_corrected_postQC_scTransformed_clustered_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
