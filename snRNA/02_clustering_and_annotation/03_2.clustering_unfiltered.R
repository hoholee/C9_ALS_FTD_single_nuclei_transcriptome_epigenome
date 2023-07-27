# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)
library(harmony)
library(ggrastr)

## read the sctransform-normalized seurat object
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBender_corrected_noQCFilters_scTransformed_SeuratV4_object.rds")

## run PCA
data_obj_sub <- RunPCA(data_obj_sub, features = VariableFeatures(object = data_obj_sub), npcs = 50, seed.use = 666)
# quick check on PCs
print(data_obj_sub[["pca"]], dims = 1:5, nfeature = 5)
VizDimLoadings(data_obj_sub, dims = 1:4, reduction = "pca", nfeatures = 10)
ggsave("./plots/PCA_all_cellBender_corrected_noQCFilters_scTransformed_cells_top4_PCs_loadings.pdf",
       device = cairo_pdf(),
       width = 8, height = 4, useDingbats = F
)

DimPlot(data_obj_sub, reduction = "pca", ncol = 6, raster = TRUE)
ggsave("./plots/PCA_all_cellBender_corrected_noQCFilters_scTransformed_cells_PCs_dotplot.pdf",
       device = cairo_pdf(),
       width = 8, height = 4, useDingbats = F
)
DimPlot(data_obj_sub, reduction = "pca", split.by = "orig.ident", ncol = 6)
ggsave("./plots/PCA_all_cellBender_corrected_noQCFilters_scTransformed_cells_PCs_dotplot_splitBySample.pdf",
       device = cairo_pdf(),
       width = 13, height = 8, useDingbats = F
)

pdf(file = "./plots/PCA_all_cellBender_corrected_noQCFilters_scTransformed_cells_PCs_heatmap.pdf", width = 13, height = 8)
DimHeatmap(data_obj_sub, dims = 1:6, cells = 500, balanced = TRUE, nfeatures = 100)
dev.off()

## determine the 'dimensionality' of the dataset (select how many PCs to include in the downstream analysis)
# jackStarw method is not appropriate for scTransfromed data; use elbow plot instead
# data_obj_sub <- JackStraw(data_obj_sub, num.replicate = 100, dims = 50)
# data_obj_sub <- ScoreJackStraw(data_obj_sub, dims = 1:50)
# JackStrawPlot(data_obj_sub, dims = 1:50)

ElbowPlot(data_obj_sub, ndims = 50)
ggsave("./plots/PCA_all_cellBender_corrected_noQCFilters_scTransformed_cells_elbowPlot.pdf",
       device = cairo_pdf(),
       width = 6, height = 4, useDingbats = F
)
pca_chosen <- 30

# run harmony
data_obj_sub <- RunHarmony(data_obj_sub,
       group.by.vars = c("orig.ident", "subject", "seq_batch"),
       dims.use = 1:pca_chosen,
       assay.use = "SCT",
       kmeans_init_nstart = 20,
       kmeans_init_iter_max = 5000
)

## clustering
data_obj_sub <- FindNeighbors(data_obj_sub,
       reduction = "harmony",
       dims = 1:pca_chosen,
       k.param = 30
)
# use the Leiden algorithm
# took ~1.5h to finish one resolution
data_obj_sub <- FindClusters(data_obj_sub,
       # resolution = c(0.5, 0.8, 1, 1.2, 1.4, 2, 3),
       resolution = 1,
       # algorithm = 1, # 1 = original Louvain
       algorithm = 4, # 4 = Leiden
       method = "igraph", # modify from "matrix" to "igraph" for large dataset
       random.seed = 666,
       verbose = TRUE
)
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
       seed.use = 666
)

# save Seurat object
saveRDS(data_obj_sub, file = "./seurat_objects/snRNA_cellBender_corrected_noQCFilters_scTransformed_clustered_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
