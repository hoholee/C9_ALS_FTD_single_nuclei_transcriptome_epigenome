# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load librar
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)
library(harmony)

## read the sctransform-normalized seurat object
selected_celltype <- "Inh"

data_obj_sub <- readRDS(paste0("./seurat_objects/snRNA_cellBender_corrected_dataset_SeuratV4_object_subset_", selected_celltype, "_scTransform_renormalized.rds"))

## run PCA
data_obj_sub <- RunPCA(data_obj_sub, features = VariableFeatures(object = data_obj_sub), npcs = 50, seed.use = 666)
# quick check on PCs
print(data_obj_sub[["pca"]], dims = 1:5, nfeature = 5)

VizDimLoadings(data_obj_sub, dims = 1:4, reduction = "pca", nfeatures = 10)
ggsave(paste0("./plots/PCA_subset_", selected_celltype, "_cellBender_corrected_QCed_scTransformed_cells_top4_PCs_loadings.pdf"), device = cairo_pdf(),
       width = 8, height = 4, useDingbats = F)

DimPlot(data_obj_sub, reduction = "pca", shuffle = TRUE, raster = TRUE)
ggsave(paste0("./plots/PCA_subset_", selected_celltype, "_cellBender_corrected_QCed_scTransformed_cells_PCs_dotplot.pdf"), device = cairo_pdf(),
       width = 8, height = 4, useDingbats = F)

DimPlot(data_obj_sub, reduction = "pca", split.by = "orig.ident", ncol = 6)
ggsave(paste0("./plots/PCA_subset_", selected_celltype, "_cellBender_corrected_QCed_scTransformed_cells_PCs_dotplot_splitBySample.pdf"), device = cairo_pdf(),
       width = 13, height = 8, useDingbats = F)

pdf(file = paste0("./plots/PCA_subset_", selected_celltype, "_cellBender_corrected_QCed_scTransformed_cells_PCs_heatmap.pdf"),
    width = 13, height = 8)
DimHeatmap(data_obj_sub, dims = 1:6, cells = 500, balanced = TRUE, nfeatures = 100)
dev.off()

## determine the 'dimensionality' of the dataset (select how many PCs to include in the downstream analysis)
# jackStarw method is not appropriate for scTransfromed data; use elbow plot instead
# data_obj_sub <- JackStraw(data_obj_sub, num.replicate = 100, dims = 50)
# data_obj_sub <- ScoreJackStraw(data_obj_sub, dims = 1:50)
# JackStrawPlot(data_obj_sub, dims = 1:50)

ElbowPlot(data_obj_sub, ndims = 50)
ggsave(paste0("./plots/PCA_subbset_", selected_celltype, "_cellBender_corrected_QCed_scTransformed_cells_elbowPlot.pdf"),
       device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

pca_chosen <- 30
k_chosen <- 30
r_chosen <- 0.5
# r_chosen <- 1
# r_chosen <- 3
a_chosen <- 1.5 
b_chosen <- 2

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
                              k.param = k_chosen)
# use the Leiden algorithm
# switch back to single worker to avoid serialize error
# plan("multiprocess", workers = 1)
# took ~1.5h to finish one resolution
data_obj_sub <- FindClusters(data_obj_sub,
                             resolution = r_chosen,
                             algorithm = 4, # 4 = Leiden
                             method = "matrix", # modify from "matrix" to "igraph" for large dataset
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
                        a = a_chosen,
                        b = b_chosen,
                        seed.use = 666)


DimPlot(data_obj_sub, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, raster = F, label = T)
DimPlot(data_obj_sub, reduction = "umap", group.by = "rna_anno_1stRound_level_2", shuffle = TRUE, raster = F, label = T)

g1 <- WhichCells(data_obj_sub, expression = seurat_clusters == "12")
DimPlot(data_obj_sub, label=T,  cells.highlight= list(g1), cols.highlight = c("darkred"), cols= "grey", raster = FALSE)

to_plot_resolution <- r_chosen

p <- DimPlot(data_obj_sub, reduction = "umap",
             group.by = paste0("SCT_snn_res.", to_plot_resolution),
             label = TRUE, pt.size = 0.2, shuffle = T, raster = F) +
  ggtitle(paste0(selected_celltype, ", PCs: ", pca_chosen,
                 ", r: ", to_plot_resolution,
                 ", k: ", k_chosen, ", a: ", a_chosen, ", b: ", b_chosen))
ggsave(paste0("./plots/umap_cellBender_corrected_", selected_celltype, "_subclustering_colored_by_clusters_r_", to_plot_resolution,".pdf"),
       plot = p, device = cairo_pdf(), width = 6, height = 4, dpi = 300) 

## building a tree relating the ‘average’ cell from each cluster
# data_obj_sub <- BuildClusterTree(data_obj_sub,
#                                  dims = c(1:pca_chosen),
#                                  reorder = FALSE)
# pdf(file = paste0("clustering_tree_", selected_celltype, ".pdf"), width = 8, height = 5)
# PlotClusterTree(data_obj_sub)
# dev.off()

## check batch effect on UMAP
# p <- DimPlot(data_obj_sub, reduction = "umap", group.by = "subject", pt.size = 0.2)
# p <- AugmentPlot(p, width = 6, height = 4, dpi = 300)
# ggsave(paste0("umap_colored_by_subject_", selected_celltype, ".pdf"), plot = p, device = cairo_pdf(), width = 6, height = 4, dpi = 300) 
# 
# p <- DimPlot(data_obj_sub, reduction = "umap", group.by = "region", pt.size = 0.2)
# p <- AugmentPlot_mod(p, width = 6, height = 4, dpi = 300)
# ggsave(paste0("umap_colored_by_region_", selected_celltype, ".pdf"), plot = p, device = cairo_pdf(), width = 6, height = 4, dpi = 300) 
# 
# p <- DimPlot(data_obj_sub, reduction = "umap", group.by = "disease", pt.size = 0.2)
# p <- AugmentPlot_mod(p, width = 6, height = 4, dpi = 300)
# ggsave(paste0("umap_colored_by_disease_", selected_celltype, ".pdf"), plot = p, device = cairo_pdf(), width = 6, height = 4, dpi = 300) 
# 
# p <- DimPlot(data_obj_sub, reduction = "umap", group.by = "sex", pt.size = 0.2)
# p <- AugmentPlot_mod(p, width = 6, height = 4, dpi = 300)
# ggsave(paste0("umap_colored_by_sex_", selected_celltype, ".pdf"), plot = p, device = cairo_pdf(), width = 6, height = 4, dpi = 300) 
# 
# p <- DimPlot(data_obj_sub, reduction = "umap", group.by = "seq_batch", pt.size = 0.2)
# p <- AugmentPlot_mod(p, width = 6, height = 4, dpi = 300)
# ggsave(paste0("umap_colored_by_seqBatch_", selected_celltype, ".pdf"), plot = p, device = cairo_pdf(), width = 6, height = 4, dpi = 300) 
# 
# umap_nCount <- FeaturePlot(data_obj_sub, "nCount_RNA")
# umap_nFeature <- FeaturePlot(data_obj_sub, "nFeature_RNA")
# umap_percent_mt <- FeaturePlot(data_obj_sub, "percent_mt")
# umap_age <- FeaturePlot(data_obj_sub, "age")
# umap_PMI <- FeaturePlot(data_obj_sub, "PMI")
# umap_doubletScore <- FeaturePlot(data_obj_sub, "doublet_scores")
# 
# CombinePlots(plots = list(umap_nCount, umap_nFeature, umap_percent_mt, umap_doubletScore, umap_age, umap_PMI), ncol = 3)
# ggsave(paste0("umap_colored_by_QC_", selected_celltype, ".pdf"),
#        device = cairo_pdf(), width = 12, height = 8, dpi = 300) 

# annotation
p <- FeaturePlot(data_obj_sub, c("AQP4", "MOG", "ADARB2", 
                            "LHX6", "SOX6", "SST",
                            "PVALB", "VIP", "LAMP5"), ncol = 3)
# ggsave(paste0("umap_colored_by_markers_", selected_celltype, ".pdf"),
#        plot = p, device = cairo_pdf(), width = 12, height = 8, dpi = 300) 
  
# test marker genes between some ambigious clusters
markers_cluster_2 <- FindMarkers(data_obj_sub,
                                   ident.1 = 2,
                                   min.pct = 0.25)

markers_cluster_17 <- FindMarkers(data_obj_sub,
                                   ident.1 = 17,
                                   min.pct = 0.25)

markers_cluster_21 <- FindMarkers(data_obj_sub,
                                   ident.1 = 21,
                                   min.pct = 0.25)

markers_cluster_9 <- FindMarkers(data_obj_sub,
                                   ident.1 = 9,
                                   min.pct = 0.25)

markers_2_out <- markers_cluster_2 %>% rownames_to_column("gene") %>% as_tibble()
markers_17_out <- markers_cluster_17 %>% rownames_to_column("gene") %>% as_tibble()
markers_21_out <- markers_cluster_21 %>% rownames_to_column("gene") %>% as_tibble()
markers_9_out <- markers_cluster_9 %>% rownames_to_column("gene") %>% as_tibble()

write_tsv(markers_2_out, paste0("marker_genes_Inh_subclustering_r_", r_chosen, "_cluster_2.txt"))
write_tsv(markers_17_out, paste0("marker_genes_Inh_subclustering_r_", r_chosen, "_cluster_17.txt"))
write_tsv(markers_21_out, paste0("marker_genes_Inh_subclustering_r_", r_chosen, "_cluster_21.txt"))
write_tsv(markers_9_out, paste0("marker_genes_Inh_subclustering_r_", r_chosen, "_cluster_9.txt"))

DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
        group.by = paste0("SCT_snn_res.", to_plot_resolution),
        features = c(
                     "MOG", "MOBP", "CNP",
                     "AQP4", "APOE", "GFAP", 
                     "SNAP25", "MEF2C",
                     "ADARB2", "LHX6", "SOX6",
                     "VIP", "LAMP5", "RELN", "PAX6", "NDNF", "SV2C", "CPLX3", "CHST9", "FAM19A1",
                     "DACH2", "CLSTN2", "FLT1", "ABI3BP", "ZBTB20", "CXCL14", "KIT",
                     "PVALB", "SST", "SCUBE3", "PIEZO2", "CUX2", "MAML3", "SGCZ",
                     "CDH12", "NPY", "KLHL14", "SPHKAP", "MAFB", "MYBPC1", "HSPA8",  "NRGN", "ITM2C",
                     "CADPS2", "SV2B", "CA10",
                     "LINC02240", "LINC01470", "DCHS2", "TAC3",
                     "PDE1A", "PAWR", "GRM1", "TANC1", "EDNRA", "EPB41L4A", "GPC5", "SCML4", "PTGER3", "KCNG1", "SEMA3C"
                     )
        ) +
  RotatedAxis()
ggsave(paste0("./plots/dotplot_markers_", selected_celltype, "_r_", to_plot_resolution, ".pdf"),
       device = cairo_pdf(), width = 12, height = 6, dpi = 300)

# add level 3 annotation
meta_data <- data_obj_sub@meta.data %>% rownames_to_column("cell_id") %>% as_tibble()

meta_data_annotated <- meta_data %>% 
  mutate(
    rna_anno_level_3 = case_when(
      seurat_clusters == 22 ~ "DB_Inh_Oligo",
      seurat_clusters == 13 ~ "DB_Inh_Oligo_Astro",
      seurat_clusters == 23 ~ "Inh_SST_NPY",
      seurat_clusters == 15 ~ "Inh_SST_KLHL14",
      seurat_clusters == 17 ~ "Inh_SST_GPC5",
      seurat_clusters == 2 ~ "Inh_SST_EDNRA",
      seurat_clusters == 3 ~ "Inh_SST_CDH12",
      seurat_clusters == 1 ~ "Inh_PVALB_CUX2",
      seurat_clusters == 5 ~ "Inh_PVALB_PIEZO2",
      seurat_clusters == 19 ~ "Inh_PVALB_MYBPC1",
      seurat_clusters == 10 ~ "Inh_PVALB_SCUBE3",
      seurat_clusters == 11 ~ "Inh_LAMP5_CHST9",
      seurat_clusters == 4 ~ "Inh_LAMP5_CPLX3_Rosehip",
      seurat_clusters == 14 ~ "Inh_LAMP5_NDNF",
      seurat_clusters == 18 ~ "Inh_ADARB2_PAX6",
      seurat_clusters == 9 ~ "Inh_ADARB2_SEMA3C",
      seurat_clusters == 20 ~ "Inh_ADARB2_LINC01470",
      seurat_clusters == 21 ~ "Inh_ADARB2_SCML4",
      seurat_clusters == 8 ~ "Inh_VIP_ABI3BP",
      seurat_clusters == 7 ~ "Inh_VIP_CLSTN2",
      seurat_clusters == 6 ~ "Inh_VIP_FLT1",
      seurat_clusters == 12 ~ "Inh_VIP_ZBTB20",
      seurat_clusters == 16 ~ "Inh_VIP_DACH2"
      )
  )

data_obj_sub$rna_anno_level_3 <- meta_data_annotated$rna_anno_level_3

# save annotation and clusters
meta_out <- meta_data_annotated %>%
  select(cell_id, starts_with("rna_anno_"), seurat_clusters) %>% 
  rename(seurat_sub_clusters = seurat_clusters)
write_tsv(meta_out, paste0("meta_annotated_", selected_celltype, "_sub_clusters.txt"))

## rename idents
Idents(data_obj_sub) <- "rna_anno_level_3"

# save Seurat object
saveRDS(data_obj_sub, file = paste0("./seurat_objects/snRNA_cellBender_corrected_dataset_SeuratV4_object_subset_", selected_celltype, "_scTransform_renormalized_annotated.rds"))

## log sessionInfo
sessionInfo()
