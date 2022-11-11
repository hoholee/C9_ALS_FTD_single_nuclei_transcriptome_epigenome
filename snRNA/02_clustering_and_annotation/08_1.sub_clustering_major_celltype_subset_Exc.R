# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)
library(harmony)

## read the sctransform-normalized seurat object
selected_celltype <- "Exc"

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
# r_chosen <- 0.5
# r_chosen <- 1
r_chosen <- 1.5
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

g1 <- WhichCells(data_obj_sub, expression = seurat_clusters == "21")
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
FeaturePlot(data_obj_sub, c("CUX2", "RORB", "FOXP2",
                            "TSHZ2", "THEMIS", "TSHZ1",
                            "MOG", "AQP4", "APOE"), ncol = 3)
# ggsave(paste0("umap_colored_by_markers_", selected_celltype, ".pdf"),
#        plot = p, device = cairo_pdf(), width = 12, height = 8, dpi = 300) 
   
# test marker genes between some ambiguous clusters
markers_cluster_2 <- FindMarkers(data_obj_sub,
                                   ident.1 = 2,
                                   min.pct = 0.25)

markers_cluster_5 <- FindMarkers(data_obj_sub,
                                   ident.1 = 5,
                                   min.pct = 0.25)

markers_cluster_8 <- FindMarkers(data_obj_sub,
                                   ident.1 = 8,
                                   min.pct = 0.25)

markers_cluster_18 <- FindMarkers(data_obj_sub,
                                   ident.1 = 18,
                                   min.pct = 0.25)

markers_cluster_7 <- FindMarkers(data_obj_sub,
                                   ident.1 = 7,
                                   min.pct = 0.25)

markers_cluster_21 <- FindMarkers(data_obj_sub,
                                   ident.1 = 21,
                                   min.pct = 0.25)

markers_2_out <- markers_cluster_2 %>% rownames_to_column("gene") %>% as_tibble()
markers_5_out <- markers_cluster_5 %>% rownames_to_column("gene") %>% as_tibble()
markers_8_out <- markers_cluster_8 %>% rownames_to_column("gene") %>% as_tibble()
markers_18_out <- markers_cluster_18 %>% rownames_to_column("gene") %>% as_tibble()
markers_7_out <- markers_cluster_7 %>% rownames_to_column("gene") %>% as_tibble()
markers_21_out <- markers_cluster_21 %>% rownames_to_column("gene") %>% as_tibble()

write_tsv(markers_2_out, paste0("marker_genes_Exc_subclustering_r_", r_chosen, "_cluster_2.txt"))
write_tsv(markers_5_out, paste0("marker_genes_Exc_subclustering_r_", r_chosen, "_cluster_5.txt"))
write_tsv(markers_8_out, paste0("marker_genes_Exc_subclustering_r_", r_chosen, "_cluster_8.txt"))
write_tsv(markers_18_out, paste0("marker_genes_Exc_subclustering_r_", r_chosen, "_cluster_18.txt"))
write_tsv(markers_7_out, paste0("marker_genes_Exc_subclustering_r_", r_chosen, "_cluster_7.txt"))
write_tsv(markers_21_out, paste0("marker_genes_Exc_subclustering_r_", r_chosen, "_cluster_21.txt"))

DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
        group.by = paste0("SCT_snn_res.", to_plot_resolution),
        features = c("MOG", "MOBP", "MBP", "CNP", "PLP1", 
                     "AQP4", "APOE",
                     "SNAP25", "MEF2C", "SLC17A7",
                     "CUX2", "LAMP5", "BMPR1B", "SERPINE2", "LINC00507", "PDGFD",
                     "CCBE1",
                     "COL5A2", "SEMA3A", "LRRK1", "LRRC2", "PRKG2", "GLIS3",
                     "TSHZ1", "SASH1",
                     "SV2C", "SYT2",
                     "ACVR1C", "BDNF", "ATP7B",
                     "RORB", "COBLL1", "PLCH1", "PRSS12", "CCDC68",
                     "COL22A1", "OTOGL",
                     "AC109466.1", "AC008415.1",
                     "FOXP2", "AC092422.1", "ALDH1A1", "GABRG1", "AL109930.1", "LINC02232",
                     "THEMIS", "DACH1", "CRH", "NPFFR2",
                     "LNX2", "TRABD2A", "SLIT3",
                     "RPRM", "SEMA3E",
                     "PDZRN4", "ANKRD30B", "LINC00299",
                     "TRPC5", "ATP10A", "LTBP1", "SMYD1",
                     "TOX", "TLE4", "FEZF2", "HTR2C", "IFNG-AS1", "KCNIP1",
                     "MDFIC", "PCOLCE2", "PCSK5",
                     "SULF1", "SEMA5A",
                     "ADRA1A", "ERG", "MYO16", "POU3F1", "NEFH", "CRYM",
                     "ITGA8", "KCNK2",
                     "GNB4", "CFLAR",
                     "GRIK1", "SYT10", "LINC00326", "CNGB1", "ADGRF5", "GRB14", "LINC02296", "GRIN3A", "ADAMTSL1", "TLL1"
                     )
        ) +
  # scale_color_gradientn(colours = colorRampPalette(scico(100, palette = "acton", direction = -1))(100), limit = c(0, 3), oob = squish) +
  RotatedAxis()
# ggsave(paste0("./plots/dotplot_markers_", selected_celltype, "_r_", to_plot_resolution, ".pdf"),
#        device = cairo_pdf(), width = 12, height = 6, dpi = 300)

# add level 3 annotation
meta_data <- data_obj_sub@meta.data %>% rownames_to_column("cell_id") %>% as_tibble()

meta_data_annotated <- meta_data %>% 
  mutate(
    rna_anno_level_3 = case_when(
      seurat_clusters == 19 ~ "DB_Exc_Oligo",
      seurat_clusters == 24 ~ "DB_Exc_Astro",
      seurat_clusters == 21 ~ "DB_Exc_Inh",
      seurat_clusters == 7 ~ "Exc_TUBB2A", # TUBB2A, TUBB4B, TUBA4A, (low number of genes detected)
      seurat_clusters == 1 ~ "Exc_L2_IT_CUX2_PDGFD", # CUX2, LAMP5, BMPR1B, SERPINE2, LINC00507, PDGFD
      seurat_clusters %in% c(4, 5) ~ "Exc_L2_IT_CUX2_CCBE1", # CUX2, LAMP5, BMPR1B, SERPINE2, LINC00507, CCBE1
      seurat_clusters == 2 ~ "Exc_L2_IT_CUX2_LRRC2", # CUX2, COL5A2, SEMA3A, LRRK1, LRRC2, PRKG2, LINC00507, GLIS3, LINC00326
      seurat_clusters == 9 ~ "Exc_L2_IT_CUX2_SV2C", # CUX2, SV2C, SYT2, BMPR1B
      seurat_clusters == 12 ~ "Exc_L2-3_IT_RORB_PRSS12", # CUX2, RORB, PLCH1, COBLL1, COL5A2, PLCH1, PRSS12, CCDC68
      seurat_clusters == 11 ~ "Exc_L3_IT_RORB_OTOGL", # RORB, PLCH1, COBLL1, COL22A1, OTOGL
      seurat_clusters == 6 ~ "Exc_L3-5_IT_RORB_GABRG1", # RORB, COBLL1, FOXP2, AC092422.1, ALDH1A1, GABRG1, AL109930.1, LINC02232 
      seurat_clusters == 17 ~ "Exc_L5_IT_RORB_NPFFR2", # RORB, THEMIS, DACH1, COBLL1, CRH, BDNF, NPFFR2
      seurat_clusters == 8 ~ "Exc_L3-5_IT_RORB_GRIN3A", # RORB, DACH1, LRRK1, FOXP2, LNX2, TRABD2A, SLIT3, AC008415.1, GRIN3A
      seurat_clusters == 18 ~ "Exc_L3-5_IT_RORB_ADAMTSL1", # RORB, DACH1, LRRK1, FOXP2, LNX2, TRABD2A, AC008415.1, ADAMTSL1
      seurat_clusters == 13 ~ "Exc_L3-5_IT_RORB_RPRM", # RORB, DACH1, LRRK1, FOXP2, LNX2, RPRM, SEMA3E
      seurat_clusters == 3 ~ "Exc_L6_IT_THEMIS_LINC00299", # THEMIS, FOXP2, PDZRN4, ANKRD30B, LINC00299
      seurat_clusters == 23 ~ "Exc_L6_IT_THEMIS_CFLAR", # THEMIS, PDZRN4, ANKRD30B, TRPC5, GNB4, CFLAR, 
      seurat_clusters == 14 ~ "Exc_L5-6_IT_THEMIS_SMYD1", # THEMIS, PDZRN4, ANKRD30B, TRPC5, ATP10A, LTBP1, SMYD1
      seurat_clusters == 15 ~ "Exc_L5-6_NP_FOXP2_HTR2C", # FOXP2, PDZRN4, TLE4, TOX, FEZF2, HTR2C, IFNG-AS1, ITGA8, KCNIP1, LHFPL3, NXPH2, ROR1, TLL1, SLC24A3, CD36, ADAMTS12
      seurat_clusters == 16 ~ "Exc_L6b_TLE4_MDFIC", # FOXP2, PDZRN4, TLE4, MDFIC, PCOLCE2, PCSK5, ITGA8
      seurat_clusters == 20 ~ "Exc_L6b_TLE4_KCNK2", # FOXP2, PDZRN4, TLE4, PCOLCE2, PCSK5, SULF1, KCNK2
      seurat_clusters == 10 ~ "Exc_L6_CT_TLE4_SEMA5A", # FOXP2, PDZRN4, TLE4, TOX, SULF1, SEMA5A
      seurat_clusters == 22 ~ "Exc_L5_ET_FEZF2_ADRA1A" # FEZF2, ADRA1A, ERG, MYO16,  
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
