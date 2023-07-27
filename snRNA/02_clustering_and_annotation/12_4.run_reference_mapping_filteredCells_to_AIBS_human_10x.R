# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)

filtered_cells_obj <- readRDS("./seurat_objects/snRNA_cellBender_corrected_QC_filteredCells_scTransformed_SeuratV4_object.rds")
aibs_obj <- readRDS("/cndd2/junhao/AIBS_human_M1_10x/AIBS_human_M1_10x_Seurat_obj_scTransformed_PCA50_clustered.rds")

anchors <- FindTransferAnchors(
    reference = aibs_obj,
    query = filtered_cells_obj,
    normalization.method = "SCT",
    query.assay = "RNA",
    dims = 1:50
)
predictions <- TransferData(anchorset = anchors, refdata = aibs_obj$cluster_label, dims = 1:50)
filtered_cells_obj <- AddMetaData(filtered_cells_obj, metadata = predictions)

saveRDS(filtered_cells_obj, "./seurat_objects/snRNA_cellBender_corrected_QC_filteredCells_addprediction_MappedToAIBS_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
