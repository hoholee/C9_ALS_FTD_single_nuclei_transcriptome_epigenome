# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)

data_obj <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_cleanV2_SeuratV4_object.rds")
filtered_cells_obj <- readRDS("./seurat_objects/snRNA_cellBender_corrected_QC_filteredCells_scTransformed_SeuratV4_object.rds")

anchors <- FindTransferAnchors(
    reference = data_obj,
    query = filtered_cells_obj,
    normalization.method = "SCT",
    query.assay = "RNA",
    dims = 1:50
)
predictions <- TransferData(anchorset = anchors, refdata = data_obj$rna_anno_2ndRound_level_2, dims = 1:50)
filtered_cells_obj <- AddMetaData(filtered_cells_obj, metadata = predictions)

saveRDS(filtered_cells_obj, "./seurat_objects/snRNA_cellBender_corrected_QC_filteredCells_addprediction_MappedToManuscriptClusters_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
