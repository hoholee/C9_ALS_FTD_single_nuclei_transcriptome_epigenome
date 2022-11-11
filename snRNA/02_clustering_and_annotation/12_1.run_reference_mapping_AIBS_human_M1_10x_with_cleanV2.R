# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(pheatmap)
library(viridis)
library(scico)

data_obj <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_cleanV2_SeuratV4_object.rds")
aibs_obj <- readRDS("/cndd2/junhao/AIBS_human_M1_10x/AIBS_human_M1_10x_Seurat_obj_scTransformed_PCA50_clustered.rds")

anchors <- FindTransferAnchors(reference = aibs_obj,
                               query = data_obj,
                               normalization.method = "SCT",
                               query.assay = "RNA",
                               dims = 1:50)
predictions <- TransferData(anchorset = anchors, refdata = aibs_obj$cluster_label, dims = 1:50)
data_obj <- AddMetaData(data_obj, metadata = predictions)

saveRDS(data_obj, "./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_cleanV2_addAIBSprediction_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
