# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)

# read in full dataset with 1st round annotation
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_SeuratV4_object.rds")

# exclude NK cells, Exc_TUBB2A and ambiguous cells and save to another Seurat object
data_obj_sub2 <- subset(data_obj_sub, subset = rna_anno_2ndRound_level_3 %in% c("NK", "Exc_TUBB2A", "Ambiguous"), invert = TRUE)
saveRDS(data_obj_sub2, file = "./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_cleanV2_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
