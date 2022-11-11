# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)

# read in data with 1st round annotations
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_1stRoundAnnotated_SeuratV4_object.rds")

## subset using `rna_anno_1stRound_level_1` in metadata
data_obj_Inh <- subset(x = data_obj_sub, subset = rna_anno_1stRound_level_1 == "Inh_neuron")
data_obj_Exc <- subset(x = data_obj_sub, subset = rna_anno_1stRound_level_1 == "Exc_neuron")
data_obj_NonNeu <- subset(x = data_obj_sub, subset = rna_anno_1stRound_level_1 == "Non_neuron")

# save Seurat object
saveRDS(data_obj_Inh, "./seurat_objects/snRNA_cellBender_corrected_SeuratV4_object_subset_Inh.rds")
saveRDS(data_obj_Exc, "./seurat_objects/snRNA_cellBender_corrected_SeuratV4_object_subset_Exc.rds")
saveRDS(data_obj_NonNeu, "./seurat_objects/snRNA_cellBender_corrected_SeuratV4_object_subset_NonNeu.rds")

## log sessionInfo
sessionInfo()
