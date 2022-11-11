# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)

## read the seurat object subset by major cell types in 1st round annotations
clusters <- c("Exc", "Inh", "NonNeu")

run_sctransform <- function(selected_cluster){

  message(paste0("Processing cluster: ", selected_cluster))
  data_obj_sub <- readRDS(paste0("./seurat_objects/snRNA_cellBender_corrected_SeuratV4_object_subset_", selected_cluster, ".rds"))

  ## use the sctransform framework
  data_obj_sub <- SCTransform(data_obj_sub,
                              do.correct.umi = TRUE,
                              variable.features.n = 3000,
                              vars.to.regress = c("percent_mt", "seq_batch"),
                              do.scale = FALSE,
                              do.center = TRUE,
                              seed.use = 666,
                              verbose = TRUE)

  saveRDS(data_obj_sub, file = paste0("./seurat_objects/snRNA_cellBender_corrected_dataset_SeuratV4_object_subset_", selected_cluster, "_scTransform_renormalized.rds"))
}

walk(clusters, run_sctransform)

## log sessionInfo
sessionInfo()
