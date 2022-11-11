# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)

## read the post-QC seurat object
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBender_corrected_postQC_SeuratV4_object.rds")

## use the sctransform framework
data_obj_sub <- SCTransform(data_obj_sub,
                            do.correct.umi = TRUE,
                            variable.features.n = 3000,
                            vars.to.regress = c("percent_mt", "seq_batch"),
                            do.scale = FALSE,
                            do.center = TRUE,
                            seed.use = 666,
                            verbose = TRUE)

# save Seurat object
saveRDS(data_obj_sub, file = "./seurat_objects/snRNA_cellBender_postQC_scTransformed_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
