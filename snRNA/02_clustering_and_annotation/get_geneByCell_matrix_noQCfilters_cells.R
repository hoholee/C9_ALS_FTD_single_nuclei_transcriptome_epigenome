# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)

data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBender_corrected_noQCFilters_scTransformed_clustered_SeuratV4_object.rds")

meta_data <- data_obj_sub@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

# get the raw counts from the `RNA` slot in the Seurat object
mat <- data_obj_sub@assays$RNA@counts
all.equal(unname(colSums(mat)), meta_data$nCount_RNA)

# get the SCTransformed "raw" counts from the `SCT` slot in the Seurat object
mat_SCT <- data_obj_sub@assays$SCT@counts
all.equal(unname(colSums(mat_SCT)), meta_data$nCount_SCT)

# compute CPM by individual cell
mat_cell_cpm <- mat
mat_cell_cpm@x <- 10^6 * mat_cell_cpm@x / rep.int(colSums(mat_cell_cpm), diff(mat_cell_cpm@p))

mat_SCT_cell_cpm <- mat_SCT
mat_SCT_cell_cpm@x <- 10^6 * mat_SCT_cell_cpm@x / rep.int(colSums(mat_SCT_cell_cpm), diff(mat_SCT_cell_cpm@p))

# save output as rds object
saveRDS(mat, "snRNA_noQCfilters_geneByCell_dgCMatrix_RNA_raw_count.rds")
saveRDS(mat_cell_cpm, "snRNA_noQCfilters_geneByCell_dgCMatrix_RNA_CPM.rds")

saveRDS(mat_SCT, "snRNA_noQCfilters_geneByCell_dgCMatrix_SCTransformed_raw_count.rds")
saveRDS(mat_SCT_cell_cpm, "snRNA_noQCfilters_geneByCell_dgCMatrix_SCTransformed_CPM.rds")

## log sessionInfo
sessionInfo()
