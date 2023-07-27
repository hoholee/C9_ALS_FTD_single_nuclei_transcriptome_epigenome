# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)

data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_cleanV2_SeuratV4_object.rds")

meta_data <- data_obj_sub@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

# read clean metadata for manuscript
# the cell order are the same as in the Seurat object
meta_data_clean <- read_tsv("./metadata_all_cells_clean_for_manuscript.txt")

# update metadata in the Seurat object, mask Donor ID
data_obj_sub$rna_anno_1stRound_level_1 <- NULL
data_obj_sub$rna_anno_1stRound_level_2 <- NULL
data_obj_sub$rna_anno_1stRound_level_3 <- NULL
data_obj_sub$rna_anno_2ndRound_level_1 <- NULL
data_obj_sub$rna_anno_2ndRound_level_2 <- NULL
data_obj_sub$rna_anno_2ndRound_level_3 <- NULL

data_obj_sub$annotation_cell_class <- meta_data_clean$annotation_cell_class
data_obj_sub$annotation_major_cell_type <- meta_data_clean$annotation_major_cell_type
data_obj_sub$annotation_cell_subtype <- meta_data_clean$annotation_cell_subtype

data_obj_sub$region <- meta_data_clean$brain_region
data_obj_sub$disease <- meta_data_clean$diagnosis
data_obj_sub$subject <- meta_data_clean$donor_id

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
saveRDS(mat, "snRNA_geneByCell_dgCMatrix_RNA_raw_count_clean_for_manuscript.rds")
saveRDS(mat_cell_cpm, "snRNA_geneByCell_dgCMatrix_RNA_CPM_clean_for_manuscript.rds")

saveRDS(mat_SCT, "snRNA_geneByCell_dgCMatrix_SCTransformed_raw_count_clean_for_manuscript.rds")
saveRDS(mat_SCT_cell_cpm, "snRNA_geneByCell_dgCMatrix_SCTransformed_CPM_clean_for_manuscript.rds")

# save Seurat object
saveRDS(data_obj_sub, file = "./seurat_objects/snRNA_SeuratV4_object_for_manuscript.rds")

## log sessionInfo
sessionInfo()
