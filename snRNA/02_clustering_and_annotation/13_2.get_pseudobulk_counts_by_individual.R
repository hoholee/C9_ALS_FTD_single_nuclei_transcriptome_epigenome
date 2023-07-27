# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(tictoc)

# read in annotated data
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_clean_SeuratV4_object.rds")

meta_data <- data_obj_sub@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble() %>%
  mutate(
    group_level2 = paste(region, disease, subject, rna_anno_2ndRound_level_2, sep = " "),
    group_level3 = paste(region, disease, subject, rna_anno_2ndRound_level_3, sep = " ")
  )

# get the raw counts from the `RNA` slot in the Seurat object
mat <- data_obj_sub@assays$RNA@counts
all.equal(unname(colSums(mat)), meta_data$nCount_RNA)

# get the SCTransformed "raw" counts from the `SCT` slot in the Seurat object
mat_SCT <- data_obj_sub@assays$SCT@counts
all.equal(unname(colSums(mat_SCT)), meta_data$nCount_SCT)

# create pseudobulk counts
# pool level-2, split by brain region, cell type and disease, and indiviudal
# get pooled library size for raw counts and SCT adjusted counts
pool_lib_size_RNA <- meta_data %>%
  group_by(group_level2) %>%
  summarise(pooled_lib_size = sum(nCount_RNA)) %>%
  ungroup()
pool_lib_size_lookup_RNA <- pool_lib_size_RNA$pooled_lib_size
names(pool_lib_size_lookup_RNA) <- pool_lib_size_RNA$group_level2

pool_lib_size_SCT <- meta_data %>%
  group_by(group_level2) %>%
  summarise(pooled_lib_size = sum(nCount_SCT)) %>%
  ungroup()
pool_lib_size_lookup_SCT <- pool_lib_size_SCT$pooled_lib_size
names(pool_lib_size_lookup_SCT) <- pool_lib_size_SCT$group_level2

# get pooled counts (pseudobulk)
# make model matrix
mm <- model.matrix(~ 0 + meta_data$group_level2)
colnames(mm) <- colnames(mm) %>%
  str_replace("meta_data\\$group_level2", "")

## multiple row of mat (gene) onto column of the model matrix (cell-type annotation)
## `mat_summary_mm` will be the total counts of each gene by cluster, by region, by disease and by individual
## convert the output to dgCMatrix format
tic()
mat_summary_mm <- mat %*% mm
mat_summary_mm <- as(mat_summary_mm, "dgCMatrix")
toc()

tic()
mat_summary_mm_SCT <- mat_SCT %*% mm
mat_summary_mm_SCT <- as(mat_summary_mm_SCT, "dgCMatrix")
toc()

# compute CPM
mat_cpm <- mat_summary_mm
all.equal(sort(colSums(mat_cpm)), sort(pool_lib_size_lookup_RNA))
mat_cpm@x <- 10^6 * mat_cpm@x / rep.int(colSums(mat_cpm), diff(mat_cpm@p))

mat_cpm_SCT <- mat_summary_mm_SCT
all.equal(sort(colSums(mat_cpm_SCT)), sort(pool_lib_size_lookup_SCT))
mat_cpm_SCT@x <- 10^6 * mat_cpm_SCT@x / rep.int(colSums(mat_cpm_SCT), diff(mat_cpm_SCT@p))

# save output as rds object
saveRDS(mat_summary_mm, "snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_RNA_raw_count.rds")
saveRDS(mat_cpm, "snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_RNA_CPM.rds")
saveRDS(mat_summary_mm_SCT, "snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_SCT_raw_count.rds")
saveRDS(mat_cpm_SCT, "snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_SCT_CPM.rds")

## log sessionInfo
sessionInfo()
