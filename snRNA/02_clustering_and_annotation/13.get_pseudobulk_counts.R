# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)

data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_clean_SeuratV4_object.rds")

meta_data <- data_obj_sub@meta.data %>% 
  rownames_to_column("cell_id") %>%
  as_tibble() %>% 
  mutate(group_level2 = paste(region, disease, rna_anno_2ndRound_level_2, sep = " "),
         group_level3 = paste(region, disease, rna_anno_2ndRound_level_3, sep = " "))

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
saveRDS(mat, "snRNA_geneByCell_dgCMatrix_RNA_raw_count.rds")
saveRDS(mat_cell_cpm, "snRNA_geneByCell_dgCMatrix_RNA_CPM.rds")

saveRDS(mat_SCT, "snRNA_geneByCell_dgCMatrix_SCTransformed_raw_count.rds")
saveRDS(mat_SCT_cell_cpm, "snRNA_geneByCell_dgCMatrix_SCTransformed_CPM.rds")

# todo: fix the issue of differences in gene annotation!
# compute TPM by individual cell
# read in bed file to get gene length (gene body, not transcript length in this case)
gene_anno <- read_tsv("Homo_sapiens.GRCh38.93.filtered.bed",
                      col_names = c("chr", "start", "end", "id", "score", "strand"),
                      col_types = "ciicic") %>% 
  mutate(id_bak = id) %>% 
  separate(id_bak, c("gene_id", "gene_name_unique"), sep = "_") %>% 
  mutate(gene_length = end - start)
# read a lookup table to handle the issue of repeat gene names and the ".x" issue
gene_id_lookup <- read_tsv("/cndd2/junhao/ALS_FTD_singleCell/integration/latest/data/raw/snrna_raw.gene")

# save the order of gene in the matrix
mat_cell_gene_order <- tibble(gene_name = rownames(mat)) %>% 
  rownames_to_column("mat_idx") %>% 
  mutate(mat_idx = as.numeric(mat_idx))

gene_anno_cell_mod <- gene_anno %>% 
  left_join(gene_id_lookup, by = "id") %>% 
  left_join(mat_cell_gene_order, by = "gene_name") %>% 
  filter(!is.na(mat_idx)) %>% 
  arrange(mat_idx)

mat_cell_tpm <- mat / gene_anno_cell_mod$gene_length * 1000
mat_cell_tpm@x <- 10^6 * mat_cell_tpm@x / rep.int(colSums(mat_cell_tpm), diff(mat_cell_tpm@p))

saveRDS(mat_cell_tpm, "snRNA_geneByCell_dgCMatrix_SCTransformed_TPM.rds")
saveRDS(mat_cell_tpm, "snRNA_geneByCell_dgCMatrix_SCTransformed_TPM.rds")

# create pseudobulk counts
# pool level-2, split by brain region, cell type and disease
# todo: further split by individual
# get pooled library size
pool_lib_size <- meta_data %>%
  group_by(group_level2) %>%
  summarise(pooled_lib_size = sum(nCount_SCT)) %>% 
  ungroup()
pool_lib_size_lookup <- pool_lib_size$pooled_lib_size
names(pool_lib_size_lookup) <- pool_lib_size$group_level2 

# get pooled counts (pseudobulk)
# make model matrix
mm <- model.matrix(~ 0 + meta_data$group_level2)
colnames(mm) <- colnames(mm) %>%
  str_replace("meta_data\\$group_level2", "")

## multiple row of mat (gene) onto column of the model matrix (cell-type annotation)
## `mat_summary_mm` will be the total counts of each gene by cluster, by region, by disease
## convert the output to dgCMatrix format
tic()
mat_summary_mm <- mat %*% mm
mat_summary_mm <- as(mat_summary_mm, "dgCMatrix")
toc()

# compute CPM
mat_cpm <- mat_summary_mm
all.equal(colSums(mat_cpm), pool_lib_size_lookup)
mat_cpm@x <- 10^6 * mat_cpm@x / rep.int(colSums(mat_cpm), diff(mat_cpm@p))

# compute TPM
# save the order of gene in the matrix
mat_gene_order <- tibble(gene_name = rownames(mat_summary_mm)) %>% 
  rownames_to_column("mat_idx") %>% 
  mutate(mat_idx = as.numeric(mat_idx))

gene_anno_mod <- gene_anno %>% 
  left_join(gene_id_lookup, by = "id") %>% 
  left_join(mat_gene_order, by = "gene_name") %>% 
  filter(!is.na(mat_idx)) %>% 
  arrange(mat_idx)

mat_tpm <- mat_summary_mm / gene_anno_mod$gene_length * 1000
mat_tpm@x <- 10^6 * mat_tpm@x / rep.int(colSums(mat_tpm), diff(mat_tpm@p))

# save output as rds object
saveRDS(mat_summary_mm, "snRNA_pseudobulk_by_region_disease_annoLevel2_dgCMatrix_raw_count.rds")
saveRDS(mat_cpm, "snRNA_pseudobulk_by_region_disease_annoLevel2_dgCMatrix_CPM.rds")
saveRDS(mat_tpm, "snRNA_pseudobulk_by_region_disease_annoLevel2_dgCMatrix_TPM.rds")

# as plain text
mat_count_df <- mat_summary_mm %>% as_tibble()
mat_cpm_df <- mat_cpm %>% as_tibble()
mat_tpm_df <- mat_tpm %>% as_tibble()

res_cout <- gene_anno_mod %>% bind_cols(mat_count_df)
res_cpm <- gene_anno_mod %>% bind_cols(mat_cpm_df)
res_tpm <- gene_anno_mod %>% bind_cols(mat_tpm_df)

write_tsv(res_cout, "snRNA_pseudobulk_by_region_disease_annoLevel2_raw_count.txt")
write_tsv(res_cpm, "snRNA_pseudobulk_by_region_disease_annoLevel2_CPM.txt")
write_tsv(res_tpm, "snRNA_pseudobulk_by_region_disease_annoLevel2_TPM.txt")

## log sessionInfo
sessionInfo()
