# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(scales)
library(scico)
library(tictoc)
library(shades)

# read in annotated data
# data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_clean_SeuratV4_object.rds")
# meta_data <- data_obj_sub@meta.data

# load the gene x cell CPM matrix
cpm_mat <- readRDS("snRNA_geneByCell_dgCMatrix_RNA_CPM.rds")

meta_data <- read_tsv("metadata_all_cells_2nd_round_annotations.txt") %>% 
  filter(rna_anno_2ndRound_level_3 != "Ambiguous") %>% 
  mutate(group = paste(region, disease, rna_anno_2ndRound_level_2, sep = "__"))

all.equal(colnames(cpm_mat), meta_data$cell_id)

group <- meta_data$group %>% as_factor()
names(group) <- meta_data$cell_id

mm <- model.matrix(~0 + group)
colnames(mm) <- levels(group)

cpm_mat_log2 <- cpm_mat
cpm_mat_log2@x <- log2(cpm_mat_log2@x + 1)

cpm_log2_sum <- cpm_mat_log2 %*% mm

group_count <- meta_data %>% count(group)
count_vec <- group_count$n
names(count_vec) <- group_count$group
count_vec <- count_vec[colnames(cpm_log2_sum)]

all.equal(colnames(cpm_log2_sum), names(count_vec))

cpm_log2_avg <- cpm_log2_sum / count_vec[col(cpm_log2_sum)]

res <- cpm_log2_avg %>%
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  as_tibble()

write_tsv(res, "avg_log2CPM_by_region_disease_level2_celltypes.txt")

## log sessionInfo
sessionInfo()
