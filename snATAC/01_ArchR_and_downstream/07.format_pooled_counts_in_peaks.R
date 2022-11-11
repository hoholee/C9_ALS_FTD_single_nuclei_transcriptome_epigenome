library(tidyverse)
library(ArchR)
library(Seurat)
library(tictoc)
set.seed(666)

addArchRThreads(threads = 8)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_addPseudobulkRep_testDAR/")

peak_se <- readRDS("test_peak_matrix_MCX_Astro_ALS_vs_Control_useAstroPeaksOnly_reproducibility6_SEobj.rds")

# extract count matrix
peak_mat <- assay(peak_se)

# extract peaks
peak_df <- peak_se@rowRanges %>% 
  as_tibble() %>% 
  mutate(peak_id = paste(seqnames, start, end, sep = "_"))

# add peak name
rownames(peak_mat) <- peak_df$peak_id

# extract cell metadata
meta_data <- colData(peak_se) %>%
  as_tibble() %>% 
  mutate(sample2 = Sample) %>% 
  separate(sample2, c("brain_region", "disease", "subject")) %>% 
  mutate(group = paste0(Cluster_sep, "_subject", subject))

## make model matrix
mm <- model.matrix(~ 0 + meta_data$group)
colnames(mm) <- colnames(mm) %>% 
  str_replace("meta_data\\$group", "")

## multiple row of mat (gene) onto column of the model matrix (cell-type annotation)
## `mat_summary_mm` will be the total counts of each gene by cluster, by region, by disease, by individual
tic()
mat_summary_mm <- peak_mat %*% mm
toc()

## save to file
count_mat_out <- mat_summary_mm %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("peak") %>% 
  as_tibble()

write_tsv(count_mat_out, "peak_Astro_ALS_vs_Control_atac_count_by_major_clusters_by_region_by_disease_by_individual.txt")

# save pooled library size (nFrags*2)
pool_lib_size <- meta_data %>%
  group_by(group) %>%
  summarise(pooled_lib_size = sum(nFrags)) %>% 
  mutate(pooled_lib_size = pooled_lib_size * 2)
write_tsv(pool_lib_size, "peak_Astro_ALS_vs_Control_pooled_libSize_by_major_clusters_by_region_by_disease_by_individual.txt")

