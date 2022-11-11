# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)

# read in full dataset with 1st round annotation 
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_1stRoundAnnotated_SeuratV4_object.rds")
meta_data <- data_obj_sub@meta.data %>% rownames_to_column("cell_id") %>% as_tibble()

# read level 3 annotations from sub-clustering results
exc_anno <- read_tsv("meta_annotated_Exc_sub_clusters.txt") %>% select(cell_id, rna_anno_level_3)
inh_anno <- read_tsv("meta_annotated_Inh_sub_clusters.txt") %>% select(cell_id, rna_anno_level_3)
nonneu_anno <- read_tsv("meta_annotated_NonNeu_sub_clusters.txt") %>% select(cell_id, rna_anno_level_3)

df_anno <- bind_rows(exc_anno, inh_anno, nonneu_anno)

meta_data_mod <- meta_data %>% 
  left_join(df_anno, by = "cell_id") %>% 
  rename(rna_anno_1stRound_level_3 = rna_anno_level_3)

data_obj_sub$rna_anno_1stRound_level_3 <- meta_data_mod$rna_anno_1stRound_level_3

meta_data_final <- meta_data_mod %>% 
  mutate(
    rna_anno_2ndRound_level_3 = case_when(
      grepl("DB_", rna_anno_1stRound_level_3) ~ "Ambiguous",
      TRUE ~ rna_anno_1stRound_level_3
    ),
    rna_anno_2ndRound_level_2 = case_when(
      grepl("^NK", rna_anno_2ndRound_level_3) ~ "NK_cell",
      grepl("Astro_", rna_anno_2ndRound_level_3) ~ "Astro",
      grepl("Endo", rna_anno_2ndRound_level_3) ~ "Endo",
      grepl("Micro", rna_anno_2ndRound_level_3) ~ "Micro",
      grepl("OPC", rna_anno_2ndRound_level_3) ~ "OPC",
      grepl("Oligo", rna_anno_2ndRound_level_3) ~ "Oligo",
      grepl("VLMC", rna_anno_2ndRound_level_3) ~ "VLMC",
      grepl("Inh_VIP", rna_anno_2ndRound_level_3) ~ "Inh_VIP",
      grepl("Inh_SST", rna_anno_2ndRound_level_3) ~ "Inh_SST",
      grepl("Inh_PVALB", rna_anno_2ndRound_level_3) ~ "Inh_PVALB",
      grepl("Inh_LAMP5", rna_anno_2ndRound_level_3) ~ "Inh_LAMP5",
      grepl("Inh_ADARB2", rna_anno_2ndRound_level_3) ~ "Inh_ADARB2_Other",
      grepl("Exc_L6", rna_anno_2ndRound_level_3) ~ "Exc_deep",
      grepl("Exc_L5", rna_anno_2ndRound_level_3) ~ "Exc_deep",
      grepl("Exc_L3", rna_anno_2ndRound_level_3) ~ "Exc_intermediate",
      grepl("Exc_L2-3", rna_anno_2ndRound_level_3) ~ "Exc_intermediate",
      grepl("Exc_L2_", rna_anno_2ndRound_level_3) ~ "Exc_superficial",
      grepl("Exc_TUBB2A", rna_anno_2ndRound_level_3) ~ "Exc_unknown",
      grepl("Ambiguous", rna_anno_2ndRound_level_3) ~ "Ambiguous"
    ),
    rna_anno_2ndRound_level_1 = case_when(
      grepl("Exc_", rna_anno_2ndRound_level_2) ~ "Exc_neuron",
      grepl("Inh_", rna_anno_2ndRound_level_2) ~ "Inh_neuron",
      grepl("Ambiguous", rna_anno_2ndRound_level_2) ~ "Ambiguous",
      TRUE ~ "Non_neuron"
    )
  )

data_obj_sub$rna_anno_2ndRound_level_1 <- meta_data_final$rna_anno_2ndRound_level_1
data_obj_sub$rna_anno_2ndRound_level_2 <- meta_data_final$rna_anno_2ndRound_level_2
data_obj_sub$rna_anno_2ndRound_level_3 <- meta_data_final$rna_anno_2ndRound_level_3

Idents(data_obj_sub) <- "rna_anno_2ndRound_level_3"

# save Seurat object
saveRDS(data_obj_sub, file = "./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_SeuratV4_object.rds")

# exclude NK cells and ambiguous cells and save to another Seurat object
data_obj_sub2 <- subset(data_obj_sub, subset = rna_anno_2ndRound_level_3 %in% c("NK_cell", "Ambiguous"), invert = TRUE)
saveRDS(data_obj_sub2, file = "./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_clean_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
