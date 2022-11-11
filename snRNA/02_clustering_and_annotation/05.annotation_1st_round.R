# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)

# read in clustered data, leiden clustering r = 1
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBender_corrected_postQC_scTransformed_clustered_SeuratV4_object.rds")

meta_data <- data_obj_sub@meta.data %>% 
  rownames_to_column("cell_id") %>% 
  as_tibble()

# Dot plot for annotation
DotPlot(data_obj_sub, cols = c("lightgrey", "#C73D4C"),
        features = c("GAD1", "GAD2", "ADARB1", "ADARB2", "SLC17A7",
                     "LHX6", "SOX6", "SNAP25", "MEF2C", "MBP", "MAG", "MOG", "CNP", "MOBP", "PLP1", "OPALIN", "ENPP6",
                     "SLC6A1", "SATB2", "GRIN2B", "RBFOX3", "CUX2", "FOXP2",
                     "TOX", "HTR2C", "PCSK5", "SULF1", "ADRA1A",
                     "THEMIS", "SMYD1", "LINC00299", "SEMA3E", "NPFFR2", "TRABD2A", "GABRG1",
                     "OTOGL", "CCDC68", "ATP7B", "SV2C", "SASH1", "SEMA3A", "CCBE1", "PDGFD",
                     "RELN", "BDNF", "RORB", "PVALB", "SCUBE3", "SST", "LAMP5", "VIP",
                     "TLE4", "PDZRN4", "TSHZ2", "PAX6", "CSF1R", "C1QB", "AQP4", "AQP1", "GRIA1", "PLCG1",
                     "ACTA2", "MYH11", "TBX18", "COLEC12", "VWF", "CYP1B1", "SOX10", "CEMIP", "P2RY14", "SLC4A4",
                     "SPINK6", "BARHL2", "CCDC155", "GABRA6", "FAT2", "PRR35", "CDH15",
                     "CRTAM", "ZIC5", "CBLN3", "SLC22A31",
                     "LDLRAP1", "GPR17", "VCAN", "NEU4", "BMP4", "CRYM", "NEFH", "POU3F1")) +
  RotatedAxis() +
 theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")

# add 1st round level 1 and level 2 annotations
meta_data_annotated <- meta_data %>% 
  mutate(
    rna_anno_1stRound_level_1 = case_when(
      seurat_clusters %in% c(14, 16, 21, 26) ~ "Inh_neuron",
      seurat_clusters %in% c(6, 10, 17, 18, 20, 24, 25, 28) ~ "Exc_neuron",
      seurat_clusters %in% c(2, 15, 19, 11, 1, 3, 4, 5, 8, 9, 12, 13, 23, 27, 7, 22, 29) ~ "Non_neuron",
      ),
    rna_anno_1stRound_level_2 = case_when(
      seurat_clusters %in% c(1, 3, 4, 5, 8, 9, 12, 13, 23, 27) ~ "Oligo",
      seurat_clusters %in% c(2, 15) ~ "Astro",
      seurat_clusters == 7 ~ "OPC",
      seurat_clusters == 11 ~ "Micro",
      seurat_clusters == 19 ~ "Endo",
      seurat_clusters %in% c(22) ~ "VLMC",
      seurat_clusters %in% c(14) ~ "Inh_MGE_SST",
      seurat_clusters %in% c(26) ~ "Inh_MGE_PVALB",
      seurat_clusters %in% c(16) ~ "Inh_CGE_VIP",
      seurat_clusters %in% c(21) ~ "Inh_CGE_LAMP5",
      seurat_clusters %in% c(6) ~ "Exc_superficial",
      seurat_clusters %in% c(10) ~ "Exc_intermediate",
      seurat_clusters %in% c(17, 18, 20, 24, 25, 28) ~ "Exc_deep",
      seurat_clusters == 29 ~ "NK"
      )
  )

data_obj_sub$rna_anno_1stRound_level_1 <- meta_data_annotated$rna_anno_1stRound_level_1
data_obj_sub$rna_anno_1stRound_level_2 <- meta_data_annotated$rna_anno_1stRound_level_2

# save Seurat object
saveRDS(data_obj_sub, file = "./seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_1stRoundAnnotated_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
