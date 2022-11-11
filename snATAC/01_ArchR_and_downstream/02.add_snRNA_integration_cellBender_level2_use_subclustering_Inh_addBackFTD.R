# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 1)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "./Dracheva_ALS_FTD_addRNA_Neu_Inh_subset")

seRNA <- readRDS("/cndd2/junhao/ALS_FTD_singleCell/snRNA_postCellBender/seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_cleanV2_SeuratV4_object.rds")
seRNA_sub <- subset(seRNA, subset = rna_anno_2ndRound_level_1 %in% c("Inh_neuron"))

rm(seRNA)

DefaultAssay(seRNA_sub) <- "RNA"

# Unconstrained Integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_level2",
  reducedDims = "IterativeLSI",
  seRNA = seRNA_sub,
  addToArrow = FALSE,
  groupRNA = "rna_anno_2ndRound_level_2",
  nameCell = "predictedCell_Un_level2",
  nameGroup = "predictedGroup_Un_level2",
  nameScore = "predictedScore_Un_level2"
)

pal <- paletteDiscrete(values = seRNA_sub$rna_anno_2ndRound_level_2)


p1 <- plotEmbedding(
  proj,
  colorBy = "cellColData",
  name = "predictedGroup_Un_level2",
  pal = pal
)

plotPDF(p1, name = "Plot-UMAP-RNA-Integration-Inh-level2.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# save
saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_addRNA_Neu_Inh_subset", load = FALSE)

## log sessionInfo
sessionInfo()
