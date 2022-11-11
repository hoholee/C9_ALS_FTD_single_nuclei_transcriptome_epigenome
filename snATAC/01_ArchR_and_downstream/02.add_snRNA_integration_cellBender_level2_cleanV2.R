# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 1)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_rmLowQualitySample/")

seRNA <- readRDS("/cndd2/junhao/ALS_FTD_singleCell/snRNA_postCellBender/seurat_objects/snRNA_cellBneder_corrected_postQC_scTransformed_clustered_2ndRoundAnnotated_cleanV2_SeuratV4_object.rds")

DefaultAssay(seRNA) <- "RNA"
# Unconstrained Integration

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "rna_anno_2ndRound_level_2",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

pal <- paletteDiscrete(values = seRNA$rna_anno_2ndRound_level_2)


p1 <- plotEmbedding(
  proj,
  colorBy = "cellColData",
  name = "predictedGroup_Un",
  pal = pal
)

plotPDF(p1, name = "Plot-UMAP-RNA-Integration-Full-level2-cleanV2.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# save
saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_addRNA_cellBender_level2_cleanV2", load = FALSE)

## log sessionInfo
sessionInfo()
