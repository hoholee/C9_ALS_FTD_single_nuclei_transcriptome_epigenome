# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 16)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep")

# subset the ArchR proj, remove ambiguous clusters
index_keep <- grep("_Mixed_", proj$atac_anno_level_2, invert = TRUE)
proj_sub <- proj[index_keep, ]

# find Marker Peaks (one level-2 cell type vs the rest)
marker_all <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "atac_anno_level_2",
  normBy = "ReadsInTSS",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 500,
  scaleTo = 10^4,
  k = 100,
  bufferRatio = 0.8,
  binarize = TRUE
)

saveRDS(marker_all, "snATAC_marker_peaks_by_annoLevel2_clean.rds")

# session info
sessionInfo()