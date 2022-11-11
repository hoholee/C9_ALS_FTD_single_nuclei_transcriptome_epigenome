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

# find Marker Peaks (one group vs the rest)
# group defined by snATAC_anno_level_2, region and disease
marker_all <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "group",
  normBy = "ReadsInTSS",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 500,
  scaleTo = 10^4,
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE
)

saveRDS(marker_all, "snATAC_marker_peaks_by_annoLevel2_region_disease.rds")

# session info
sessionInfo()