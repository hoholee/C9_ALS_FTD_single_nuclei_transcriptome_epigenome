library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 8)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_addPseudobulkRep/")
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_addRNA_addPseudobulkRep_testDAR", load = TRUE)

# add peak sets of interest
# MCX Control Oligo vs. Astro
peak_set1 <- readRDS("./Dracheva_ALS_FTD_addRNA_addPseudobulkRep/PeakCalls/Astro_MCX_ALS-reproduciblePeaks.gr.rds")
peak_set2 <- readRDS("./Dracheva_ALS_FTD_addRNA_addPseudobulkRep/PeakCalls/Astro_MCX_Control-reproduciblePeaks.gr.rds")

# filter peaks by reproducibility
peak_set1 <- peak_set1[peak_set1$Reproducibility >= 6]
peak_set2 <- peak_set2[peak_set2$Reproducibility >= 6]

groupPeaks <- GRangesList(peak_set1, peak_set2)
names(groupPeaks) <- c("Astro_MCX_ALS", "Astro_MCX_Control")
unionPeaks <- unlist(groupPeaks)
unionPeaks <- nonOverlappingGR(unionPeaks, by = "groupScoreQuantile", decreasing = TRUE)

# add peak set
proj <- addPeakSet(ArchRProj = proj, peakSet = unionPeaks, force = TRUE)

# add peak matrix
proj <- addPeakMatrix(
  ArchRProj = proj,
  binarize = FALSE,
  ceiling = 4,
  force = TRUE
)

saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_addRNA_addPseudobulkRep_testDAR", load = FALSE)


# test disease vs control
marker_ALS_vs_Control_astro <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Cluster_sep",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Astro_MCX_ALS",
  bgdGroups = "Astro_MCX_Control",
  maxCells = 500,
  scaleTo = 10^4,
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE
)

pma <- markerPlot(seMarker = marker_ALS_vs_Control_astro, name = "Astro_MCX_ALS", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
pma

pv <- markerPlot(seMarker = marker_ALS_vs_Control_astro, name = "Astro_MCX_ALS", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv

# plotPDF(pma, pv, name = "MCX_Astro_ALS_vs_Control_Markers_MA_Volcano_useAstroPeaksOnly", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
plotPDF(pma, pv, name = "MCX_Astro_ALS_vs_Control_Markers_MA_Volcano_useAstroPeaksOnly_reproducibility6", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

# get peak matrix
peak_mat <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix"
)
saveRDS(peak_mat, "test_peak_matrix_MCX_Astro_ALS_vs_Control_useAstroPeaksOnly_reproducibility6_SEobj.rds")
