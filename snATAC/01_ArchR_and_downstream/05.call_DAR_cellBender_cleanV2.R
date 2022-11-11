library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 8)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_addPseudobulkRep/")

# pairwise test
# sanity checks to compare between cell types

marker_oligo_vs_astro <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Cluster_sep",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Oligo_MCX_Control",
  bgdGroups = "Astro_MCX_Control",
  maxCells = 500,
  scaleTo = 10^4,
  k = 100, 
  bufferRatio = 0.8,
  binarize = FALSE
)

marker_list <- getMarkers(marker_oligo_vs_astro, cutOff = "FDR <= 0.01 & Log2FC >= 1")
marker_list
marker_list_gr <- getMarkers(marker_oligo_vs_astro, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

pma <- markerPlot(seMarker = marker_oligo_vs_astro, name = "Oligo_MCX_Control", cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1", plotAs = "MA")
pma

pv <- markerPlot(seMarker = marker_oligo_vs_astro, name = "Oligo_MCX_Control", cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv

plotPDF(pma, pv, name = "MCX_Control_Oligo_vs_Astro_Markers_MA_Volcano", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

# plot tracks

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Cluster_sep", 
  # region = marker_list_gr$Oligo_MCX_Control[1],
  geneSymbol = c("TJAP1"),
  features =  marker_list_gr$Oligo_MCX_Control,
  upstream = 50000,
  downstream = 50000
)
# grid::grid.draw(p)
grid::grid.draw(p$TJAP1)
plotPDF(p, name = "Plot-Tracks-With-Features_MCX_Control_Oligo_vs_Astro_DAR_top1", width = 6, height = 8, ArchRProj = proj, addDOC = FALSE)

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

plotPDF(pma, pv, name = "MCX_Astro_ALS_vs_Control_Markers_MA_Volcano", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

marker_list_gr <- getMarkers(marker_ALS_vs_Control_astro, cutOff = "Log2FC >= 2", returnGR = TRUE)

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Cluster_sep", 
  region = marker_list_gr$Astro_MCX_ALS[3],
  # geneSymbol = c("TJAP1"),
  features =  marker_list_gr$Astro_MCX_ALS,
  upstream = 50000,
  downstream = 50000
)
grid::grid.draw(p)
# grid::grid.draw(p$TJAP1)



# todo: test in each cell type
# ...


# find Marker Peaks (one group vs the rest)
marker_all <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Cluster_sep",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 500,
  scaleTo = 10^4,
  k = 100, 
  bufferRatio = 0.8,
  binarize = FALSE
)
saveRDS(marker_all, "test_marker_all.rds")
# 
# heatmap_peaks <- markerHeatmap(
#   seMarker = marker_all,
#   cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
#   transpose = TRUE
# )
# draw(heatmap_peaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotMarkerHeatmap(
  seMarker = marker_all,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
  transpose = TRUE,
  nLabel = 1
)

# get bigwigs
getGroupBW(
  ArchRProj = proj,
  groupBy = "Cluster_sep",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)
