library(tidyverse)
library(ArchR)
set.seed(666)

addArchRThreads(threads = 16)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_rmLowQualitySample/")

getAvailableMatrices(proj)

proj <- filterDoublets(ArchRProj = proj)

proj <- addIterativeLSI(ArchRProj = proj,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI",
                        iterations = 2,
                        clusterParams = list(
                          resolution = c(2),
                          sampleCells = 10000,
                          maxClusters = 6,
                          n.start = 10
                        ),
                        varFeatures = 25000,
                        dimsToUse = 1:30)

proj <- addClusters(input = proj,
                    reducedDims = "IterativeLSI",
                    name = "Clusters",
                    seed = 666,
                    method = "Seurat",
                    knnAssign = 10,
                    nOutlier = 5,
                    maxClusters = 50,
                    resolution = 2,
                    # force = TRUE
                    )

proj <- addUMAP(ArchRProj = proj,
                reducedDims = "IterativeLSI",
                name = "UMAP",
                nNeighbors = 60,
                minDist = 0.1,
                metric = "cosine",
                seed = 666,
                # force = TRUE
                )

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

# ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


proj <- addImputeWeights(proj)

markerGenes  <- c(
  "AQP4", "APOE", "AQP1", "GRIA1", "PLCG1",  # Astro
  "C1QB", "P2RY12", # Micro
  "VWF", "PECAM1", "CLDN5", # Endo
  "COLEC12", "TBX18", "CYP1B1", # VLMC
  "ACTA2", "MYH11", "TAGLN", # Pericytes
  "PDGFRA", "VCAN", # OPC
  "MOG", "CNP", "MBP", "MAG", "MOBP", # Oligo 
  "SNAP25", "SLC17A7", "SATB2", "CUX2", "TSHZ2", "RORB", "FOXP2", "TLE4", # Exc
  "ADARB2", "VIP", "LAMP5", "CPLX3", "CXCL14", "NDNF", "RELN", "LHX6", "PVALB", "SST", "NPY", "SOX6" # Inh
)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

# p$SNAP25

#Rearrange for grid plotting
# p2 <- lapply(p, function(x){
#   x + guides(color = FALSE, fill = FALSE) + 
#     theme_ArchR(baseSize = 6.5) +
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#     theme(
#       axis.text.x=element_blank(), 
#       axis.ticks.x=element_blank(), 
#       axis.text.y=element_blank(), 
#       axis.ticks.y=element_blank()
#     )
# })
# do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)

# grid::grid.newpage()
# grid::grid.draw(p$CD14)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

# ArchRBrowser(ArchRProj = proj)

# need to double-check the important notes on saving a ArchRProject
proj <- saveArchRProject(ArchRProj = proj)
