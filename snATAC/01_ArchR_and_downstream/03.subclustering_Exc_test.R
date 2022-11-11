library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 1)

addArchRGenome("hg38")

# proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA/")
#
# meta_data <- getCellColData(proj) %>%
#   as.data.frame() %>%
#   rownames_to_column("cell_id") %>%
#   as_tibble()
#
# # Exc
# meta_data_exc <- meta_data %>%
#   filter(grepl("Exc", predictedGroup_Un))
#
# proj_exc <- subsetArchRProject(
#   ArchRProj = proj,
#   cells = meta_data_exc$cell_id,
#   outputDirectory = "Dracheva_ALS_FTD_addRNA_Exc_subset_test",
#   dropCells = TRUE,
#   force = TRUE
# )

proj_exc <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_Exc_subset_test")

getAvailableMatrices(proj_exc)


proj_exc <- addIterativeLSI(ArchRProj = proj_exc,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI",
                        iterations = 4,
                        clusterParams = list(
                          resolution = c(0.4),
                          sampleCells = 10000,
                          maxClusters = 6,
                          n.start = 10
                        ),
                        varFeatures = 15000,
                        dimsToUse = 1:30,
                        force = TRUE)

proj_exc <- addClusters(input = proj_exc,
                    reducedDims = "IterativeLSI",
                    name = "Clusters",
                    seed = 666,
                    method = "Seurat",
                    knnAssign = 10,
                    nOutlier = 5,
                    maxClusters = 50,
                    resolution = 2,
                    force = TRUE
                    )

proj_exc <- addUMAP(ArchRProj = proj_exc,
                reducedDims = "IterativeLSI",
                name = "UMAP",
                nNeighbors = 60,
                minDist = 0.01,
                metric = "cosine",
                seed = 666,
                force = TRUE
                )

p1 <- plotEmbedding(ArchRProj = proj_exc, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = proj_exc, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters-ExcSubset-test.pdf",
        ArchRProj = proj_exc, addDOC = FALSE, width = 5, height = 5)

proj_exc <- addImputeWeights(proj_exc)

markerGenes  <- c(
  "SNAP25", "SLC17A7", "SATB2", "CUX2", "TSHZ2", "RORB", "FOXP2", "TLE4" # Exc
)

p <- plotEmbedding(
  ArchRProj = proj_exc,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj_exc)
)

plotPDF(plotList = p,
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
        ArchRProj = proj_exc,
        addDOC = FALSE, width = 5, height = 5)


proj_exc <- saveArchRProject(ArchRProj = proj_exc)
