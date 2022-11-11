library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 1)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA/")

meta_data <- getCellColData(proj) %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

# Neunons
meta_data_neu <- meta_data %>%
  filter(grepl("Exc", predictedGroup_Un) | grepl("Inh", predictedGroup_Un))

proj_neu <- subsetArchRProject(
  ArchRProj = proj,
  cells = meta_data_neu$cell_id,
  outputDirectory = "Dracheva_ALS_FTD_addRNA_Neu_subset_test",
  dropCells = TRUE,
  force = TRUE
)

getAvailableMatrices(proj_neu)


proj_neu <- addIterativeLSI(ArchRProj = proj_neu,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI",
                        iterations = 4,
                        clusterParams = list(
                          resolution = c(0.1, 0.2, 0.4),
                          sampleCells = 10000,
                          maxClusters = 6,
                          n.start = 10
                        ),
                        varFeatures = 15000,
                        dimsToUse = 1:30,
                        force = TRUE)

proj_neu <- addClusters(input = proj_neu,
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

proj_neu <- addUMAP(ArchRProj = proj_neu,
                reducedDims = "IterativeLSI",
                name = "UMAP",
                nNeighbors = 60,
                minDist = 0.01,
                metric = "cosine",
                seed = 666,
                force = TRUE
                )

p1 <- plotEmbedding(ArchRProj = proj_neu, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = proj_neu, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters-ExcSubset-test.pdf",
        ArchRProj = proj_neu, addDOC = FALSE, width = 5, height = 5)

proj_neu <- addImputeWeights(proj_neu)

markerGenes  <- c(
  "SNAP25", "SLC17A7", "SATB2", "CUX2", "TSHZ2", "RORB", "FOXP2", "TLE4", # Exc
  "GAD1", "GAD2", "ADARB2", "VIP", "LAMP5", "CPLX3", "CXCL14", "NDNF", "RELN", "LHX6", "PVALB", "SST", "NPY", "SOX6" # Inh
)

p <- plotEmbedding(
  ArchRProj = proj_neu,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj_neu)
)

plotPDF(plotList = p,
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
        ArchRProj = proj_neu,
        addDOC = FALSE, width = 5, height = 5)


proj_neu <- saveArchRProject(ArchRProj = proj_neu)
