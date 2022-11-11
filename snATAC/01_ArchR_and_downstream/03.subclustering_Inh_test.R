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
# # Inh
# meta_data_inh <- meta_data %>%
#   filter(grepl("Inh", predictedGroup_Un))
#
# proj_inh <- subsetArchRProject(
#   ArchRProj = proj,
#   cells = meta_data_inh$cell_id,
#   outputDirectory = "Dracheva_ALS_FTD_addRNA_Inh_subset_test",
#   dropCells = TRUE,
#   force = TRUE
# )

proj_inh <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_Inh_subset_test")

getAvailableMatrices(proj_inh)


proj_inh <- addIterativeLSI(ArchRProj = proj_inh,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI",
                        iterations = 4,
                        clusterParams = list(
                          resolution = c(0.1, 0.2, 0.4),
                          sampleCells = 10000,
                         # maxClusters = 6,
                          n.start = 10
                        ),
                        varFeatures = 15000,
                        dimsToUse = 1:30,
                        force = TRUE)

proj_inh <- addClusters(input = proj_inh,
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

proj_inh <- addUMAP(ArchRProj = proj_inh,
                reducedDims = "IterativeLSI",
                name = "UMAP",
                nNeighbors = 60,
                minDist = 0.01,
                metric = "cosine",
                seed = 666,
                force = TRUE
                )

p1 <- plotEmbedding(ArchRProj = proj_inh, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = proj_inh, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters-InhSubset-test.pdf",
        ArchRProj = proj_inh, addDOC = FALSE, width = 5, height = 5)

proj_inh <- addImputeWeights(proj_inh)

markerGenes  <- c(
  "SNAP25",
  "GAD1", "GAD2", "ADARB2", "VIP", "LAMP5", "CPLX3", "CXCL14", "NDNF", "RELN", "LHX6", "PVALB", "SST", "NPY", "SOX6" # Inh
)

p <- plotEmbedding(
  ArchRProj = proj_inh,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj_inh)
)

plotPDF(plotList = p,
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
        ArchRProj = proj_inh,
        addDOC = FALSE, width = 5, height = 5)


proj_inh <- saveArchRProject(ArchRProj = proj_inh)
