# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 16)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_cellBender_level2_cleanV2/")

# read meta data with annotations
meta_data_anno <- read_tsv("./metadata_merged_addAnno.txt")

# add group by level 2 annotations, brain region and disease
# consider doing the same thing with level 3 annotations?
meta_data_add_group <- meta_data_anno %>%
  mutate(sample = Sample) %>%
  separate(sample, c("region", "disease", "subject"), sep = "_") %>%
  mutate(group = paste(atac_anno_level_2, region, disease, sep = "_"))

# retrieve meta data from ArchR object
meta_data <- getCellColData(proj) %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

# check if cell IDs and cluster IDs (from full data set clustering)
# match between the two meta data tables
message("Checking if cell IDs and clusters IDs match between two metadata tables")
all.equal(meta_data$cell_id, meta_data_add_group$cell_id)
# equivelent to the above line
all.equal(getCellNames(proj), meta_data_add_group$cell_id)
all.equal(meta_data$Clusters, meta_data_add_group$Clusters_full)

# update meta data in the ArchR object
# drop previous annotations
proj@cellColData$Clusters <- NULL
proj@cellColData$predictedCell_Un <- NULL
proj@cellColData$predictedGroup_Un <- NULL
proj@cellColData$predictedScore_Un <- NULL

# add new meta data
to_add_columns <- colnames(meta_data_add_group)[!colnames(meta_data_add_group) %in% colnames(meta_data)]

for (i in seq_along(to_add_columns)) {
  proj <- addCellColData(
    ArchRProj = proj,
    data = pull(meta_data_add_group, to_add_columns[i]),
    name = to_add_columns[i],
    cells = getCellNames(proj),
    force = TRUE
  )
}

# add pseudo-bulk replicates
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "group",
  minCells = 50,
  maxCells = 500,
  minReplicates = 6,
  maxReplicates = 28,
  sampleRatio = 0.8
)

# save
saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep", load = FALSE)

# call peaks
proj <- loadArchRProject(path = "./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep")
addArchRThreads(threads = 16)

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "group",
  peakMethod = "Macs2",
  pathToMacs2 = "/cndd/junhao/anaconda3/envs/macs2_latest/bin/macs2", # use a conda version of macs2 that works
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 25,
  excludeChr = c("chrM", "chrY"),
  genomeSize = 2.7e9,
  shift = -75,
  extsize = 150,
  cutOff = 0.1,
  additionalParams = "--nomodel --nolambda",
  extendSummits = 250,
  promoterRegion = c(2000, 100)
)

saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep", load = FALSE)

# add peak matrix
proj <- addPeakMatrix(
  ArchRProj = proj,
  binarize = FALSE,
  ceiling = 4
)

saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep", load = FALSE)

# session info
sessionInfo()