# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 8)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep")

meta <- getCellColData(proj) %>%
  as_tibble() %>%
  mutate(group_individual = paste(atac_anno_level_2, region, disease, subject, sep = "_"))

proj <- addCellColData(
  proj,
  data = meta$group_individual,
  name = "group_individual",
  cells = getCellNames(proj)
)

# get bigwigs
getGroupBW(
  ArchRProj = proj,
  groupBy = "group_individual",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

# session info
sessionInfo()
