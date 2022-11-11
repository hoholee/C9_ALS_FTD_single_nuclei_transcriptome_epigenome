library(tidyverse)
library(ArchR)
set.seed(666)

addArchRThreads(threads = 20)
addArchRGenome("hg38")

# setup input files
inputFiles <- list.files("./atac_fragments", pattern = ".gz$", full.names = TRUE)
names(inputFiles) <- inputFiles %>% str_replace("_atac_fragments.tsv.gz", "") %>% str_replace("./atac_fragments/", "")

# remove subject 908 and low quality samples
rm_idx <- which(names(inputFiles) %in% c("MCX_FTD_908", "MCX_FTD_36", "MCX_FTD_54", "MCX_FTD_61",
                                         "mFCX_Control_945", "mFCX_FTD_36", "mFCX_FTD_61", "mFCX_FTD_908"))
inputFiles_clean <- inputFiles[-rm_idx]

# consider to create a custom ArchRGenome with the gene annotation used in the snRNA analysis
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles_clean,
  sampleNames = names(inputFiles_clean),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Dracheva_ALS_FTD_rmLowQualitySample",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)
# save ArchRProject to a backup path before further processing
saveArchRProject(ArchRProj = proj, outputDirectory = "Dracheva_ALS_FTD_rmLowQualitySample", load = FALSE)
