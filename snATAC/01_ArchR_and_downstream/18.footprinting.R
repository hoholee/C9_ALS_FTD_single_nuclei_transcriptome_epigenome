# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(viridis)
library(scico)

set.seed(666)
addArchRThreads(threads = 16)
addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep")

motifPositions <- getPositions(proj)

selected_motifs <- c(
    "AC0242|FOSL/JUND|bZIP",
    "AC0239|BACH/NFE|bZIP",
    "AC0244|MAFG/MAFF|bZIP"
)

seFoot <- getFootprints(
    ArchRProj = proj,
    positions = motifPositions[selected_motifs],
    groupBy = "group",
    useGroups = c("Astro_MCX_ALS", "Astro_MCX_Control"),
    flank = 250,
    minCells = 25
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "Subtract",
    plotName = "Footprints-Subtract-Bias",
    addDOC = FALSE,
    smoothWindow = 5
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "Divide",
    plotName = "Footprints-Divide-Bias",
    addDOC = FALSE,
    smoothWindow = 5
)

selected_motifs <- c(
    "AC0092|ZNF|C2H2_ZF",
    "AC0606|SNAI/ZEB|C2H2_ZF",
    "AC0610|TCF/ASCL|bHLH"
)

seFoot <- getFootprints(
    ArchRProj = proj,
    positions = motifPositions[selected_motifs],
    groupBy = "group",
    useGroups = c("Exc_superficial_MCX_ALS", "Exc_superficial_MCX_Control"),
    flank = 250,
    minCells = 25
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "Subtract",
    plotName = "Footprints-Subtract-Bias-Exc-superficial",
    addDOC = FALSE,
    smoothWindow = 5
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "Divide",
    plotName = "Footprints-Divide-Bias-Exc-superficial",
    addDOC = FALSE,
    smoothWindow = 5
)

selected_motifs <- c(
    "AC0261|SOX/SRY|Sox",
    "AC0227|SPI/BCL11A|Ets"
)

seFoot <- getFootprints(
    ArchRProj = proj,
    positions = motifPositions[selected_motifs],
    groupBy = "group",
    useGroups = c("Oligo_MCX_ALS", "Oligo_MCX_Control", "Micro_MCX_ALS", "Micro_MCX_Control"),
    flank = 250,
    minCells = 25
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "Subtract",
    plotName = "Footprints-Subtract-Bias-celltype",
    addDOC = FALSE,
    smoothWindow = 5
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "Divide",
    plotName = "Footprints-Divide-Bias-celltype",
    addDOC = FALSE,
    smoothWindow = 5
)


selected_motifs <- c(
    "AC0242|FOSL/JUND|bZIP"
)

seFoot <- getFootprints(
    ArchRProj = proj,
    positions = motifPositions[selected_motifs],
    groupBy = "group",
    useGroups = c("OPC_MCX_ALS", "OPC_MCX_Control"),
    flank = 250,
    minCells = 25
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "Subtract",
    plotName = "Footprints-Subtract-Bias-OPC",
    addDOC = FALSE,
    smoothWindow = 5
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "Divide",
    plotName = "Footprints-Divide-Bias-OPC",
    addDOC = FALSE,
    smoothWindow = 5
)
# seGroupMotif <- getGroupSE(
#     ArchRProj = proj,
#     useMatrix = "vierstra_motif_archetypeMatrix",
#     groupBy = "group"
# )

# seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames == "z", ]

# rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
#   rowMaxs(assay(seZ) - assay(seZ)[,x])
# }) %>% Reduce("cbind", .) %>% rowMaxs

# corGSM_MM <- correlateMatrices(
#     ArchRProj = proj,
#     useMatrix1 = "GeneScoreMatrix",
#     useMatrix2 = "vierstra_motif_archetypeMatrix",
#     reducedDims = "IterativeLSI"
# )

# session info
sessionInfo()