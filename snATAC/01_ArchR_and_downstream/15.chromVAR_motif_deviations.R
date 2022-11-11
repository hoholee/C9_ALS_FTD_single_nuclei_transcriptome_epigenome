# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(666)

addArchRThreads(threads = 16)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep")

# add back impute weights
proj <- addImputeWeights(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    corCutOff = 0.75,
    dimsToUse = NULL,
    scaleDims = NULL,
    td = 3,
    ka = 4,
    sampleCells = 5000,
    nRep = 2,
    k = 15,
    epsilon = 1,
    useHdf5 = TRUE,
    randomSuffix = FALSE,
    threads = getArchRThreads(),
    seed = 666,
    verbose = TRUE,
    logFile = createLogFile("addImputeWeights")
)

# check if motif annotations were added in the ArchR object
proj@peakAnnotation %>% names()

# add backgroud peaks
# sample peaks based on similarity in GC-content and number of fragments across all samples
# using the Mahalanobis distance
proj <- addBgdPeaks(
    ArchRProj = proj,
    nIterations = 50,
    w = 0.1,
    binSize = 50,
    seed = 666,
    method = "chromVAR",
    force = TRUE
)

# compute per-cell deviations accross all of our motif annotations
proj <- addDeviationsMatrix(
    ArchRProj = proj,
    peakAnnotation = "vierstra_motif_archetype",
    out = c("z", "deviations"),
    binarize = FALSE,
    force = TRUE
)

# plot variability of the deviations
# which is the standard deviation of the z scores computed above across all cell/samples for a set of peaks
plot_var_dev <- getVarDeviations(
    ArchRProj = proj,
    name = "vierstra_motif_archetypeMatrix",
    plot = TRUE,
    n = 10
)

plotPDF(plot_var_dev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

# save the variability output
var_dev <- getVarDeviations(
    ArchRProj = proj,
    name = "vierstra_motif_archetypeMatrix",
    plot = FALSE
)

var_dev_df <- var_dev %>% as_tibble()
write_tsv(var_dev_df, "./vierstra_archetype_motif_chromVAR_variability_ranks.txt")

motif_id <- c("AC0259", "AC0227", "AC0622", "AC0610", "AC0606", "AC0261", "AC0489", "AC0405", "AC0024", "AC0525")
marker_motifs <- getFeatures(proj, select = paste(motif_id, collapse = "|"), useMatrix = "vierstra_motif_archetypeMatrix")
marker_motifs

marker_motifs_z <- grep("z:", marker_motifs, value = TRUE)
marker_motifs_z

p <- plotGroups(
    ArchRProj = proj,
    groupBy = "group",
    colorBy = "vierstra_motif_archetypeMatrix",
    name = marker_motifs_z,
    imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(seq_along(p), function(x) {
    if (x != 1) {
        p[[x]] + guides(color = FALSE, fill = FALSE) +
            theme_ArchR(baseSize = 8) +
            theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank()
            ) + ylab("")
    } else {
        p[[x]] + guides(color = FALSE, fill = FALSE) +
            theme_ArchR(baseSize = 8) +
            theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
            theme(
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank()
            ) + ylab("")
    }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))), p2))

plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 12, height = 12, ArchRProj = proj, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "vierstra_motif_archetypeMatrix",
    name = sort(marker_motifs_z),
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(p, function(x) {
    x + guides(color = FALSE, fill = FALSE) +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
})
do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation-UMAP", width = 12, height = 12, ArchRProj = proj, addDOC = FALSE)

# save ArchR object
saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep", load = FALSE)

# session info
sessionInfo()