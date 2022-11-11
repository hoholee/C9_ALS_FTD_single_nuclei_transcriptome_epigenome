# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)
library(viridis)
library(circlize)
set.seed(666)

addArchRThreads(threads = 16)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep")

# read meta data with annotations
# meta_data_anno <- read_tsv("./metadata_merged_addAnno.txt")

marker_all <- readRDS("snATAC_marker_peaks_by_annoLevel2_region_disease.rds")

# subset samples from the marker object
sample_index <- grep("__", colnames(marker_all), invert = TRUE)
marker_all_sub <- marker_all[, sample_index]

heatmapPeaks <- plotMarkerHeatmap(
    seMarker = marker_all_sub,
    cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
    transpose = TRUE,
    nLabel = 1,
    # returnMatrix = TRUE
    returnMatrix = FALSE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap4", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

marker_all <- readRDS("snATAC_marker_peaks_by_annoLevel2_clean.rds")

heatmapPeaks <- plotMarkerHeatmap(
    seMarker = marker_all,
    cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
    transpose = TRUE,
    nLabel = 1,
    returnMatrix = FALSE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap5", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)


mat_heatmapPeaks <- plotMarkerHeatmap(
    seMarker = marker_all,
    cutOff = "FDR < 0.05 & Log2FC > 0.5",
    transpose = TRUE,
    nLabel = 1,
    returnMatrix = TRUE
)

color_fun <- colorRamp2(seq(-2, 2, length.out = 10), viridis(n = 10, option = "D"))
mat_row_order <- c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep",
    "Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST",
    "Astro", "Micro", "Oligo", "OPC"
)

pdf(("test_marker_peaks_heatmap_clean.pdf"), width = 8, height = 8)
ht <- Heatmap(
    matrix = mat_heatmapPeaks,
    col = color_fun,
    name = "Column Z-scores",
    na_col = "grey",
    # rect_gp = gpar(col = "white", lwd = 2),
    # border_gp = gpar(col = "black", lty = 2),
    color_space = "LAB",
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    row_order = mat_row_order,
    #   clustering_distance_rows = "euclidean",
    #   clustering_method_rows = "ward.D2",
    #   row_dend_reorder = TRUE,
    #   show_row_dend = TRUE,
    cluster_columns = FALSE,
    # clustering_distance_columns = "euclidean",
    # clustering_method_columns = "ward.D2",
    # column_dend_reorder = TRUE,
    # show_column_dend = FALSE,
    #   column_order = mat_col_order_2,
    #   row_km = 7,
    #   row_km_repeats = 100,
    # row_split = 8,
    # cluster_row_slices = FALSE,
    #   top_annotation = mat_col_anno_2,
    #   right_annotation = mat_row_anno_2,
    use_raster = TRUE,
    raster_by_magick = TRUE,
    raster_magick_filter = "Lanczos",
    raster_device = "CairoPNG",
    width = unit(3, "inches"),
    height = unit(2, "inches"),
    #   heatmap_width = unit(4, "inches"),
    #   heatmap_height = unit(4, "inches"),
)
draw(ht)
dev.off()

# save ArchR object
# saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep", load = FALSE)
# session info
sessionInfo()