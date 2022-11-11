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

# read meta data with annotations
# meta_data_anno <- read_tsv("./metadata_merged_addAnno.txt")

marker_all <- readRDS("snATAC_marker_peaks_by_annoLevel2_region_disease.rds")

heatmapPeaks <- markerHeatmap(
    seMarker = marker_all,
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    transpose = TRUE,
    nLabel = 1
)

# draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

heatmapPeaks <- markerHeatmap(
    seMarker = marker_all,
    cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
    transpose = TRUE,
    nLabel = 1
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap2", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

# add motif annotaions
# Use Jeff Vierstra's Non-redundant clustering Archetypes (v2beta)
proj <- addMotifAnnotations(
    ArchRProj = proj,
    motifSet = "vierstra",
    collection = "archetype",
    annoName = "vierstra_motif_archetype",
    cutOff = 5e-05,
    width = 7
)

# test motif enrichment in all marker peaks
enrich_motifs <- peakAnnoEnrichment(
    seMarker = marker_all,
    ArchRProj = proj,
    peakAnnotation = "vierstra_motif_archetype",
    cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
    background = "all" # consider changing this to "bgdPeaks"
)

saveRDS(enrich_motifs, "enriched_Vierstra_motifs_in_snATAC_marker_peaks_by_annoLevel2_region_disease.rds")

# plot heatmap of enriched motifs
heatmap_EM <- plotEnrichHeatmap(enrich_motifs, n = 10, transpose = TRUE)
plotPDF(heatmap_EM, name = "Motifs-Enriched-Marker-Heatmap", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

heatmap_EM <- plotEnrichHeatmap(enrich_motifs, n = 10, transpose = TRUE, clusterCols = FALSE)
plotPDF(heatmap_EM, name = "Motifs-Enriched-Marker-Heatmap3", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

# ArchR ATAC collection
proj <- addArchRAnnotations(
    ArchRProj = proj,
    db = "ArchR",
    collection = "ATAC"
)

enrich_ATAC <- peakAnnoEnrichment(
    seMarker = marker_all,
    ArchRProj = proj,
    peakAnnotation = "ATAC",
    cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
    background = "all" # consider changing this to "bgdPeaks"
)

saveRDS(enrich_motifs, "enriched_ArchR_ATACcollection_in_snATAC_marker_peaks_by_annoLevel2_region_disease.rds")

heatmap_ATAC <- plotEnrichHeatmap(enrich_ATAC, n = 10, transpose = TRUE, clusterCols = FALSE)
plotPDF(heatmap_ATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

# ArchR Encode TFBS collection
proj <- addArchRAnnotations(
    ArchRProj = proj,
    db = "ArchR",
    collection = "EncodeTFBS"
)

enrich_Encode <- peakAnnoEnrichment(
    seMarker = marker_all,
    ArchRProj = proj,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
    background = "all" # consider changing this to "bgdPeaks"
)

saveRDS(enrich_Encode, "enriched_ArchR_EncodeTFBScollection_in_snATAC_marker_peaks_by_annoLevel2_region_disease.rds")

heatmap_encodeTFBS <- plotEnrichHeatmap(enrich_Encode, n = 10, transpose = TRUE, clusterCols = FALSE)
plotPDF(heatmap_encodeTFBS, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

# ArchR Cistrome TFBS collection
proj <- addArchRAnnotations(
    ArchRProj = proj,
    db = "ArchR",
    collection = "CistromeTFBS"
)

enrich_Cistrome <- peakAnnoEnrichment(
    seMarker = marker_all,
    ArchRProj = proj,
    peakAnnotation = "CistromeTFBS",
    cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
    background = "all" # consider changing this to "bgdPeaks"
)

saveRDS(enrich_Cistrome, "enriched_ArchR_CistromeTFBScollection_in_snATAC_marker_peaks_by_annoLevel2_region_disease.rds")

heatmap_CistromeTFBS <- plotEnrichHeatmap(enrich_Cistrome, n = 10, transpose = TRUE, clusterCols = FALSE)
plotPDF(heatmap_CistromeTFBS, name = "CistromeTFBS-Enriched-Marker-Heatmap", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

# ArchR Codex collection
proj <- addArchRAnnotations(
    ArchRProj = proj,
    db = "ArchR",
    collection = "Codex"
)

enrich_Codex <- peakAnnoEnrichment(
    seMarker = marker_all,
    ArchRProj = proj,
    peakAnnotation = "Codex",
    cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
    background = "all" # consider changing this to "bgdPeaks"
)

saveRDS(enrich_Codex, "enriched_ArchR_CODEXcollection_in_snATAC_marker_peaks_by_annoLevel2_region_disease.rds")

heatmap_Codex <- plotEnrichHeatmap(enrich_Codex, n = 10, transpose = TRUE, clusterCols = FALSE)
plotPDF(heatmap_Codex, name = "CODEX-Enriched-Marker-Heatmap", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

# Custom enrichment test with our H3K27ac ChIP-seq peaks
H3K27ac_peaks <- c(
    H3K27ac_NeuN = "./H3K27ac_ChIPseq_peaks/H3K27ac_ChIPseq_peaks_ALS_vs_Control_in_NeuN.bed",
    H3K27ac_Astro = "./H3K27ac_ChIPseq_peaks/H3K27ac_ChIPseq_peaks_ALS_vs_Control_in_Astro.bed",
    H3K27ac_MG = "./H3K27ac_ChIPseq_peaks/H3K27ac_ChIPseq_peaks_ALS_vs_Control_in_MG.bed",
    H3K27ac_Oligo = "./H3K27ac_ChIPseq_peaks/H3K27ac_ChIPseq_peaks_ALS_vs_Control_in_Oligo.bed"
)

proj <- addPeakAnnotations(
    ArchRProj = proj,
    regions = H3K27ac_peaks,
    name = "H3K27ac_ChIPseq"
)

enrich_H3K27ac <- peakAnnoEnrichment(
    seMarker = marker_all,
    ArchRProj = proj,
    peakAnnotation = "H3K27ac_ChIPseq",
    cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
    background = "all" # consider changing this to "bgdPeaks"
)

saveRDS(enrich_H3K27ac, "enriched_H3K27ac_ChIPseqPeaks_in_snATAC_marker_peaks_by_annoLevel2_region_disease.rds")

heatmap_H3K27ac <- plotEnrichHeatmap(enrich_H3K27ac, n = 10, transpose = TRUE, clusterCols = FALSE)
plotPDF(heatmap_H3K27ac, name = "H3K27ac-Enriched-Marker-Heatmap", width = 20, height = 15, ArchRProj = proj, addDOC = FALSE)

# save ArchR object
saveArchRProject(ArchRProj = proj, outputDirectory = "./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep", load = FALSE)
# session info
sessionInfo()