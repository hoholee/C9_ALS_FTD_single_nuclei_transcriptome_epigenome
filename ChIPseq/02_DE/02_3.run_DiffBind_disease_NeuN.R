# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(DiffBind)

bam_path <- "/cndd3/Public_Datasets/Dracheva_ALS_ChIPseq/bam/"
peak_path <- "/cndd2/junhao/ALS_FTD_singleCell/ChIPseq/peaks_bed/"

subject_id <- c("113", "114", "115", "128", "129", "130",
                "302", "303", "309", "311", "312", "320", "322", "325", "328")

subject_sex <- c("F", "F", "F", "M", "M", "M",
                 "F", "F", "F", "F", "F", "M", "M", "M", "M")
names(subject_sex) <- subject_id

# disease effect: ALS vs. Control in NeuN
## read in peakSets with sampleSheet

sampleSheet <- tibble(
  SampleID = paste0(subject_id, "-NeuN"),
  Tissue = "MCX",
  Factor = subject_sex[subject_id], # sex as Factor
  Condition = if_else(grepl("^1", SampleID), "ALS", "Control"), # disease as Condition
  Replicate = c(1:6, 1:9),
  bamReads = paste0(bam_path, SampleID, ".merged.bam"),
  Peaks = paste0(peak_path, SampleID, "_peaks_broadPeak_keepChr.bed"),
  PeakCaller = "bed",
  PeakFormat = "bed",
  ScoreCol = 5
)

dba <- dba(sampleSheet = sampleSheet, bLowerScoreBetter = FALSE)
dba$config$singleEnd <- FALSE

## plot correlation matrix
pdf(file = "DiffBind_disease_ALS_vs_Control_in_NeuN_raw_cor_matrix.pdf")
plot(dba)
dev.off()

## apply blacklist
peakdata <- dba.show(dba)$Intervals
dba <- dba.blacklist(dba, blacklist = DBA_BLACKLIST_HG38, greylist = FALSE)
peakdata.BL <- dba.show(dba)$Intervals
peakdata - peakdata.BL

## check number of consensus peaks as a function of number of samples detected
rate <- dba.overlap(dba, mode = DBA_OLAP_RATE)

pdf(file = "DiffBind_disease_ALS_vs_Control_in_NeuN_consensus_peaks_across_samples.pdf")
plot(rate, type = 'b', xlab = "# peaksets", ylab = "# common peaks")
dev.off()

## Generate consensus for each condition
consensus_per_cond <- dba.peakset(dba, consensus = DBA_CONDITION, minOverlap = 2)

## Examine sample groups consensus peaks and overlap
dba.show(consensus_per_cond, consensus_per_cond$masks$Consensus)

pdf(file = "DiffBind_disease_ALS_vs_Control_in_NeuN_consensus_peaks_by_group.pdf")
dba.plotVenn(consensus_per_cond, consensus_per_cond$masks$Consensus)
dev.off()

## Retrieve union of consensus peaks
consensus <- dba.peakset(consensus_per_cond,
                         bRetrieve = TRUE,
                         peaks = consensus_per_cond$masks$Consensus,
                         minOverlap = 1)

## count reads
dba <- dba.count(dba,
                 peaks = consensus,
                 score = DBA_SCORE_NORMALIZED,
                 summits = TRUE,
                 bRemoveDuplicates = TRUE,
                 bUseSummarizeOverlaps = TRUE,
                 bParallel = TRUE)

pdf(file = "DiffBind_disease_ALS_vs_Control_in_NeuN_consensus_peaks_cor_matrix.pdf")
plot(dba)
dev.off()

## normalization
dba <- dba.normalize(dba,
                     method = DBA_ALL_METHODS,
                     normalize = DBA_NORM_NATIVE,
                     library = DBA_LIBSIZE_BACKGROUND,
                     background = TRUE)

## set up contrast
## add sex as a covariance
dba <- dba.contrast(dba,
                    design = "~Factor + Condition",
                    minMembers = 2)

## run the DE test
dba <- dba.analyze(dba,
                   method = DBA_ALL_METHODS,
                   bParallel = TRUE)
dba

save(dba, file = "diffBind_H3K27ac_ALS_vs_Control_in_NeuN_addCovSex_dba.rda")

pdf(file = "DiffBind_disease_ALS_vs_Control_in_NeuN_addCovSex_results.pdf", onefile = TRUE)
dba.plotHeatmap(dba, contrast = 1, method = DBA_EDGER)
dba.plotHeatmap(dba, contrast = 1, method = DBA_DESEQ2)
dba.plotVolcano(dba, contrast = 1, method = DBA_EDGER, bFlip = TRUE)
dba.plotVolcano(dba, contrast = 1, method = DBA_DESEQ2, bFlip = TRUE)
dba.plotMA(dba, contrast = 1, method = DBA_EDGER, bFlip = TRUE)
dba.plotMA(dba, contrast = 1, method = DBA_DESEQ2, bFlip = TRUE)
dba.plotBox(dba, contrast = 1, method = DBA_EDGER)
dba.plotBox(dba, contrast = 1, method = DBA_DESEQ2)
dba.plotHeatmap(dba, contrast = 1, method = DBA_EDGER, correlations = FALSE, scale = "row")
dba.plotHeatmap(dba, contrast = 1, method = DBA_DESEQ2, correlations = FALSE, scale = "row")
dba.plotVenn(dba, contrast = 1, bDB = TRUE, method = DBA_ALL_METHODS)

dba.plotHeatmap(dba, contrast = 2, method = DBA_EDGER)
dba.plotHeatmap(dba, contrast = 2, method = DBA_DESEQ2)
dba.plotVolcano(dba, contrast = 2, method = DBA_EDGER, bFlip = TRUE)
dba.plotVolcano(dba, contrast = 2, method = DBA_DESEQ2, bFlip = TRUE)
dba.plotMA(dba, contrast = 2, method = DBA_EDGER, bFlip = TRUE)
dba.plotMA(dba, contrast = 2, method = DBA_DESEQ2, bFlip = TRUE)
dba.plotBox(dba, contrast = 2, method = DBA_EDGER)
dba.plotBox(dba, contrast = 2, method = DBA_DESEQ2)
dba.plotHeatmap(dba, contrast = 2, method = DBA_EDGER, correlations = FALSE, scale = "row")
dba.plotHeatmap(dba, contrast = 2, method = DBA_DESEQ2, correlations = FALSE, scale = "row")
dba.plotVenn(dba, contrast = 2, bDB = TRUE, method = DBA_ALL_METHODS)
dev.off()

dba_edgeR_df <- dba.report(dba,
                           contrast = 1,
                           bCounts = TRUE,
                           method = DBA_EDGER,
                           th = 1,
                           bFlip = TRUE,
                           bCalled = TRUE,
                           bCalledDetail = TRUE,
                           bNormalized = TRUE,
                          ) %>%
  as_tibble()

dba_deseq2_df <- dba.report(dba,
                           contrast = 1,
                           bCounts = TRUE,
                           method = DBA_DESEQ2,
                           th = 1,
                           bFlip = TRUE,
                           bCalled = TRUE,
                           bCalledDetail = TRUE,
                           bNormalized = TRUE,
                          ) %>%
  as_tibble()

dba_edgeR_df2 <- dba.report(dba,
                           contrast = 2,
                           bCounts = TRUE,
                           method = DBA_EDGER,
                           th = 1,
                           bFlip = TRUE,
                           bCalled = TRUE,
                           bCalledDetail = TRUE,
                           bNormalized = TRUE,
                          ) %>%
  as_tibble()

dba_deseq2_df2 <- dba.report(dba,
                           contrast = 2,
                           bCounts = TRUE,
                           method = DBA_DESEQ2,
                           th = 1,
                           bFlip = TRUE,
                           bCalled = TRUE,
                           bCalledDetail = TRUE,
                           bNormalized = TRUE,
                          ) %>%
  as_tibble()

write_tsv(dba_edgeR_df, "diffBind_H3K27ac_ALS_vs_Control_in_NeuN_addCovSex_sexTerm_res_edgeR.txt")
write_tsv(dba_deseq2_df, "diffBind_H3K27ac_ALS_vs_Control_in_NeuN_addCovSex_sexTerm_res_DEseq2.txt")
write_tsv(dba_edgeR_df2, "diffBind_H3K27ac_ALS_vs_Control_in_NeuN_addCovSex_diseaseTerm_res_edgeR.txt")
write_tsv(dba_deseq2_df2, "diffBind_H3K27ac_ALS_vs_Control_in_NeuN_addCovSex_diseaseTerm_res_DEseq2.txt")
write_tsv(dba.show(dba) %>% as_tibble(), "diffBind_H3K27ac_ALS_vs_Control_in_NeuN_newConsensus_FRiP.txt")

# dba_edgeR_gr <- dba.report(dba,
#                            contrast = 2,
#                            bCounts = TRUE,
#                            method = DBA_EDGER,
#                            th = 0.05,
#                            bFlip = TRUE,
#                            bCalled = TRUE,
#                            bCalledDetail = TRUE,
#                            bNormalized = TRUE,
#                            bDB = TRUE,
#                            bAll = FALSE,
#                            bGain = TRUE,
#                            bLoss = TRUE
#                           )
#
# dba_deseq2_gr <- dba.report(dba,
#                            contrast = 2,
#                            bCounts = TRUE,
#                            method = DBA_DESEQ2,
#                            th = 0.05,
#                            bFlip = TRUE,
#                            bCalled = TRUE,
#                            bCalledDetail = TRUE,
#                            bNormalized = TRUE,
#                            bDB = TRUE,
#                            bAll = FALSE,
#                            bGain = TRUE,
#                            bLoss = TRUE
#                           )
#
# profiles_edgeR <- dba.plotProfile(dba, sites = dba_edgeR_gr,
#                                   # merge = NULL,
#                                   # scores = "Fold",
#                                   scores = "width",
#                                   absScores = TRUE, normalize = TRUE)
# dba.plotProfile(profiles_edgeR)
#
# profiles_deseq2 <- dba.plotProfile(dba, sites = dba_deseq2_gr,
#                                    # merge = NULL,
#                                    scores = "Fold", absScores = TRUE, normalize = TRUE)
# dba.plotProfile(profiles_deseq2)



## log sessionInfo
sessionInfo()
