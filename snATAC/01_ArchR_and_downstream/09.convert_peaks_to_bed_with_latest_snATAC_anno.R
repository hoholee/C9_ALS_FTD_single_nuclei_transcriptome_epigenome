# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 16)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep")

# union peaks (iterative merged)
union_peaks <- getPeakSet(proj) %>%
  as_tibble()
union_peaks_bed <- union_peaks %>%
  select(seqnames, start, end) %>%
  mutate(start = start - 1L)

write_tsv(union_peaks, "union_peaks_cellBender_latest_snATAC_anno.txt")
write_tsv(union_peaks_bed, "union_peaks_cellBender_latest_snATAC_anno.bed", col_names = FALSE)

# get reproducible peaks set for each group
meta_data <- getCellColData(proj) %>% as_tibble()

group_list <- meta_data %>%
  distinct(group) %>%
  filter(!grepl("_Mixed__", group)) %>%
  filter(!group %in% c("Inh_LAMP5_MCX_FTD")) %>%
  arrange(group)

convert_rds_to_bed <- function(selected_group) {
  message(selected_group)
  peak <- readRDS(paste0("./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep/PeakCalls/", selected_group, "-reproduciblePeaks.gr.rds")) %>%
    as_tibble()
  peak_bed <- peak %>%
    select(seqnames, start, end) %>%
    mutate(start = start - 1L)
  write_tsv(peak, paste0("./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep/PeakCalls/", selected_group, "_reproduciblePeaks.txt"))
  write_tsv(peak_bed, paste0("./Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep/PeakCalls/", selected_group, "_reproduciblePeaks.bed"),
    col_names = FALSE
  )
}

walk(group_list$group, convert_rds_to_bed)

# session info
sessionInfo()