# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
library(tictoc)
set.seed(666)

addArchRThreads(threads = 8)

addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_cellBender_level2_cleanV2_addPseudobulkRep/")

## add impute weights
proj <- addImputeWeights(proj)

# meta_data_ori <- getCellColData(proj) %>% as_tibble()
meta_data <- read_tsv("./metadata_merged_addAnno.txt")

# add grouping by disease, region and `atac_anno_level_2`
meta_data_mod <- meta_data %>%
  mutate(sample = Sample) %>%
  separate(sample, c("region", "disease", "subject"), sep = "_") %>%
  mutate(group = paste(atac_anno_level_2, region, disease, sep = "_"))

proj <- addCellColData(
  ArchRProj = proj,
  data = meta_data_mod$group,
  name = "group",
  cells = getCellNames(proj)
)

call_DE_gene_score <- function(cell_group_1, cell_group_2) {
  message(paste0("Handling: ", cell_group_1, " vs ", cell_group_2))

  cell_group_1_mod <- cell_group_1 %>% str_replace_all(" ", "_")
  cell_group_2_mod <- cell_group_2 %>% str_replace_all(" ", "_")

  DE_gene_score <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "group",
    useGroups = cell_group_1,
    bgdGroups = cell_group_2,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    # normBy = "ReadsInTSS",
    k = 100,
    bufferRatio = 0.8
  )

  DE_df <- getMarkers(DE_gene_score, cutOff = "FDR <= 1")[[1]] %>%
    as_tibble() %>%
    mutate(
      group_1 = cell_group_1,
      group_2 = cell_group_2
    )

  saveRDS(DE_gene_score, paste0("./differential_gene_scores_by_snATAC_anno/DE_geneScores_", cell_group_1_mod, "_vs_", cell_group_2_mod, "_SEobj.rds"))
  write_tsv(DE_df, paste0("./differential_gene_scores_by_snATAC_anno/DE_geneScores_", cell_group_1_mod, "_vs_", cell_group_2_mod, "_res.txt"))
}

# test run and sanity check
call_DE_gene_score("Astro_MCX_Control", "Oligo_MCX_Control")

todo_celltypes <- meta_data_mod %>%
  distinct(atac_anno_level_2) %>%
  filter(!grepl("Mixed", atac_anno_level_2)) %>%
  arrange(atac_anno_level_2)
todo_group_1 <- paste0(todo_celltypes$atac_anno_level_2, "_MCX_ALS")
todo_group_2 <- paste0(todo_celltypes$atac_anno_level_2, "_MCX_Control")

walk2(todo_group_1, todo_group_2, call_DE_gene_score)

todo_group_1 <- paste0(todo_celltypes$atac_anno_level_2, "_mFCX_ALS")
todo_group_2 <- paste0(todo_celltypes$atac_anno_level_2, "_mFCX_Control")

walk2(todo_group_1, todo_group_2, call_DE_gene_score)

# log session info
sessionInfo()