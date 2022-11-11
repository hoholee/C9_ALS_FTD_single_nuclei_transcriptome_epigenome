# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
set.seed(666)

addArchRThreads(threads = 8)

addArchRGenome("hg38")

# full dataset, snRNA label transferred with level 2
proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_cellBender_level2_cleanV2_addPseudobulkRep/")

metadata_full_level_2 <- getCellColData(proj) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()

umap_full <- getEmbedding(proj, embedding = "UMAP", returnDF = TRUE) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    dplyr::rename(
        UMAP_1 = "IterativeLSI#UMAP_Dimension_1",
        UMAP_2 = "IterativeLSI#UMAP_Dimension_2"
    )

write_tsv(metadata_full_level_2, "metadata_full_dataset_add_snRNA_labelTransferred_level2.txt")
write_tsv(umap_full, "umap_full_dataset.txt")

# full dataset, snRNA label transferred with level 3
proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_cellBender_level3_cleanV2")

metadata_full_level_3 <- getCellColData(proj) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()

write_tsv(metadata_full_level_3, "metadata_full_dataset_add_snRNA_labelTransferred_level3.txt")

# sub-clustering
# Exc-Neurons
proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_Neu_Exc_subset")
metadata_Exc <- getCellColData(proj) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()

umap_Exc <- getEmbedding(proj, embedding = "UMAP", returnDF = TRUE) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    dplyr::rename(
        UMAP_1 = "IterativeLSI#UMAP_Dimension_1",
        UMAP_2 = "IterativeLSI#UMAP_Dimension_2"
    )

write_tsv(metadata_Exc, "metadata_Exc_add_snRNA_labelTransferred.txt")
write_tsv(umap_Exc, "umap_Exc.txt")

# Inh-Neurons
proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_Neu_Inh_subset")
metadata_Inh <- getCellColData(proj) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()

umap_Inh <- getEmbedding(proj, embedding = "UMAP", returnDF = TRUE) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    dplyr::rename(
        UMAP_1 = "IterativeLSI#UMAP_Dimension_1",
        UMAP_2 = "IterativeLSI#UMAP_Dimension_2"
    )

write_tsv(metadata_Inh, "metadata_Inh_add_snRNA_labelTransferred.txt")
write_tsv(umap_Inh, "umap_Inh.txt")

# Glia
proj <- loadArchRProject(path = "Dracheva_ALS_FTD_addRNA_Glia_subset")

metadata_Glia <- getCellColData(proj) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()

umap_Glia <- getEmbedding(proj, embedding = "UMAP", returnDF = TRUE) %>%
    rownames_to_column("cell_id") %>%
    as_tibble() %>%
    dplyr::rename(
        UMAP_1 = "IterativeLSI#UMAP_Dimension_1",
        UMAP_2 = "IterativeLSI#UMAP_Dimension_2"
    )

write_tsv(metadata_Glia, "metadata_Glia_add_snRNA_labelTransferred.txt")
write_tsv(umap_Glia, "umap_Glia.txt")

# merge all clusters, snRNA transferred labels, prediciton scores and UMAP coordinates from full dataset and sub-clustering
# remove `ReadsInPeaks` and `FRIP` as these depend on peak calling, which in turn depends on how the cells are grouped

metadata_full_level_2_mod <- metadata_full_level_2 %>%
    select(-ReadsInPeaks, -FRIP) %>%
    dplyr::rename(
        predictedCell_from_snRNA_full_level2 = predictedCell_Un,
        predictedGroup_from_snRNA_full_level2 = predictedGroup_Un,
        predictedScore_from_snRNA_full_level2 = predictedScore_Un,
        Clusters_full = Clusters,
        Clusters_major_from_full_level2 = Cluster_major,
        cell_group_by_majorClusters_region_disease = Cluster_sep
    )

metadata_full_level_3_mod <- metadata_full_level_3 %>%
    select(cell_id, predictedCell_Un, predictedGroup_Un, predictedScore_Un) %>%
    dplyr::rename(
        predictedCell_from_snRNA_full_level3 = predictedCell_Un,
        predictedGroup_from_snRNA_full_level3 = predictedGroup_Un,
        predictedScore_from_snRNA_full_level3 = predictedScore_Un,
    )

metadata_Exc_mod <- metadata_Exc %>%
    select(
        cell_id, Clusters,
        predictedCell_Un_level2, predictedGroup_Un_level2, predictedScore_Un_level2,
        predictedCell_Un, predictedGroup_Un, predictedScore_Un,
    ) %>%
    dplyr::rename(
        Clusters_sub_Exc = Clusters,
        predictedCell_from_snRNA_Exc_level2 = predictedCell_Un_level2,
        predictedGroup_from_snRNA_Exc_level2 = predictedGroup_Un_level2,
        predictedScore_from_snRNA_Exc_level2 = predictedScore_Un_level2,
        predictedCell_from_snRNA_Exc_level3 = predictedCell_Un,
        predictedGroup_from_snRNA_Exc_level3 = predictedGroup_Un,
        predictedScore_from_snRNA_Exc_level3 = predictedScore_Un,
    )

metadata_Inh_mod <- metadata_Inh %>%
    select(
        cell_id, Clusters,
        predictedCell_Un_level2, predictedGroup_Un_level2, predictedScore_Un_level2,
        predictedCell_Un, predictedGroup_Un, predictedScore_Un,
    ) %>%
    dplyr::rename(
        Clusters_sub_Inh = Clusters,
        predictedCell_from_snRNA_Inh_level2 = predictedCell_Un_level2,
        predictedGroup_from_snRNA_Inh_level2 = predictedGroup_Un_level2,
        predictedScore_from_snRNA_Inh_level2 = predictedScore_Un_level2,
        predictedCell_from_snRNA_Inh_level3 = predictedCell_Un,
        predictedGroup_from_snRNA_Inh_level3 = predictedGroup_Un,
        predictedScore_from_snRNA_Inh_level3 = predictedScore_Un,
    )

metadata_Glia_mod <- metadata_Glia %>%
    select(
        cell_id, Clusters,
        predictedCell_Un_level2, predictedGroup_Un_level2, predictedScore_Un_level2,
        predictedCell_Un, predictedGroup_Un, predictedScore_Un,
    ) %>%
    dplyr::rename(
        Clusters_sub_Glia = Clusters,
        predictedCell_from_snRNA_Glia_level2 = predictedCell_Un_level2,
        predictedGroup_from_snRNA_Glia_level2 = predictedGroup_Un_level2,
        predictedScore_from_snRNA_Glia_level2 = predictedScore_Un_level2,
        predictedCell_from_snRNA_Glia_level3 = predictedCell_Un,
        predictedGroup_from_snRNA_Glia_level3 = predictedGroup_Un,
        predictedScore_from_snRNA_Glia_level3 = predictedScore_Un,
    )

umap_full_mod <- umap_full %>%
    dplyr::rename(
        full_UMAP_1 = UMAP_1,
        full_UMAP_2 = UMAP_2
    )

umap_Exc_mod <- umap_Exc %>%
    dplyr::rename(
        Exc_UMAP_1 = UMAP_1,
        Exc_UMAP_2 = UMAP_2
    )

umap_Inh_mod <- umap_Inh %>%
    dplyr::rename(
        Inh_UMAP_1 = UMAP_1,
        Inh_UMAP_2 = UMAP_2
    )

umap_Glia_mod <- umap_Glia %>%
    dplyr::rename(
        Glia_UMAP_1 = UMAP_1,
        Glia_UMAP_2 = UMAP_2
    )

merged_metadata <- metadata_full_level_2_mod %>%
    left_join(metadata_full_level_3_mod, by = "cell_id") %>%
    left_join(metadata_Exc_mod, by = "cell_id") %>%
    left_join(metadata_Inh_mod, by = "cell_id") %>%
    left_join(metadata_Glia_mod, by = "cell_id") %>%
    left_join(umap_full_mod, by = "cell_id") %>%
    left_join(umap_Exc_mod, by = "cell_id") %>%
    left_join(umap_Inh_mod, by = "cell_id") %>%
    left_join(umap_Glia_mod, by = "cell_id")

write_tsv(merged_metadata, "metadata_merged.txt")

# log session info
sessionInfo()