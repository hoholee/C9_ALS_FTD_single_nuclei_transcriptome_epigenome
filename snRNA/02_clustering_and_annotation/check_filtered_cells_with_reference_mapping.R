# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)
library(harmony)
library(ggrastr)

## read seurat object
data_obj_map_to_manuscript <- readRDS("./seurat_objects/snRNA_cellBender_corrected_QC_filteredCells_addprediction_MappedToManuscriptClusters_SeuratV4_object.rds")
data_obj_map_to_AIBS <- readRDS("./seurat_objects/snRNA_cellBender_corrected_QC_filteredCells_addprediction_MappedToAIBS_SeuratV4_object.rds")

meta_data_map_to_manuscript <- data_obj_map_to_manuscript@meta.data %>%
       rownames_to_column("cell_id") %>%
       as_tibble()

meta_data_map_to_AIBS <- data_obj_map_to_AIBS@meta.data %>%
       rownames_to_column("cell_id") %>%
       as_tibble()

p <- meta_data_map_to_manuscript %>% ggplot(aes(prediction.score.max))
p +
       geom_density() +
       facet_wrap(~predicted.id, scales = "free") +
       theme_bw(base_size = 10, base_family = "Helvetica") +
       theme(
              strip.background = element_blank()
       )

p <- meta_data_map_to_manuscript %>% ggplot(aes(prediction.score.max))
p +
       stat_ecdf() +
       facet_wrap(~predicted.id, scales = "free") +
       theme_bw(base_size = 10, base_family = "Helvetica") +
       theme(
              strip.background = element_blank()
       )

p <- meta_data_map_to_AIBS %>% ggplot(aes(prediction.score.max))
p +
       geom_density() +
       facet_wrap(~predicted.id, scales = "free") +
       theme_bw(base_size = 10, base_family = "Helvetica") +
       theme(
              strip.background = element_blank()
       )

## log sessionInfo
sessionInfo()
