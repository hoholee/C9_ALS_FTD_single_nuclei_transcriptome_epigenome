# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)
library(sctransform)
library(scales)
library(harmony)
## read the sctransform-normalized seurat object
data_obj_sub <- readRDS("./seurat_objects/snRNA_cellBender_corrected_dataset_SeuratV4_object_subset_Exc_scTransform_renormalized_annotated.rds")

# use the Leiden algorithm
# switch back to single worker to avoid serialize error
# plan("multiprocess", workers = 1)
# took ~1.5h to finish one resolution

leiden_resolution_list <- rev(c(0.05, 0.1, 0.3, 0.5, 0.8, 1, 1.8, 2, 2.5, 3, 3.5))

get_multiple_resolution_clusters <- function(selected_r){
	message(paste0("Handling resolution: ", selected_r))
	data_obj_mod <- FindClusters(data_obj_sub,
				     resolution = selected_r,
				     algorithm = 4, # 4 = Leiden
				     method = "matrix", # modify from "matrix" to "igraph" for large dataset
				     random.seed = 666,
				     verbose = TRUE)
	res_meta <- data_obj_mod@meta.data %>%
	  rownames_to_column("cell_id") %>%
	  as_tibble() %>%
	  select(cell_id, seurat_clusters) %>%
	  rename(seurat_NonNeu_subclusters = seurat_clusters) %>%
	  mutate(leiden_res = selected_r)
	write_tsv(res_meta, paste0("snRNA_subclustering_Exc_cellBender_corrected_leiden_clusters_res_", selected_r, ".txt"))
}

walk(leiden_resolution_list, get_multiple_resolution_clusters)

## log sessionInfo
sessionInfo()
