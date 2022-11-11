# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(Seurat)

## set data path and sample info
data_dir <- "/cndd2/junhao/ALS_FTD_singleCell/run_cellBender_on_raw_snRNA/"
sample_meta <- read_tsv(paste0(data_dir, "sample_metadata.txt"), col_types = "cccccddc")

# exclude samples from subject 908, which turns out to be not FTD patient
# also exclude sample mFCX_FTD_54, which failed in cellBender pre-processing (and also other QC)
sample_meta_sub <- sample_meta %>% filter(subject != "908", sample_id != "mFCX_FTD_54")

## make a function to load the corrected count data from cellBender
read_count_data_for_multiple_samples <- function(selected_sample){
   message(paste0("Reading data from: ", selected_sample, "..."))
   # load the count data from cellBender generated h5 files
   data <- Read10X_h5(filename = paste0(data_dir, selected_sample, "_cellBender_corrected_filtered.h5"),
                      use.names = TRUE, unique.features = TRUE)
   # initialize the Seurat object
   obj <- CreateSeuratObject(counts = data,
                             project = selected_sample,
                             assay = "RNA",
                             min.cells = 3,
                             min.features = 200)
   # add metadata to the Seurat object
   meta_selected <- sample_meta_sub %>% filter(sample_id == selected_sample)
   obj$region <- meta_selected$region
   obj$disease <- meta_selected$disease
   obj$subject <- meta_selected$subject
   obj$sex <- meta_selected$sex
   obj$age <- meta_selected$age
   obj$PMI <- meta_selected$PMI
   obj$seq_batch <- meta_selected$seq_batch

   obj
}


## read data
obj_list <- map(sample_meta_sub$sample_id, read_count_data_for_multiple_samples)

## merge into a big Seurat object
data_obj <- merge(obj_list[[1]],
                  y = obj_list[2:length(obj_list)],
                  add.cell.ids = sample_meta_sub$sample_id,
                  project = "Dracheva_ALS_FTD_cellBender_corrected")

## add the percentage of reads that map to the mitochondrial genome as one of the QC metrics
## edit the regrex for MT genes if needed
data_obj[["percent_mt"]] <- PercentageFeatureSet(data_obj, pattern = "^MT-")
saveRDS(data_obj, file = "./seurat_objects/snRNA_cellBender_corrected_SeuratV4_object.rds")

## log sessionInfo
sessionInfo()
