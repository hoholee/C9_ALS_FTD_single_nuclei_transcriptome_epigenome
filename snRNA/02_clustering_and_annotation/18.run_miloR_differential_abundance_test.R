# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(miloR)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(tidyverse)
library(patchwork)

data_obj_sub <- readRDS("./seurat_objects/snRNA_postQC_scTransformed_clustered_2ndRoundAnnotated_clean_SeuratV4_object.rds")

data_obj_sub <- FindNeighbors(data_obj_sub,
                              reduction = "harmony",
                              dims = 1:30,
                              k.param = 30)
data_sce <- as.SingleCellExperiment(data_obj_sub)
data_milo <- Milo(data_sce)
# miloR::graph(data_milo) <- miloR::graph(buildFromAdjacency(as(data_obj_sub@graphs$SCT_snn, "dgCMatrix"), k = 30))

data_milo <- buildGraph(data_milo, k = 30, d = 30, reduced.dim = "HARMONY")

data_milo <- makeNhoods(data_milo, prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = "HARMONY")

plotNhoodSizeHist(data_milo)

data_milo <- countCells(data_milo, meta.data = data.frame(colData(data_milo)), sample = "orig.ident")

head(nhoodCounts(data_milo))

data_design <- data.frame(colData(data_milo))[, c("orig.ident", "disease", "region", "sex", "seq_batch")]
data_design <- distinct(data_design)
data_design$group <- paste(data_design$disease, data_design$region, sep = "_")
rownames(data_design) <- data_design$orig.ident
data_design

# this step took ~3h to finish
data_milo <- calcNhoodDistance(data_milo, d = 30)

# da_results <- testNhoods(data_milo, design = ~ 0+ disease, design.df = data_design,
#                          model.contrasts = "diseaseFTD - diseaseControl")
# da_results <- testNhoods(data_milo, design = ~ 0 + disease + region + sex + seq_batch, design.df = data_design,
#                          model.contrasts = "diseaseALS - diseaseControl")
da_results <- testNhoods(data_milo, design = ~ 0 + group, design.df = data_design,
                         model.contrasts = "groupALS_MCX - groupControl_MCX")
da_results <- testNhoods(data_milo, design = ~ 0 + group, design.df = data_design,
                         model.contrasts = "groupALS_mFCX - groupControl_mFCX")
da_results <- testNhoods(data_milo, design = ~ 0 + group, design.df = data_design,
                         model.contrasts = "groupFTD_MCX - groupControl_MCX")
da_results <- testNhoods(data_milo, design = ~ 0 + group, design.df = data_design,
                         model.contrasts = "groupFTD_mFCX - groupControl_mFCX")
# add other covariants
da_results <- testNhoods(data_milo, design = ~ 0 + group +sex + seq_batch, design.df = data_design,
                         model.contrasts = "groupFTD_MCX - groupControl_MCX")

da_results %>%
  arrange(SpatialFDR) %>%
  head()

data_milo <- buildNhoodGraph(data_milo)

plotUMAP(data_milo, colour_by = "rna_anno_2ndRound_level_2") + plotNhoodGraphDA(data_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

da_results <- annotateNhoods(data_milo, da_results, coldata_col = "rna_anno_2ndRound_level_2")
head(da_results)

ggplot(da_results, aes(rna_anno_2ndRound_level_2_fraction)) + geom_histogram(bins=50)

da_results$celltype <- ifelse(da_results$rna_anno_2ndRound_level_2_fraction < 0.7, "Mixed", da_results$rna_anno_2ndRound_level_2)

plotDAbeeswarm(da_results, group.by = "celltype")

da_results <- annotateNhoods(data_milo, da_results, coldata_col = "rna_anno_2ndRound_level_3")
head(da_results)

ggplot(da_results, aes(rna_anno_2ndRound_level_3_fraction)) + geom_histogram(bins=50)

da_results$celltype <- ifelse(da_results$rna_anno_2ndRound_level_3_fraction < 0.7, "Mixed", da_results$rna_anno_2ndRound_level_3)

plotDAbeeswarm(da_results, group.by = "celltype")

## log sessionInfo
sessionInfo()
