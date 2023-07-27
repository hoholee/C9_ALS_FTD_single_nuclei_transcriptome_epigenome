# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(DESeq2)
library(variancePartition)

# read in pseudobulk grouping metadata
sample_info <- read_tsv("./pseudobulk_group_metadata.tsv")

# read the pseudobulk counts genereated from 13_2.get_pseudobulk_counts_by_individual.R
mat_summary_mm <- readRDS("snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_RNA_raw_count.rds")
mat_summary_mm_SCT <- readRDS("snRNA_pseudobulk_by_region_disease_individual_annoLevel2_dgCMatrix_SCT_raw_count.rds")

run_PCA_by_group <- function(selected_cell_type, selected_region, selected_diseases) {
  selected_samples <- sample_info %>%
    filter(
      cell_type == selected_cell_type,
      region == selected_region,
      disease %in% selected_diseases
    ) %>%
    pull(group)

  mat_sub <- mat_summary_mm[, selected_samples]

  sample_info_df <- sample_info %>%
    filter(group %in% selected_samples) %>%
    mutate(group = factor(group, levels = selected_samples)) %>%
    arrange(group) %>%
    as.data.frame() %>%
    column_to_rownames("group")

  dds <- DESeqDataSetFromMatrix(
    countData = mat_sub,
    colData = sample_info_df,
    design = ~ disease + sex
  )

  # Perform median ratio normalization
  dds_norm <- estimateSizeFactors(dds)
  norm_counts <- counts(dds_norm, normalized = TRUE)

  # Identify constant or zero-variance genes
  zero_variance_genes <- which(apply(norm_counts, 1, var) == 0)

  # Remove constant or zero-variance genes from the data
  norm_counts_filtered <- norm_counts[-zero_variance_genes, ]
  norm_counts_filtered_log <- log1p(norm_counts_filtered)

  # Perform PCA
  pca_res <- prcomp(t(norm_counts_filtered_log), scale = TRUE, center = TRUE)

  # Convert PCA results and metadata into a data frame
  data.frame(pca_res$x, colData(dds_norm)) %>% as_tibble()
}

cell_type_list <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

res <- map_dfr(
  cell_type_list,
  run_PCA_by_group,
  selected_region = "MCX",
  selected_diseases = c("ALS", "Control")
) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_list))

# Create PCA plot with ellipses
res %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = disease, shape = sex), size = 1) +
  # stat_ellipse(level = 0.95, type = "t", linetype = "dashed") +
  facet_wrap(~cell_type) +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )

# Display PCA plot with ellipses
pca_plot

# Transform counts for data visualization
rld <- rlog(dds_norm, blind = TRUE)

# Plot PCA
DESeq2::plotPCA(rld, intgroup = "disease", ntop = 3000, returnData = TRUE)


# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds_norm) > 1) >= 0.5 * ncol(dds_norm)

quantLog <- log2(fpm(dds_norm)[isexpr, ] + 1)

# define formula
# form <- ~disease + sex + age + (1|subject)
form <- ~ disease + sex + age
form <- ~ disease + sex + age + cell_type

# run variance partitioning
varPart <- fitExtractVarPartModel(quantLog, form, sample_info_df)
# varPart <- fitExtractVarPartModel(quantLog, form, sample_info)

# violin plot
plotVarPart(sortCols(varPart), label.angle = 60)

# sample_info_full <- sample_info %>%
#   mutate(group = factor(group, levels = colnames(mat_summary_mm))) %>%
#   arrange(group) %>%
#   as.data.frame() %>%
#   column_to_rownames("group")

# Create a DESeqDataSet object
# dds <- DESeqDataSetFromMatrix(
#   countData = mat_summary_mm + 1,
#   colData = sample_info_full,
#   design = ~ disease + sex
# )
## log sessionInfo
sessionInfo()
