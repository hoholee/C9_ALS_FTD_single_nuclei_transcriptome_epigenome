# Load library
library(ggbiplot)
library(tidyverse)
library(edgeR)
library(magrittr)
library(viridis)
library(gridExtra)
library(ggrepel)

# Read in data and tidy up
# remove subject 908
exp_data <- read_tsv("allSamples_rsem_genes_results_counts_annotated.txt") %>% 
  select(-starts_with("908"))
tpm_data <- read_tsv("allSamples_rsem_genes_results_TPM_annotated.txt") %>% 
  select(-starts_with("908"))

sample_meta <- read_tsv("sample_metadata.txt", col_types = "ccccdcc") %>% 
  filter(subject_id != "908")

gene_annotation <- read_tsv("/cndd2/junhao/genome/hg38/star_rsem_idx/gencode.v35.gene.annotation.txt")

# Run PCA

# Define ggscreeplot

ggscreeplot <- function(pcobj, type = c('pev', 'cev')) 
{
  type <- match.arg(type)
  d <- pcobj$sdev^2
  yvar <- switch(type, 
                 pev = d / sum(d), 
                 cev = cumsum(d) / sum(d))
  
  yvar.lab <- switch(type,
                     pev = 'proportion of explained variance',
                     cev = 'cumulative proportion of explained variance')
  
  df <- data.frame(PC = factor(1:length(d)), yvar = yvar)
  
  ggplot(data = df, aes(x = PC, y = yvar, group=1)) + 
    xlab('principal components') + ylab(yvar.lab) +
    geom_point() + geom_path() + theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
}

# all samples
tpm_data_rm0 <- tpm_data %>% mutate(RowSum = rowSums(select(., -starts_with("gene")))) %>% filter(RowSum>0)
# remove mitochondria and ribosomal genes
# tpm_data_rm0 <- tpm_data %>% filter(!(grepl("^Mt_", geneType) | grepl("rRNA", geneType))) %>%
#   mutate(RowSum = rowSums(select(., -starts_with("gene")))) %>% filter(RowSum>0)
t_tpmData <- tpm_data_rm0 %>% select(-starts_with("gene"), -RowSum)
t_tpmData <- log10(t_tpmData + 1) %>% t()

colnames(t_tpmData) <- tpm_data_rm0$geneID
pca_tpm <- prcomp(t_tpmData, center = T, scale. = T)
ggscreeplot(pca_tpm)
ggsave(file = "./plots/pca/PCA_allSamples_log10TPM_screePlot.pdf", device = cairo_pdf(), width = 10, height = 4, useDingbats = F)

allPCcombn <- combn(3,2)
plist <- apply(allPCcombn,2,function(x){
  choices <- x
  p <- ggbiplot(
    pca_tpm,
    pc.biplot = T,
    var.axes = F,
    # labels = sample_meta$sample,
    groups = sample_meta$cell_type,
    # groups = sample_meta$disease,
    ellipse = T,
    choices = choices
  ) + theme_bw(base_size = 16) + 
    scale_color_manual(values = c("#377EB8","#4DAF4A","#E41A1C")) +
    # scale_color_brewer(palette = "Set1") +
    theme(legend.position = "top") +
    coord_cartesian()
  p
})
plist[["ncol"]] <- 3
do.call(grid.arrange, plist)

p <- ggbiplot(
  pca_tpm,
  pc.biplot = T,
  var.axes = F,
  groups = sample_meta$cell_type,
  ellipse = T,
  choices = c(1, 2)
) + theme_bw(base_size = 10) + 
  scale_color_manual(values = c("#377EB8","#4DAF4A","#E41A1C")) +
  theme(legend.position = "top", panel.grid.minor = element_blank()) +
  coord_cartesian()
p
ggsave("./plots/pca/PCA_allSamples_log10TPM_PC1vsPC2_colorByCelltype.pdf", p, device = cairo_pdf(),
       width = 6, height = 4, useDingbats = F)

p <- ggbiplot(
  pca_tpm,
  pc.biplot = T,
  var.axes = F,
  groups = sample_meta$disease,
  ellipse = T,
  choices = c(1, 2)
) + theme_bw(base_size = 10) + 
  scale_color_manual(values = c("#377EB8","#4DAF4A","#E41A1C")) +
  theme(legend.position = "top", panel.grid.minor = element_blank()) +
  coord_cartesian()
p
ggsave("./plots/pca/PCA_allSamples_log10TPM_PC1vsPC2_colorByDisease.pdf", p, device = cairo_pdf(),
       width = 6, height = 4, useDingbats = F)

p <- ggbiplot(
  pca_tpm,
  pc.biplot = T,
  var.axes = F,
  groups = sample_meta$region,
  ellipse = F,
  choices = c(1, 2)
) + theme_bw(base_size = 10) + 
  scale_color_manual(values = c("#377EB8","#4DAF4A","#E41A1C")) +
  theme(legend.position = "top", panel.grid.minor = element_blank()) +
  coord_cartesian()
p
ggsave("./plots/pca/PCA_allSamples_log10TPM_PC1vsPC2_colorByRegion.pdf", p, device = cairo_pdf(),
       width = 6, height = 4, useDingbats = F)

p <- ggbiplot(
  pca_tpm,
  pc.biplot = T,
  var.axes = F,
  groups = sample_meta$sex,
  ellipse = F,
  choices = c(1, 2)
) + theme_bw(base_size = 10) + 
  scale_color_manual(values = c("#377EB8","#4DAF4A","#E41A1C")) +
  theme(legend.position = "top", panel.grid.minor = element_blank()) +
  coord_cartesian()
p
ggsave("./plots/pca/PCA_allSamples_log10TPM_PC1vsPC2_colorBySex.pdf", p, device = cairo_pdf(),
       width = 6, height = 4, useDingbats = F)

## subsets
# neurons

neurons_subset <- sample_meta %>% filter(cell_type == "Neurons")
# neurons_subset <- sample_meta %>% filter(cell_type == "Oligodendrocytes")
# neurons_subset <- sample_meta %>% filter(cell_type == "Other_glias")

tpm_data_neurons_rm0 <- tpm_data %>%
  select(starts_with("gene"), one_of(neurons_subset$sample_id)) %>%
  mutate(RowSum = rowSums(select(., -starts_with("gene")))) %>% filter(RowSum>0)
t_tpmData_neurons <- tpm_data_neurons_rm0 %>% select(-starts_with("gene"), -RowSum)
t_tpmData_neurons <- log10(t_tpmData_neurons + 1) %>% t()

colnames(t_tpmData_neurons) <- tpm_data_neurons_rm0$geneID
pca_tpm_neurons <- prcomp(t_tpmData_neurons, center = T, scale. = T)
ggscreeplot(pca_tpm_neurons)

# check order
all(colnames(tpm_data_neurons_rm0)[c(4:41)]==neurons_subset$sample_id)

allPCcombn <- combn(6,2)
plist <- apply(allPCcombn,2,function(x){
  choices <- x
  p <- ggbiplot(
    pca_tpm_neurons,
    pc.biplot = T,
    var.axes = F,
    # groups = neurons_subset$sex,
    # groups = neurons_subset$region,
    groups = neurons_subset$disease,
    ellipse = T,
    choices = choices
  ) + theme_bw(base_size = 16) + 
    scale_color_manual(values = c("#377EB8","#4DAF4A","#E41A1C")) +
    # scale_color_brewer(palette = "Set1") +
    theme(legend.position = "top") +
    coord_cartesian()
  p
})
plist[["ncol"]] <- 3
do.call(grid.arrange, plist)

p <- ggbiplot(
  pca_tpm_neurons,
  pc.biplot = T,
  var.axes = F,
  groups = neurons_subset$disease,
  ellipse = F,
  choices = c(1, 2)
) + theme_bw(base_size = 10) + 
  scale_color_manual(values = c("#377EB8","#4DAF4A","#E41A1C")) +
  theme(legend.position = "top", panel.grid.minor = element_blank()) +
  coord_cartesian()
p
# check XIST

exp_XIST <- tpm_data %>% 
  filter(geneName == "XIST") %>% 
  gather(sample_id, TPM, -starts_with("gene")) %>% 
  left_join(sample_meta)

p <- exp_XIST %>% ggplot(aes(sex, log10(TPM +1)))

p + geom_point()

# check neuronal markers

exp_SNAP25 <- tpm_data %>% 
  # filter(geneName == "SNAP25") %>% 
  filter(geneName == "SOX10") %>% 
  gather(sample_id, TPM, -starts_with("gene")) %>% 
  left_join(sample_meta)

p <- exp_SNAP25 %>% ggplot(aes(cell_type, log10(TPM +1)))

p + geom_point()

# check mitochondria and ribosomal genes
library(ggrepel)

tpm_data_rm0 %>%
  select(-RowSum) %>% 
  filter(grepl("Mt_", geneType)) %>%
  gather(sample_id, TPM, -starts_with("gene")) %>% 
  left_join(sample_meta) %>% 
  # ggplot(aes(cell_type, log10(TPM+1))) +
  ggplot(aes(cell_type, TPM)) +
  geom_point() +
  geom_text_repel(aes(label = sample_id)) +
  facet_wrap(~geneName) +
  theme_bw(base_size = 10, base_family = "Helvetica")

tpm_data_rm0 %>%
  select(-RowSum) %>% 
  filter(grepl("rRNA", geneType)) %>%
  gather(sample_id, TPM, -starts_with("gene")) %>% 
  left_join(sample_meta) %>% 
  ggplot(aes(cell_type, log10(TPM+1))) +
  geom_point() +
  geom_text_repel(aes(label = sample_id), data = . %>% filter(TPM > 100)) +
  facet_wrap(~geneName, nrow = 3) +
  theme_bw(base_size = 10, base_family = "Helvetica")

exp_data_rm0 <- exp_data %>% mutate(RowSum = rowSums(select(., -starts_with("gene")))) %>% filter(RowSum>0)

exp_data_rm0 %>%
  select(-RowSum) %>% 
  filter(grepl("Mt_", geneType)) %>%
  gather(sample_id, raw_counts, -starts_with("gene")) %>% 
  left_join(sample_meta) %>% 
  # ggplot(aes(cell_type, raw_counts)) +
  ggplot(aes(cell_type, log10(raw_counts + 1))) +
  geom_point() +
  geom_text_repel(aes(label = sample_id), data = . %>% filter(raw_counts > 500)) +
  facet_wrap(~geneName) +
  theme_bw(base_size = 10, base_family = "Helvetica")

exp_data_rm0 %>%
  select(-RowSum) %>% 
  filter(grepl("rRNA", geneType)) %>%
  gather(sample_id, raw_counts, -starts_with("gene")) %>% 
  left_join(sample_meta) %>% 
  ggplot(aes(cell_type, log10(raw_counts + 1))) +
  geom_point() +
  geom_text_repel(aes(label = sample_id), data = . %>% filter(raw_counts > 500)) +
  facet_wrap(~geneName, scales = "free_y") +
  theme_bw(base_size = 10, base_family = "Helvetica")

## check ALS subtypes with marker genes in Tam et al. 2019 Cell Reports

als_subtype_markers <- read_tsv("als_subtypes_genes_Tam_2019_CellRep.txt")
als_samples <- sample_meta %>% 
  # filter(disease == "ALS", cell_type == "Neurons", region == "motor_cortex")
  # filter(disease == "ALS", cell_type == "Neurons", region == "mid_frontal_cortex")
  # filter(disease == "ALS", cell_type == "Other_glias")
  # filter(disease == "ALS", cell_type == "Neurons")
  filter(disease == "ALS", cell_type == "Oligodendrocytes")
  # filter(disease == "ALS")

tpm_data_markers_sub <- tpm_data %>%
  filter(geneName %in% als_subtype_markers$Gene) %>% 
  select(starts_with("gene"), one_of(als_samples$sample_id)) %>% 
  mutate(RowSum = rowSums(select(., -starts_with("gene")))) %>% 
  filter(RowSum>0)
  
t_tpmData_markers_sub <- tpm_data_markers_sub %>% select(-starts_with("gene"), -RowSum)
t_tpmData_markers_sub <- log10(t_tpmData_markers_sub + 1) %>% t()

colnames(t_tpmData_markers_sub) <- tpm_data_markers_sub$geneID
pca_tpm_markers_sub <- prcomp(t_tpmData_markers_sub, center = T, scale. = T)
ggscreeplot(pca_tpm_markers_sub)

allPCcombn <- combn(3,2)
plist <- apply(allPCcombn,2,function(x){
  choices <- x
  p <- ggbiplot(
    pca_tpm_markers_sub,
    pc.biplot = T,
    var.axes = F,
    groups = als_samples$subject_id,
    # groups = als_samples$cell_type,
    # groups = als_samples$region,
    # groups = als_samples$sex,
    ellipse = F,
    choices = choices
  ) + theme_bw(base_size = 16) + 
    # scale_color_manual(values = c("#377EB8","#4DAF4A","#E41A1C")) +
    # scale_color_brewer(palette = "Set1") +
    theme(legend.position = "top") +
    coord_cartesian()
  p
})
plist[["ncol"]] <- 3
do.call(grid.arrange, plist)

library(pheatmap)
library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

anno_df <- als_samples %>% 
  select(sample_id, subject_id, sex, region) %>% 
  as.data.frame() %>% 
  column_to_rownames("sample_id")

anno_gene_df <- tpm_data_markers_sub %>% 
  select(geneID, geneName) %>% 
  left_join(als_subtype_markers, by = c("geneName" = "Gene")) %>% 
  select(geneID, Subtype) %>% 
  as.data.frame() %>% 
  column_to_rownames("geneID")

ht_df <- tpm_data_markers_sub %>% 
  select(-geneName, -geneType, -RowSum) %>% 
  as.data.frame() %>% 
  column_to_rownames("geneID")

hclust_rows <- sort_hclust(hclust(dist(ht_df), method = "ward.D2"))
hclust_cols <- hclust(dist(t(ht_df)),method = "ward.D2")


pheatmap(ht_df,
         scale = "row",
         show_rownames = FALSE,
         annotation_col = anno_df,
         annotation_row = anno_gene_df,
         cluster_rows = hclust_rows,
         cluster_cols = hclust_cols)

pheatmap(ht_df,
         scale = "row",
         show_rownames = FALSE,
         annotation_col = anno_df,
         annotation_row = anno_gene_df)

# Correlation matrix

library(Hmisc)
library(ggthemes)

sample_meta_sorted <- sample_meta %>% arrange(cell_type, region, disease, sex)
tpm_cor_matrix <- tpm_data_rm0 %>% 
  select(one_of(sample_meta_sorted$sample_id)) %>% 
  as.matrix() %>% 
  rcorr(type = "spearman")

# remove mitochondria and ribosomal genes
tpm_cor_matrix <- tpm_data_rm0 %>% 
  filter(!(grepl("^Mt_", geneType) | grepl("rRNA", geneType))) %>%
  select(one_of(sample_meta_sorted$sample_id)) %>% 
  as.matrix() %>% 
  rcorr(type = "spearman")

tpm_cor_df <- as_data_frame(tpm_cor_matrix$r)

# log TPM
tpm_data_rm0_log <- tpm_data_rm0 %>% select(one_of(sample_meta_sorted$sample_id)) %>% as.matrix()
# remove mitochondria and ribosomal genes
tpm_data_rm0_log <- tpm_data_rm0 %>% 
  filter(!(grepl("^Mt_", geneType) | grepl("rRNA", geneType))) %>%
  select(one_of(sample_meta_sorted$sample_id)) %>% 
  as.matrix()

tpm_data_rm0_log <- log10(tpm_data_rm0_log+1)
tpm_data_rm0_log %<>% rcorr(type="spearman")
tpm_cor_df <- as_data_frame(tpm_data_rm0_log$r)

tpm_cor_gg <- tpm_cor_df %>%
  mutate(sample = rownames(tpm_cor_matrix$r)) %>%
  select(sample, everything()) %>%
  gather(sampleY, cor, -sample) %>%
  mutate(sample = factor(sample, levels = sample_meta_sorted$sample_id),
         sampleY = factor(sampleY, levels = sample_meta_sorted$sample_id))

p <- ggplot(tpm_cor_gg, aes(sample, sampleY, fill = cor))


p + geom_tile() +
  scale_fill_viridis(name = "Spearman Correlaion",
                     direction = 1,
                     option = "D") +
  coord_equal() + labs(x = NULL, y = NULL) + theme_tufte(base_family = "Helvetica", base_size = 5) +
  theme(
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.spacing.x = unit(0.5, "cm"),
    panel.spacing.y = unit(0.5, "cm")
  )
ggsave(filename = "./plots/pca/correlationMatrix_allSamples_log10TPM_spearman.pdf", width = 8, height = 6)


## edgeR

## run edgeR separately for MCX and mFCX, for each cell type

selected_region <- "motor_cortex"
# selected_region <- "mid_frontal_cortex"
selected_celltype <- "Neurons"
# selected_celltype <- "Oligodendrocytes"
# selected_celltype <- "Other_glias"
# selected_sex <- "M"
selected_sex <- "F"

# sample_meta_subset <- sample_meta_sorted %>% filter(region == selected_region, cell_type == selected_celltype)

## remove outlier samples
# outlier_samples <- c("54FN", "91FM", "54FO", "54FM", "674MO")
outlier_samples <- c("54FN", "91FM", "54FO", "54FM", "674MO", "61FN", "36FN", "36MN", "91FO", "56MO")
sample_meta_subset <- sample_meta_sorted %>%
  filter(region == selected_region, cell_type == selected_celltype, !sample_id %in% outlier_samples)

# keep all genes
# exp_data_subset <- exp_data

# select chosen samples, filter out MT-genes and rRNA
# exp_data_subset <- exp_data %>%
#   filter(!(grepl("^MT-", geneName) | geneType == "rRNA"))

# only keep protein-coding genes? 
exp_data_subset <- exp_data %>%
  filter(geneType == "protein_coding") %>% 
  filter(!grepl("^MT-", geneName))

exp_data_subset_df <- exp_data_subset %>% 
  select(one_of(sample_meta_subset$sample_id)) %>%
  mutate_all(funs(as.integer(round(.)))) %>%  # round raw counts
  as.data.frame()

rownames(exp_data_subset_df) <- exp_data_subset$geneID

y <- DGEList(counts = exp_data_subset_df, group = sample_meta_subset$disease)
y
y$samples
# keep <- rowSums(cpm(y)>2) >= 6 # CPM threshold depends on the smallest library, sample number threshold is the number of samples in the smallest group
keep <- rowSums(cpm(y)>20) >= 6 # CPM threshold depends on the smallest library, sample number threshold is the number of samples in the smallest group
y <- y[keep, , keep.lib.size = FALSE]
y$samples
y <- calcNormFactors(y, method = "TMM")
y$samples


## GLM mode
# do not adjust for sex differences
design <- model.matrix(~ 0 + factor(sample_meta_subset$disease, levels = c("Control", "ALS", "FTD")))
colnames(design) <- c("Control", "ALS", "FTD")
design

y <- estimateDisp(y, design, robust = T)
fit <- glmQLFit(y, design, robust = T)
y$common.dispersion
plotBCV(y)
plotQLDisp(fit)

contrast <- makeContrasts(ALS_vs_Control = ALS - Control,
                          FTD_vs_Control = FTD - Control,
                          FTD_vs_ALS = FTD - ALS,
                          levels = design)

CPM_filtered <- cpm(y) %>% as_data_frame() %>% 
  rename_all(funs(paste0("CPM_", .))) %>%
  mutate(geneID = rownames(y$counts)) %>%
  select(geneID, everything())

# ALS vs Control
qlf_ALS_vs_Control <- glmQLFTest(fit, contrast = contrast[, "ALS_vs_Control"])
tt_ALS_vs_Control <- topTags(qlf_ALS_vs_Control, n=dim(y)[1])
results_qlf_ALS_vs_Control <- tt_ALS_vs_Control$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_ALS_vs_Control %>% filter(FDR<0.05) %>% dim()
  

# FTD vs Control
qlf_FTD_vs_Control <- glmQLFTest(fit, contrast = contrast[, "FTD_vs_Control"])
tt_FTD_vs_Control <- topTags(qlf_FTD_vs_Control, n=dim(y)[1])
results_qlf_FTD_vs_Control <- tt_FTD_vs_Control$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_FTD_vs_Control %>% filter(FDR<0.05) %>% dim()
  

# FTD vs ALS
qlf_FTD_vs_ALS <- glmQLFTest(fit, contrast = contrast[, "FTD_vs_ALS"])
tt_FTD_vs_ALS <- topTags(qlf_FTD_vs_ALS, n=dim(y)[1])
results_qlf_FTD_vs_ALS <- tt_FTD_vs_ALS$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_FTD_vs_ALS %>% filter(FDR<0.05) %>% dim()


 # model with samples in comparison only

selected_region <- "motor_cortex"
# selected_region <- "mid_frontal_cortex"
selected_celltype <- "Neurons"
# selected_celltype <- "Oligodendrocytes"
# selected_celltype <- "Other_glias"
# selected_sex <- "M"
selected_sex <- "F"

run_edgeR <- function(selected_region,
                      selected_celltype, 
                      selected_disease,
                      keep_outlier = FALSE,
                      pc_genes_only = FALSE,
                      de_mode = "glm"){
  
  sample_meta_subset <- sample_meta_sorted %>%
    filter(region == selected_region, 
           cell_type == selected_celltype,
           disease %in% selected_disease)
  
  ## remove outlier samples
  if(!keep_outlier){
    # outlier_samples <- c("54FN", "91FM", "54FO", "54FM", "674MO")
    outlier_samples <- c("54FN", "91FM", "54FO", "54FM", "674MO", "61FN", "36FN", "36MN", "91FO", "56MO")
    sample_meta_subset <- sample_meta_subset %>% filter(!sample_id %in% outlier_samples)
  }
  
  if(pc_genes_only){
    # only keep protein-coding genes? 
    exp_data_subset <- exp_data %>%
      filter(geneType == "protein_coding") %>% 
      filter(!grepl("^MT-", geneName))
  } else{
    # keep all genes
    exp_data_subset <- exp_data
    # select chosen samples, filter out MT-genes and rRNA
    # exp_data_subset <- exp_data %>%
    #   filter(!(grepl("^MT-", geneName) | geneType == "rRNA"))
  }
  
  
  exp_data_subset_df <- exp_data_subset %>% 
    select(one_of(sample_meta_subset$sample_id)) %>%
    mutate_all(funs(as.integer(round(.)))) %>%  # round raw counts
    as.data.frame()
  
  rownames(exp_data_subset_df) <- exp_data_subset$geneID
  
  y <- DGEList(counts = exp_data_subset_df, group = sample_meta_subset$disease)
  # print(y)
  # print(y$samples)
  keep <- rowSums(cpm(y)>2) >= 6 # CPM threshold depends on the smallest library, sample number threshold is the number of samples in the smallest group
  y <- y[keep, , keep.lib.size = FALSE]
  # print(y$samples)
  y <- calcNormFactors(y, method = "TMM")
  # print(y$samples)
  
  if(de_mode == "glm"){
    ## GLM mode
    # do not adjust for sex differences
    design <- model.matrix(~ 0 + factor(sample_meta_subset$disease, levels = selected_disease))
    colnames(design) <- selected_disease
    # print(design)
    
    y <- estimateDisp(y, design, robust = T)
    fit <- glmQLFit(y, design, robust = T)
    # print(y$common.dispersion)
    # plotBCV(y)
    # plotQLDisp(fit)
    contrast_args <- list(comp = paste0(selected_disease[2], "-", selected_disease[1]),
                          levels = design)
    contrast <- do.call(makeContrasts, contrast_args)
    # print(contrast)
    
    CPM_filtered <- cpm(y) %>% as_data_frame() %>% 
      rename_all(funs(paste0("CPM_", .))) %>%
      mutate(geneID = rownames(y$counts)) %>%
      select(geneID, everything())
    
    qlf <- glmQLFTest(fit, contrast = contrast[, "comp"])
    tt <- topTags(qlf, n=dim(y)[1])
    results <- tt$table %>% 
      rownames_to_column("geneID") %>% 
      as_tibble() %>%
      left_join(gene_annotation) %>% 
      select(geneID, geneName, geneType, everything()) %>%
      left_join(tpm_data %>% 
                  select(geneID, one_of(sample_meta_subset$sample_id)) %>% 
                  rename_if(is.double, funs(paste0("TPM_", .)))) %>%
      left_join(CPM_filtered)
    
  }
  
  if (de_mode == "exact"){
    y <- estimateDisp(y, robust = T)
    et <- exactTest(y, pair = selected_disease)
    tt <- topTags(et, n = dim(y)[1])
    
    CPM_filtered <- cpm(y) %>% as_data_frame() %>% 
      rename_all(funs(paste0("CPM_", .))) %>%
      mutate(geneID = rownames(y$counts)) %>%
      select(geneID, everything())
    
    results <- tt$table %>%
      rownames_to_column("geneID") %>% 
      as_tibble() %>%
      left_join(gene_annotation) %>% 
      select(geneID, geneName, geneType, everything()) %>%
      left_join(tpm_data %>% 
                  select(geneID, one_of(sample_meta_subset$sample_id)) %>% 
                  rename_if(is.double, funs(paste0("TPM_", .)))) %>%
      left_join(CPM_filtered)
    
  }
  
    # results %>%
    #   filter(FDR<0.05)
  
    # results %>%
    #   filter(FDR<0.05) %>%
    #   dim()
  
    # results %>% 
    #   mutate(comparison = paste0(selected_disease[2], "_vs_", selected_disease[1]),
    #          region = selected_region,
    #          cell_type = selected_celltype) %>% 
    #   select(comparison, region, cell_type, geneID:FDR)
    
    results %>% 
      mutate(comparison = paste0(selected_disease[2], "_vs_", selected_disease[1]),
             region = selected_region,
             cell_type = selected_celltype) %>% 
      select(geneID, geneName, logFC, logCPM, PValue, FDR,
             comparison, region, cell_type, starts_with("CPM"))
} 

# test genes found in snRNA
all_diseases <- c("ALS", "FTD", "Control")
disease_color_panel <- c("#f2a034", "#bf36ff", "#36bcff")
names(disease_color_panel) <- all_diseases

library(ggrepel)

tmp <- run_edgeR("motor_cortex","Neurons", c("Control", "ALS"), F, F)
tmp <- run_edgeR("mid_frontal_cortex","Neurons", c("Control", "ALS"), F, F)
tmp <- run_edgeR("motor_cortex","Other_glias", c("Control", "ALS"), F, F)
tmp <- run_edgeR("mid_frontal_cortex","Other_glias", c("Control", "ALS"), F, F)

test_gene <- function(selected_gene, selected_disease){
  tmp_df <- tmp %>% 
    filter(geneName == selected_gene)
  tmp_df_long <- tmp_df %>% 
    gather(sample_id, CPM, starts_with("CPM")) %>% 
    mutate(sample_id = str_replace(sample_id, "CPM_", "")) %>% 
    left_join(sample_meta_sorted) 
  p <- tmp_df_long %>% 
    ggplot(aes(disease, CPM))
  
  p + geom_boxplot(aes(color = disease), width = 0.25) +
    geom_point(aes(color = disease)) +
    geom_text_repel(aes(label = subject_id)) +
    scale_color_manual(values = disease_color_panel[selected_disease]) +
    ggtitle(paste(selected_gene, tmp_df$region[1], tmp_df$cell_type[1],
                  paste0("log2FC: ", round(tmp_df$logFC, 3)), sep = ", ")) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())
  ggsave(filename = paste0("./interest_genes/FANS_sorted_RNA_", selected_gene, "_",
                           tmp_df$cell_type[1], "_", tmp_df$region[1], "_",
                           selected_disease[1], "_vs_", selected_disease[2], 
                           "_individual_boxplot.pdf"),
         width = 6, height = 4, device = cairo_pdf(), useDingbats = FALSE)
}

# neu
test_gene("HSP90AA1", c("ALS", "Control"))
test_gene("HSP90AB1", c("ALS", "Control"))
test_gene("HSPH1", c("ALS", "Control"))
test_gene("CALM1", c("ALS", "Control"))
test_gene("CALM3", c("ALS", "Control"))
test_gene("DYNLL1", c("ALS", "Control"))
test_gene("KCNB2", c("ALS", "Control"))
test_gene("KCND3", c("ALS", "Control"))
test_gene("KCNQ3", c("ALS", "Control"))

# MgAs
test_gene("GFAP", c("ALS", "Control"))
test_gene("CD44", c("ALS", "Control"))
test_gene("CHI3L1", c("ALS", "Control"))
test_gene("SLC1A2", c("ALS", "Control"))
test_gene("SLC1A3", c("ALS", "Control"))
test_gene("CACNB2", c("ALS", "Control"))
test_gene("NRP1", c("ALS", "Control"))
test_gene("ERBB4", c("ALS", "Control"))
test_gene("NFIB", c("ALS", "Control"))
test_gene("EDN1", c("ALS", "Control"))
test_gene("SLC1A1", c("ALS", "Control"))
test_gene("ACTN4", c("ALS", "Control"))
test_gene("CTGF", c("ALS", "Control"))
test_gene("ID1", c("ALS", "Control"))
test_gene("PLCG2", c("ALS", "Control"))
test_gene("SLC2A3", c("ALS", "Control"))
test_gene("TNS1", c("ALS", "Control"))


#
res_merged <- bind_rows(
  run_edgeR("motor_cortex","Neurons", c("Control", "ALS"), F, F),
  run_edgeR("motor_cortex","Oligodendrocytes", c("Control", "ALS"), F, F),
  run_edgeR("motor_cortex","Other_glias", c("Control", "ALS"), F, F),
  
  run_edgeR("motor_cortex","Neurons", c("Control", "FTD"), F, F),
  run_edgeR("motor_cortex","Oligodendrocytes", c("Control", "FTD"), F, F),
  run_edgeR("motor_cortex","Other_glias", c("Control", "FTD"), F, F),
  
  run_edgeR("motor_cortex","Neurons", c("FTD", "ALS"), F, F),
  run_edgeR("motor_cortex","Oligodendrocytes", c("FTD", "ALS"), F, F),
  run_edgeR("motor_cortex","Other_glias", c("FTD", "ALS"), F, F),
  
  run_edgeR("mid_frontal_cortex","Neurons", c("Control", "ALS"), F, F),
  run_edgeR("mid_frontal_cortex","Oligodendrocytes", c("Control", "ALS"), F, F),
  run_edgeR("mid_frontal_cortex","Other_glias", c("Control", "ALS"), F, F),
  
  run_edgeR("mid_frontal_cortex","Neurons", c("Control", "FTD"), F, F),
  run_edgeR("mid_frontal_cortex","Oligodendrocytes", c("Control", "FTD"), F, F),
  run_edgeR("mid_frontal_cortex","Other_glias", c("Control", "FTD"), F, F),
  
  run_edgeR("mid_frontal_cortex","Neurons", c("FTD", "ALS"), F, F),
  run_edgeR("mid_frontal_cortex","Oligodendrocytes", c("FTD", "ALS"), F, F),
  run_edgeR("mid_frontal_cortex","Other_glias", c("FTD", "ALS"), F, F)
)

write_tsv(res_merged, "edgeR_res_glmTest_allExpressedGenes_rmOutlierSamples.txt")

#
run_edgeR("motor_cortex","Neurons", c("Control", "ALS"), F, T)
run_edgeR("motor_cortex","Oligodendrocytes", c("Control", "ALS"), F, T)
run_edgeR("motor_cortex","Other_glias", c("Control", "ALS"), F, T)

run_edgeR("motor_cortex","Neurons", c("Control", "FTD"), F, T)
run_edgeR("motor_cortex","Oligodendrocytes", c("Control", "FTD"), F, T)
run_edgeR("motor_cortex","Other_glias", c("Control", "FTD"), F, T)

run_edgeR("motor_cortex","Neurons", c("FTD", "ALS"), F, T)
run_edgeR("motor_cortex","Oligodendrocytes", c("FTD", "ALS"), F, T)
run_edgeR("motor_cortex","Other_glias", c("FTD", "ALS"), F, T)

run_edgeR("mid_frontal_cortex","Neurons", c("Control", "ALS"), F, T)
run_edgeR("mid_frontal_cortex","Oligodendrocytes", c("Control", "ALS"), F, T)
run_edgeR("mid_frontal_cortex","Other_glias", c("Control", "ALS"), F, T)

run_edgeR("mid_frontal_cortex","Neurons", c("Control", "FTD"), F, T)
run_edgeR("mid_frontal_cortex","Oligodendrocytes", c("Control", "FTD"), F, T)
run_edgeR("mid_frontal_cortex","Other_glias", c("Control", "FTD"), F, T)

run_edgeR("mid_frontal_cortex","Neurons", c("FTD", "ALS"), F, T)
run_edgeR("mid_frontal_cortex","Oligodendrocytes", c("FTD", "ALS"), F, T)
run_edgeR("mid_frontal_cortex","Other_glias", c("FTD", "ALS"), F, T)

run_edgeR("motor_cortex","Neurons", c("Control", "ALS"), F, T, de_mode = "exact")
run_edgeR("motor_cortex","Oligodendrocytes", c("Control", "ALS"), F, T, de_mode = "exact")
run_edgeR("motor_cortex","Other_glias", c("Control", "ALS"), F, T, de_mode = "exact")

run_edgeR("motor_cortex","Neurons", c("Control", "FTD"), F, T, de_mode = "exact")
run_edgeR("motor_cortex","Oligodendrocytes", c("Control", "FTD"), F, T, de_mode = "exact")
run_edgeR("motor_cortex","Other_glias", c("Control", "FTD"), F, T, de_mode = "exact")

run_edgeR("motor_cortex","Neurons", c("FTD", "ALS"), F, T, de_mode = "exact")
run_edgeR("motor_cortex","Oligodendrocytes", c("FTD", "ALS"), F, T, de_mode = "exact")
run_edgeR("motor_cortex","Other_glias", c("FTD", "ALS"), F, T, de_mode = "exact")

run_edgeR("mid_frontal_cortex","Neurons", c("Control", "ALS"), F, T, de_mode = "exact")
run_edgeR("mid_frontal_cortex","Oligodendrocytes", c("Control", "ALS"), F, T, de_mode = "exact")
run_edgeR("mid_frontal_cortex","Other_glias", c("Control", "ALS"), F, T, de_mode = "exact")

run_edgeR("mid_frontal_cortex","Neurons", c("Control", "FTD"), F, T, de_mode = "exact")
run_edgeR("mid_frontal_cortex","Oligodendrocytes", c("Control", "FTD"), F, T, de_mode = "exact")
run_edgeR("mid_frontal_cortex","Other_glias", c("Control", "FTD"), F, T, de_mode = "exact")

run_edgeR("mid_frontal_cortex","Neurons", c("FTD", "ALS"), F, T, de_mode = "exact")
run_edgeR("mid_frontal_cortex","Oligodendrocytes", c("FTD", "ALS"), F, T, de_mode = "exact")
run_edgeR("mid_frontal_cortex","Other_glias", c("FTD", "ALS"), F, T, de_mode = "exact")
 

# adjusting for sex differences
design <- model.matrix(~ factor(sample_meta_subset$sex) + factor(sample_meta_subset$disease, levels = c("Control", "ALS", "FTD")))
design

y <- estimateDisp(y, design, robust = T)
fit <- glmQLFit(y, design, robust = T)
y$common.dispersion
plotBCV(y)
plotQLDisp(fit)

CPM_filtered <- cpm(y) %>% as_data_frame() %>% 
  rename_all(funs(paste0("CPM_", .))) %>%
  mutate(geneID = rownames(y$counts)) %>%
  select(geneID, everything())

# ALS vs Control
qlf_ALS_vs_Control <- glmQLFTest(fit, coef = 3)
tt_ALS_vs_Control <- topTags(qlf_ALS_vs_Control, n=dim(y)[1])
results_qlf_ALS_vs_Control <- tt_ALS_vs_Control$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_ALS_vs_Control %>% filter(FDR<0.05) %>% dim()
  

# FTD vs Control
qlf_FTD_vs_Control <- glmQLFTest(fit, coef = 4)
tt_FTD_vs_Control <- topTags(qlf_FTD_vs_Control, n=dim(y)[1])
results_qlf_FTD_vs_Control <- tt_FTD_vs_Control$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_FTD_vs_Control %>% filter(FDR<0.05) %>% dim()
  

# FTD vs ALS
qlf_FTD_vs_ALS <- glmQLFTest(fit, contrast = c(0, 0, -1, 1))
tt_FTD_vs_ALS <- topTags(qlf_FTD_vs_ALS, n=dim(y)[1])
results_qlf_FTD_vs_ALS <- tt_FTD_vs_ALS$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_FTD_vs_ALS %>% filter(FDR<0.05) %>% dim()
  

# separating male and female samples

sample_meta_subset <- sample_meta_sorted %>% 
  filter(region == selected_region,
         cell_type == selected_celltype,
         sex == selected_sex)

# keep all genes
exp_data_subset <- exp_data

# select chosen samples, filter out MT-genes and rRNA
# exp_data_subset <- exp_data %>%
#   filter(!(grepl("^MT-", geneName) | geneType == "rRNA"))
exp_data_subset_df <- exp_data_subset %>% 
  select(one_of(sample_meta_subset$sample_id)) %>%
  mutate_all(funs(as.integer(round(.)))) %>%  # round raw counts
  as.data.frame()

rownames(exp_data_subset_df) <- exp_data_subset$geneID

y <- DGEList(counts = exp_data_subset_df, group = sample_meta_subset$disease)
y
y$samples
keep <- rowSums(cpm(y)>2) >= 3 # CPM threshold depends on the smallest library, sample number threshold is the number of samples in the smallest group 
y <- y[keep, , keep.lib.size = FALSE]
y$samples
y <- calcNormFactors(y, method = "TMM")
y$samples


## GLM mode
design <- model.matrix(~ 0 + factor(sample_meta_subset$disease, levels = c("Control", "ALS", "FTD")))
colnames(design) <- c("Control", "ALS", "FTD")
design

y <- estimateDisp(y, design, robust = T)
fit <- glmQLFit(y, design, robust = T)
y$common.dispersion
plotBCV(y)
plotQLDisp(fit)

contrast <- makeContrasts(ALS_vs_Control = ALS - Control,
                          FTD_vs_Control = FTD - Control,
                          FTD_vs_ALS = FTD - ALS,
                          levels = design)

CPM_filtered <- cpm(y) %>% as_data_frame() %>% 
  rename_all(funs(paste0("CPM_", .))) %>%
  mutate(geneID = rownames(y$counts)) %>%
  select(geneID, everything())

# ALS vs Control
qlf_ALS_vs_Control <- glmQLFTest(fit, contrast = contrast[, "ALS_vs_Control"])
tt_ALS_vs_Control <- topTags(qlf_ALS_vs_Control, n=dim(y)[1])
results_qlf_ALS_vs_Control <- tt_ALS_vs_Control$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_ALS_vs_Control %>% filter(FDR<0.05) %>% dim()
  

# FTD vs Control
qlf_FTD_vs_Control <- glmQLFTest(fit, contrast = contrast[, "FTD_vs_Control"])
tt_FTD_vs_Control <- topTags(qlf_FTD_vs_Control, n=dim(y)[1])
results_qlf_FTD_vs_Control <- tt_FTD_vs_Control$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_FTD_vs_Control %>% filter(FDR<0.05) %>% dim()
  

# FTD vs ALS
qlf_FTD_vs_ALS <- glmQLFTest(fit, contrast = contrast[, "FTD_vs_ALS"])
tt_FTD_vs_ALS <- topTags(qlf_FTD_vs_ALS, n=dim(y)[1])
results_qlf_FTD_vs_ALS <- tt_FTD_vs_ALS$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_FTD_vs_ALS %>% filter(FDR<0.05) %>% dim()
  


# write_tsv(results.qlf.GLUvsGABA, "edgeR_glmQLFTest_GLUvsGABA_DEgenes.txt")






## sanity check, testing DE on cell types in control samples

selected_region <- "motor_cortex"
# selected_region <- "mid_frontal_cortex"
selected_disease <- "Control"
# selected_disease <- "ALS"
# selected_disease <- "FTD"
# selected_celltype <- "Neurons"
# selected_celltype <- "Oligodendrocytes"
# selected_celltype <- "Other_glias"

sample_meta_subset <- sample_meta_sorted %>% filter(region == selected_region, disease == selected_disease)

# keep all genes
exp_data_subset <- exp_data

# select chosen samples, filter out MT-genes and rRNA
# exp_data_subset <- exp_data %>%
#   filter(!(grepl("^MT-", geneName) | geneType == "rRNA"))
exp_data_subset_df <- exp_data_subset %>% 
  select(one_of(sample_meta_subset$sample_id)) %>%
  mutate_all(funs(as.integer(round(.)))) %>%  # round raw counts
  as.data.frame()

rownames(exp_data_subset_df) <- exp_data_subset$geneID

y <- DGEList(counts = exp_data_subset_df, group = sample_meta_subset$cell_type)
y
y$samples
keep <- rowSums(cpm(y)>2) >= 6 # CPM threshold depends on the smallest library, sample number threshold is the number of samples in the smallest group 
y <- y[keep, , keep.lib.size = FALSE]
y$samples
y <- calcNormFactors(y, method = "TMM")
y$samples


## GLM mode
design <- model.matrix(~ 0 + factor(sample_meta_subset$cell_type, levels = c("Neurons", "Oligodendrocytes", "Other_glias")))
colnames(design) <- c("Neurons", "Oligodendrocytes", "Other_glias")
design

y <- estimateDisp(y, design, robust = T)
fit <- glmQLFit(y, design, robust = T)
y$common.dispersion
plotBCV(y)
plotQLDisp(fit)

contrast <- makeContrasts(N_vs_O = Neurons - Oligodendrocytes,
                          N_vs_G = Neurons - Other_glias,
                          O_vs_G = Oligodendrocytes - Other_glias,
                          levels = design)

CPM_filtered <- cpm(y) %>% as_data_frame() %>% 
  rename_all(funs(paste0("CPM_", .))) %>%
  mutate(geneID = rownames(y$counts)) %>%
  select(geneID, everything())

# Neurons vs Oligodendrocytes
qlf_Neurons_vs_Oligodendrocytes <- glmQLFTest(fit, contrast = contrast[, "N_vs_O"])
tt_Neurons_vs_Oligodendrocytes <- topTags(qlf_Neurons_vs_Oligodendrocytes, n=dim(y)[1])
results_qlf_Neurons_vs_Oligodendrocytes <- tt_Neurons_vs_Oligodendrocytes$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_Neurons_vs_Oligodendrocytes %>% filter(FDR<0.05) %>% dim()

# Neurons vs Other glias
qlf_Neurons_vs_OtherGlias <- glmQLFTest(fit, contrast = contrast[, "N_vs_G"])
tt_Neurons_vs_OtherGlias <- topTags(qlf_Neurons_vs_OtherGlias, n=dim(y)[1])
results_qlf_Neurons_vs_OtherGlias <- tt_Neurons_vs_OtherGlias$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_Neurons_vs_OtherGlias %>% filter(FDR<0.05) %>% dim()
  
# Oligodendrocytes vs Other glias
qlf_Oligodendrocytes_vs_OtherGlias <- glmQLFTest(fit, contrast = contrast[, "O_vs_G"])
tt_Oligodendrocytes_vs_OtherGlias <- topTags(qlf_Oligodendrocytes_vs_OtherGlias, n=dim(y)[1])
results_qlf_Oligodendrocytes_vs_OtherGlias <- tt_Oligodendrocytes_vs_OtherGlias$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_Oligodendrocytes_vs_OtherGlias %>% filter(FDR<0.05) %>% dim()
  

## sanity check, testing DE between replicates/individuals

# selected_region <- "motor_cortex"
selected_region <- "mid_frontal_cortex"
selected_disease <- "Control"
# selected_disease <- "ALS"
# selected_disease <- "FTD"
selected_celltype <- "Neurons"
# selected_celltype <- "Oligodendrocytes"
# selected_celltype <- "Other_glias"

sample_meta_subset <- sample_meta_sorted %>% 
  filter(region == selected_region, 
         disease == selected_disease,
         cell_type == selected_celltype)

# keep all genes
exp_data_subset <- exp_data

# select chosen samples, filter out MT-genes and rRNA
# exp_data_subset <- exp_data %>%
#   filter(!(grepl("^MT-", geneName) | geneType == "rRNA"))

exp_data_subset_df <- exp_data_subset %>% 
  select(one_of(sample_meta_subset$sample_id)) %>%
  mutate_all(funs(as.integer(round(.)))) %>%  # round raw counts
  as.data.frame()

rownames(exp_data_subset_df) <- exp_data_subset$geneID

y <- DGEList(counts = exp_data_subset_df, group = sample_meta_subset$sex)
y
y$samples
keep <- rowSums(cpm(y)>2) >= 3 # CPM threshold depends on the smallest library, sample number threshold is the number of samples in the smallest group 
y <- y[keep, , keep.lib.size = FALSE]
y$samples
y <- calcNormFactors(y, method = "TMM")
y$samples


## GLM mode
design <- model.matrix(~ 0 + factor(sample_meta_subset$sex, levels = c("M", "F")))
colnames(design) <- c("M", "F")
design

y <- estimateDisp(y, design, robust = T)
fit <- glmQLFit(y, design, robust = T)
y$common.dispersion
plotBCV(y)
plotQLDisp(fit)

contrast <- makeContrasts(M_vs_F = M - F,
                          levels = design)

CPM_filtered <- cpm(y) %>% as_data_frame() %>% 
  rename_all(funs(paste0("CPM_", .))) %>%
  mutate(geneID = rownames(y$counts)) %>%
  select(geneID, everything())

# Male vs Female
qlf_Male_vs_Female <- glmQLFTest(fit, contrast = contrast[, "M_vs_F"])
tt_Male_vs_Female <- topTags(qlf_Male_vs_Female, n=dim(y)[1])
results_qlf_Male_vs_Female <- tt_Male_vs_Female$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_qlf_Male_vs_Female %>% filter(FDR<0.05) %>% dim()


# try the glmLRT framework

selected_region <- "motor_cortex"
# selected_region <- "mid_frontal_cortex"
selected_celltype <- "Neurons"
# selected_celltype <- "Oligodendrocytes"
# selected_celltype <- "Other_glias"
# selected_sex <- "M"
selected_sex <- "F"

sample_meta_subset <- sample_meta_sorted %>% filter(region == selected_region, cell_type == selected_celltype)

# keep all genes
exp_data_subset <- exp_data

# select chosen samples, filter out MT-genes and rRNA
# exp_data_subset <- exp_data %>%
#   filter(!(grepl("^MT-", geneName) | geneType == "rRNA"))

exp_data_subset_df <- exp_data_subset %>% 
  select(one_of(sample_meta_subset$sample_id)) %>%
  mutate_all(funs(as.integer(round(.)))) %>%  # round raw counts
  as.data.frame()

rownames(exp_data_subset_df) <- exp_data_subset$geneID

y <- DGEList(counts = exp_data_subset_df, group = sample_meta_subset$disease)
y
y$samples
keep <- rowSums(cpm(y)>2) >= 6 # CPM threshold depends on the smallest library, sample number threshold is the number of samples in the smallest group
# keep <- rowSums(cpm(y)>5) >= 6 # CPM threshold depends on the smallest library, sample number threshold is the number of samples in the smallest group 
y <- y[keep, , keep.lib.size = FALSE]
y$samples
y <- calcNormFactors(y, method = "TMM")
y$samples


## GLM mode, LRT
# do not adjust for sex differences
design <- model.matrix(~ 0 + factor(sample_meta_subset$disease, levels = c("Control", "ALS", "FTD")))
colnames(design) <- c("Control", "ALS", "FTD")
design

y <- estimateDisp(y, design, robust = T)
fit <- glmFit(y, design, robust = T)
y$common.dispersion
plotBCV(y)

contrast <- makeContrasts(ALS_vs_Control = ALS - Control,
                          FTD_vs_Control = FTD - Control,
                          FTD_vs_ALS = FTD - ALS,
                          levels = design)

CPM_filtered <- cpm(y) %>% as_data_frame() %>% 
  rename_all(funs(paste0("CPM_", .))) %>%
  mutate(geneID = rownames(y$counts)) %>%
  select(geneID, everything())

# ALS vs Control
lrt_ALS_vs_Control <- glmLRT(fit, contrast = contrast[, "ALS_vs_Control"])
tt_ALS_vs_Control <- topTags(lrt_ALS_vs_Control, n=dim(y)[1])
results_lrt_ALS_vs_Control <- tt_ALS_vs_Control$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_lrt_ALS_vs_Control %>% filter(FDR<0.05) %>% dim()
  

# FTD vs Control
lrt_FTD_vs_Control <- glmLRT(fit, contrast = contrast[, "FTD_vs_Control"])
tt_FTD_vs_Control <- topTags(lrt_FTD_vs_Control, n=dim(y)[1])
results_lrt_FTD_vs_Control <- tt_FTD_vs_Control$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_lrt_FTD_vs_Control %>% filter(FDR<0.05) %>% dim()
  

# FTD vs ALS
lrt_FTD_vs_ALS <- glmLRT(fit, contrast = contrast[, "FTD_vs_ALS"])
tt_FTD_vs_ALS <- topTags(lrt_FTD_vs_ALS, n=dim(y)[1])
results_lrt_FTD_vs_ALS <- tt_FTD_vs_ALS$table %>% 
  rownames_to_column("geneID") %>% 
  as_tibble() %>%
  left_join(gene_annotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpm_data %>% select(geneID, one_of(sample_meta_subset$sample_id)) %>% rename_if(is.double, funs(paste0("TPM_", .)))) %>%
  left_join(CPM_filtered)

results_lrt_FTD_vs_ALS %>% filter(FDR<0.05) %>% dim()
  







## ANOVA-like test
## need intercept term in the design formula?

anov <- glmQLFTest(fit, contrast = contrast)
tt.anov <- topTags(anov, n=dim(y)[1])
results.anov <- tt.anov$table %>%
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>%
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GLU"), ends_with("GABA"), ends_with("OLIG"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GLU"), ends_with("GABA"), ends_with("OLIG")))
write_tsv(results.anov, "edgeR_glmQLFTest_anova_DEgenes.txt")


## Control for batch effect
design <- model.matrix(~sampleInfo$batch+sampleInfo$cellType)
# colnames(design) <- c("GABA","GLU","OLIG")
rownames(design) <- sampleInfo$sample
design

y <- DGEList(counts = expData.sub, group = sampleInfo$cellType)
keep <- rowSums(cpm(y)>2) >= 9 
y <- y[keep, , keep.lib.size = FALSE]
y <- calcNormFactors(y, method = "TMM")

y <- estimateDisp(y, design, robust = T)
fit <- glmQLFit(y, design, robust = T, coef = 5:6)

CPM.filtered <- cpm(y) %>% as_data_frame() %>% 
  rename_(.dots = setNames(names(.), paste0("CPM_", names(.)))) %>%
  mutate(geneID = rownames(y$counts)) %>%
  select(geneID, everything())

# GLU vs GABA
qlf.GLUvsGABA.controlBatch <- glmQLFTest(fit, coef = 5)
tt.GLUvsGABA.controlBatch <- topTags(qlf.GLUvsGABA.controlBatch, n=dim(y)[1])
results.qlf.GLUvsGABA.controlBatch <- tt.GLUvsGABA.controlBatch$table %>% 
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GLU"), ends_with("GABA"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GLU"), ends_with("GABA")))
  
write_tsv(results.qlf.GLUvsGABA.controlBatch, "edgeR_glmQLFTest_GLUvsGABA_DEgenes_controlBatch.txt")

# GLU vs OLIG
qlf.GLUvsOLIG.controlBatch <- glmQLFTest(fit, contrast = c(0, 0, 0, 0, 1, -1))
tt.GLUvsOLIG.controlBatch <- topTags(qlf.GLUvsOLIG.controlBatch, n=dim(y)[1])
results.qlf.GLUvsOLIG.controlBatch <- tt.GLUvsOLIG.controlBatch$table %>%
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>%
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GLU"), ends_with("OLIG"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GLU"), ends_with("OLIG")))
write_tsv(results.qlf.GLUvsOLIG.controlBatch, "edgeR_glmQLFTest_GLUvsOLIG_DEgenes_controlBatch.txt")

# GABA vs OLIG
qlf.GABAvsOLIG.controlBatch <- glmQLFTest(fit, contrast = c(0, 0, 0, 0, 0, -1))
tt.GABAvsOLIG.controlBatch <- topTags(qlf.GABAvsOLIG.controlBatch, n=dim(y)[1])
results.qlf.GABAvsOLIG.controlBatch <- tt.GABAvsOLIG.controlBatch$table %>%
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>%
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GABA"), ends_with("OLIG"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GABA"), ends_with("OLIG")))
write_tsv(results.qlf.GABAvsOLIG.controlBatch, "edgeR_glmQLFTest_GABAvsOLIG_DEgenes_controlBatch.txt")


## control for subject 
design <- model.matrix(~sampleInfo$subject+sampleInfo$cellType)
rownames(design) <- sampleInfo$sample
design

y <- DGEList(counts = expData.sub, group = sampleInfo$cellType)
keep <- rowSums(cpm(y)>2) >= 9 
y <- y[keep, , keep.lib.size = FALSE]
y <- calcNormFactors(y, method = "TMM")

y <- estimateDisp(y, design, robust = T)
fit <- glmQLFit(y, design, robust = T, coef = 10:11)

CPM.filtered <- cpm(y) %>% as_data_frame() %>% 
  rename_(.dots = setNames(names(.), paste0("CPM_", names(.)))) %>%
  mutate(geneID = rownames(y$counts)) %>%
  select(geneID, everything())

# GLU vs GABA
qlf.GLUvsGABA.controlSubject <- glmQLFTest(fit, coef = 10)
tt.GLUvsGABA.controlSubject <- topTags(qlf.GLUvsGABA.controlSubject, n=dim(y)[1])
results.qlf.GLUvsGABA.controlSubject <- tt.GLUvsGABA.controlSubject$table %>% 
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GLU"), ends_with("GABA"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GLU"), ends_with("GABA")))
  
write_tsv(results.qlf.GLUvsGABA.controlSubject, "edgeR_glmQLFTest_GLUvsGABA_DEgenes_controlSubject.txt")

# GLU vs OLIG
qlf.GLUvsOLIG.controlSubject <- glmQLFTest(fit, contrast = c(rep(0, 9), 1, -1))
tt.GLUvsOLIG.controlSubject <- topTags(qlf.GLUvsOLIG.controlSubject, n=dim(y)[1])
results.qlf.GLUvsOLIG.controlSubject <- tt.GLUvsOLIG.controlSubject$table %>%
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>%
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GLU"), ends_with("OLIG"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GLU"), ends_with("OLIG")))
write_tsv(results.qlf.GLUvsOLIG.controlSubject, "edgeR_glmQLFTest_GLUvsOLIG_DEgenes_controlSubject.txt")

# GABA vs OLIG
qlf.GABAvsOLIG.controlSubject <- glmQLFTest(fit, contrast = c(rep(0, 10), -1))
tt.GABAvsOLIG.controlSubject <- topTags(qlf.GABAvsOLIG.controlSubject, n=dim(y)[1])
results.qlf.GABAvsOLIG.controlSubject <- tt.GABAvsOLIG.controlSubject$table %>%
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>%
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GABA"), ends_with("OLIG"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GABA"), ends_with("OLIG")))
write_tsv(results.qlf.GABAvsOLIG.controlSubject, "edgeR_glmQLFTest_GABAvsOLIG_DEgenes_controlSubject.txt")

## control for subject and batch

design <- model.matrix(~sampleInfo$subject+sampleInfo$batch+sampleInfo$cellType)
rownames(design) <- sampleInfo$sample
design

y <- DGEList(counts = expData.sub, group = sampleInfo$cellType)
keep <- rowSums(cpm(y)>2) >= 9 
y <- y[keep, , keep.lib.size = FALSE]
y <- calcNormFactors(y, method = "TMM")

y <- estimateDisp(y, design, robust = T)
fit <- glmQLFit(y, design, robust = T, coef = 13:14)

CPM.filtered <- cpm(y) %>% as_data_frame() %>% 
  rename_(.dots = setNames(names(.), paste0("CPM_", names(.)))) %>%
  mutate(geneID = rownames(y$counts)) %>%
  select(geneID, everything())

# GLU vs GABA
qlf.GLUvsGABA.controlSubjectAndBatch <- glmQLFTest(fit, coef = 13)
tt.GLUvsGABA.controlSubjectAndBatch <- topTags(qlf.GLUvsGABA.controlSubjectAndBatch, n=dim(y)[1])
results.qlf.GLUvsGABA.controlSubjectAndBatch <- tt.GLUvsGABA.controlSubjectAndBatch$table %>% 
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>% 
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GLU"), ends_with("GABA"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GLU"), ends_with("GABA")))
  
write_tsv(results.qlf.GLUvsGABA.controlSubjectAndBatch, "edgeR_glmQLFTest_GLUvsGABA_DEgenes_controlSubjectAndBatch.txt")

# GLU vs OLIG
qlf.GLUvsOLIG.controlSubjectAndBatch <- glmQLFTest(fit, contrast = c(rep(0, 12), 1, -1))
tt.GLUvsOLIG.controlSubjectAndBatch <- topTags(qlf.GLUvsOLIG.controlSubjectAndBatch, n=dim(y)[1])
results.qlf.GLUvsOLIG.controlSubjectAndBatch <- tt.GLUvsOLIG.controlSubjectAndBatch$table %>%
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>%
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GLU"), ends_with("OLIG"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GLU"), ends_with("OLIG")))
write_tsv(results.qlf.GLUvsOLIG.controlSubjectAndBatch, "edgeR_glmQLFTest_GLUvsOLIG_DEgenes_controlSubjectAndBatch.txt")

# GABA vs OLIG
qlf.GABAvsOLIG.controlSubjectAndBatch <- glmQLFTest(fit, contrast = c(rep(0, 13), -1))
tt.GABAvsOLIG.controlSubjectAndBatch <- topTags(qlf.GABAvsOLIG.controlSubjectAndBatch, n=dim(y)[1])
results.qlf.GABAvsOLIG.controlSubjectAndBatch <- tt.GABAvsOLIG.controlSubjectAndBatch$table %>%
  as_data_frame() %>%
  mutate(geneID = rownames(.)) %>%
  left_join(geneAnnotation) %>%
  select(geneID, geneName, geneType, geneStatus, everything()) %>%
  left_join(tpmData.renamed %>% select(geneID, ends_with("GABA"), ends_with("OLIG"))) %>%
  left_join(CPM.filtered %>% select(geneID, ends_with("GABA"), ends_with("OLIG")))
write_tsv(results.qlf.GABAvsOLIG.controlSubjectAndBatch, "edgeR_glmQLFTest_GABAvsOLIG_DEgenes_controlSubjectAndBatch.txt")
