library(tidyverse)
library(edgeR)

count_mat <- read_tsv("peak_Astro_ALS_vs_Control_atac_count_by_major_clusters_by_region_by_disease_by_individual.txt")
pool_libsize <- read_tsv("peak_Astro_ALS_vs_Control_pooled_libSize_by_major_clusters_by_region_by_disease_by_individual.txt")

run_edger <- function(selected_group1, selected_group2){
  
  message(paste0("Working on ", selected_group1, " vs ", selected_group2))
  
  df_test <- count_mat %>% 
    select(peak, contains(selected_group1), contains(selected_group2)) %>% 
    as.data.frame() %>% 
    column_to_rownames("peak")
  
  meta_test <- tibble(sample = colnames(df_test)) %>% 
    mutate(Subject = str_extract(sample, "subject.*"),
           Subject = str_replace(Subject, "subject", ""),
           sample2 = sample) %>% 
    separate(sample2, c("cell_type", "brain_region", "disease", "subject")) %>% 
    left_join(pool_libsize, by = c("sample" = "group"))
  
  y <- DGEList(counts = df_test,
               group = meta_test$disease,
               lib.size = meta_test$pooled_lib_size,
               remove.zeros = TRUE) 
  
  keep <- rowSums(cpm(y)>1) >= 6
  # table(keep)
  y <- y[keep, , keep.lib.size = TRUE]
  # y <- y[keep, , keep.lib.size = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  # y$samples
  
  # design <- model.matrix(~ factor(meta_test$Sex) + relevel(factor(meta_test$disease), ref = "Control"))
  design <- model.matrix(~ relevel(factor(meta_test$disease), ref = "Control"))
  
  y <- estimateDisp(y, design, robust = T)
  fit <- glmQLFit(y, design, robust = T)
  message(y$common.dispersion)
  # plotBCV(y)
  # plotQLDisp(fit)
  
  cpm_filtered <- cpm(y) %>% 
    as_tibble() %>%
    rename_all(funs(paste0("CPM_", .))) %>% 
    mutate(gene = rownames(y$counts)) %>% 
    select(gene, everything())
  
  # qlf <- glmQLFTest(fit, coef = 3)
  qlf <- glmQLFTest(fit, coef = 2)
  tt <- topTags(qlf, n=dim(y)[1])
  res <- tt$table %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    left_join(cpm_filtered)
  
  message(res %>% filter(FDR < 0.05) %>% nrow())
  
  write_tsv(res, paste0("edgeR_atac_peaks_", selected_group1, "_vs_", selected_group2, ".txt"))
}

# run_edger("Astro_MCX_ALS", "_MCX_control")
# 
# param <- tibble(selected_group1 = colnames(count_mat_sub)) %>% 
#   filter(selected_group1 != "peak") %>% 
#   filter(!grepl("Mixed", selected_group1)) %>% 
#   filter(!grepl("_control", selected_group1)) %>% 
#   mutate(selected_group1 = str_replace(selected_group1, "_subject.*", "")) %>% 
#   distinct(selected_group1) %>% 
#   mutate(selected_group2 = str_replace(selected_group1, "(_ALS)|(_FTD)", "_control"))
# 
# pwalk(param, run_edger)
