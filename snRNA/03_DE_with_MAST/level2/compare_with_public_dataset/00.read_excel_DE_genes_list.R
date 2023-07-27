library(tidyverse)
library(readxl)
library(fs)

# Manolis ALS single cell DE genes list, June 2021 submission ver.
file_1 <- "./Manolis_ALS_single_cell_DE_C9_ALS_June2021_submission.xlsx"

file_1_sheets <- excel_sheets(file_1) %>% set_names()

file_1_DE <- map_dfr(
  file_1_sheets,
  ~ read_excel(path = file_1, sheet = .),
  .id = "cell_type"
) %>%
  mutate(comparison = "Manolis_C9ALS_vs_Control")

file_2 <- "./Manolis_ALS_single_cell_DE_C9_FTLD_June2021_submission.xlsx"

file_2_sheets <- excel_sheets(file_2) %>% set_names()

file_2_DE <- map_dfr(
  file_2_sheets,
  ~ read_excel(path = file_2, sheet = .),
  .id = "cell_type"
) %>%
  mutate(comparison = "Manolis_C9FTD_vs_Control")

file_3 <- "./Manolis_ALS_single_cell_DE_sporadic_ALS_June2021_submission.xlsx"

file_3_sheets <- excel_sheets(file_3) %>% set_names()

file_3_DE <- map_dfr(
  file_3_sheets,
  ~ read_excel(path = file_3, sheet = .),
  .id = "cell_type"
) %>%
  mutate(comparison = "Manolis_sALS_vs_Control")

file_4 <- "./Manolis_ALS_single_cell_DE_sporadic_FTLD_June2021_submission.xlsx"

file_4_sheets <- excel_sheets(file_4) %>% set_names()

file_4_DE <- map_dfr(
  file_4_sheets,
  ~ read_excel(path = file_4, sheet = .),
  .id = "cell_type"
) %>%
  mutate(comparison = "Manolis_sFTD_vs_Control")

res <- bind_rows(file_1_DE, file_2_DE, file_3_DE, file_4_DE)
write_tsv(res, "Manolis_ALS_single_cell_DE_genes_June2021_submission.tsv")

# Mathys, Nature 2019, AD, prefrontal cortex, no pathlogy vs. pathology
file_5 <- "./Mathys_Nat2019_AD_prefrontalCortex_single_cell_DE_noPathology_vs_Pathology.xlsx"
file_5_sheets <- excel_sheets(file_5) %>% set_names()

file_5_DE <- map_dfr(
  file_5_sheets,
  ~ read_excel(
    path = file_5,
    sheet = .,
    col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "text", "text")
  ),
  .id = "cell_type"
) %>%
  mutate(comparison = "Mathys_Nat2019_AD_prefrontalCortex_noPathology_vs_Pathology")
write_tsv(file_5_DE, "Mathys_Nat2019_AD_prefrontalCortex_single_cell_DE_noPathology_vs_Pathology.tsv")

# Morabi, Nature Genetics 2021, AD, prefrontal cortex, diagnosis DE by cell type and cluster
file_6 <- "./Morabito_NatGenet2021_AD_prefrontalCortex_single_cell_cell_type_and_cluster_diagnosis_DE.xlsx"
file_6_sheets <- excel_sheets(file_6) %>% set_names()

file_6_DE_cell_type <- read_excel(
  path = file_6,
  sheet = file_6_sheets[1]
) %>%
  mutate(comparison = "Morabito_NatGenet2021_AD_prefrontalCortex_cell_type_diagnosis_DE")

write_tsv(file_6_DE_cell_type, "Morabito_NatGenet2021_AD_prefrontalCortex_single_cell_cell_type_diagnosis_DE.tsv")

file_6_DE_cluster <- read_excel(
  path = file_6,
  sheet = file_6_sheets[2]
) %>%
  mutate(comparison = "Morabito_NatGenet2021_AD_prefrontalCortex_cluster_diagnosis_DE")

write_tsv(file_6_DE_cluster, "Morabito_NatGenet2021_AD_prefrontalCortex_single_cell_cluster_diagnosis_DE.tsv")

# Grubman, Nature 2019, AD, EntorhinalCortex
read_Grubman <- function(cell_type) {
  read_csv(
    str_glue("./Grubman_AD_single_cell_DE_{cell_type}.csv"),
    col_types = "cdd"
  ) %>%
    mutate(
      cell_type = cell_type
    )
}

Grubman_DE <- map_dfr(
  c("astrocyte", "endo", "microglia", "neuron", "oligo", "OPC"),
  read_Grubman
)

write_tsv(Grubman_DE, "Grubman_Nat2019_AD_EntorhinalCortex_single_Cell_DE.tsv")
