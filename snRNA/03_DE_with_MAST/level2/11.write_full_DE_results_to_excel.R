library(tidyverse)
library(readxl)
library(writexl)

df <- read_tsv("./MAST_res_level2_summary.txt")

df_ALS_motor <- df %>%
    filter(cond_1 == "ALS", cond_2 == "Control", region == "MCX", cell_type != "Exc_unknown") %>%
    select(-cond_1, -cond_2, -region, -num_cells_cond_1, -num_cells_cond_2) %>%
    mutate(
        cell_type = case_when(
            cell_type == "Exc_superficial" ~ "Exc_upper",
            cell_type == "Inh_ADARB2_Other" ~ "Inh_other_CGE",
            TRUE ~ cell_type
        )
    ) %>%
    rename(
        avg_log2CPM_C9ALS = avg_log2CPM_cond_1,
        avg_log2CPM_Control = avg_log2CPM_cond_2,
        avg_log2FC = avg_logFC
    ) %>%
    group_by(cell_type)

df_ALS_frontal <- df %>%
    filter(cond_1 == "ALS", cond_2 == "Control", region == "mFCX", cell_type != "Exc_unknown") %>%
    select(-cond_1, -cond_2, -region, -num_cells_cond_1, -num_cells_cond_2) %>%
    mutate(
        cell_type = case_when(
            cell_type == "Exc_superficial" ~ "Exc_upper",
            cell_type == "Inh_ADARB2_Other" ~ "Inh_other_CGE",
            TRUE ~ cell_type
        )
    ) %>%
    rename(
        avg_log2CPM_C9ALS = avg_log2CPM_cond_1,
        avg_log2CPM_Control = avg_log2CPM_cond_2,
        avg_log2FC = avg_logFC
    ) %>%
    group_by(cell_type)

df_FTD_motor <- df %>%
    filter(cond_1 == "FTD", cond_2 == "Control", region == "MCX", cell_type != "Exc_unknown") %>%
    select(-cond_1, -cond_2, -region, -num_cells_cond_1, -num_cells_cond_2) %>%
    mutate(
        cell_type = case_when(
            cell_type == "Exc_superficial" ~ "Exc_upper",
            cell_type == "Inh_ADARB2_Other" ~ "Inh_other_CGE",
            TRUE ~ cell_type
        )
    ) %>%
    rename(
        avg_log2CPM_C9FTD = avg_log2CPM_cond_1,
        avg_log2CPM_Control = avg_log2CPM_cond_2,
        avg_log2FC = avg_logFC
    ) %>%
    group_by(cell_type)

df_FTD_frontal <- df %>%
    filter(cond_1 == "FTD", cond_2 == "Control", region == "mFCX", cell_type != "Exc_unknown") %>%
    select(-cond_1, -cond_2, -region, -num_cells_cond_1, -num_cells_cond_2) %>%
    mutate(
        cell_type = case_when(
            cell_type == "Exc_superficial" ~ "Exc_upper",
            cell_type == "Inh_ADARB2_Other" ~ "Inh_other_CGE",
            TRUE ~ cell_type
        )
    ) %>%
    rename(
        avg_log2CPM_C9FTD = avg_log2CPM_cond_1,
        avg_log2CPM_Control = avg_log2CPM_cond_2,
        avg_log2FC = avg_logFC
    ) %>%
    group_by(cell_type)


list_df_ALS_motor <- df_ALS_motor %>% group_split(.keep = FALSE)
list_df_ALS_frontal <- df_ALS_frontal %>% group_split(.keep = FALSE)
list_df_FTD_motor <- df_FTD_motor %>% group_split(.keep = FALSE)
list_df_FTD_frontal <- df_FTD_frontal %>% group_split(.keep = FALSE)

names(list_df_ALS_motor) <- df_ALS_motor %>%
    group_keys() %>%
    pull(cell_type)
names(list_df_ALS_frontal) <- df_ALS_frontal %>%
    group_keys() %>%
    pull(cell_type)
names(list_df_FTD_motor) <- df_FTD_motor %>%
    group_keys() %>%
    pull(cell_type)
names(list_df_FTD_frontal) <- df_FTD_frontal %>%
    group_keys() %>%
    pull(cell_type)


write_xlsx(list_df_ALS_motor, "MAST_DE_major_cell_types_ALS_vs_Control_motor.xlsx")
write_xlsx(list_df_ALS_frontal, "MAST_DE_major_cell_types_ALS_vs_Control_frontal.xlsx")
write_xlsx(list_df_FTD_motor, "MAST_DE_major_cell_types_FTD_vs_Control_motor.xlsx")
write_xlsx(list_df_FTD_frontal, "MAST_DE_major_cell_types_FTD_vs_Control_frontal.xlsx")
