library(tidyverse)
library(ggrepel)
library(ggrastr)

# df_obs <- read_tsv("./MAST_res_level2_summary.txt")
df_obs <- read_tsv("./MAST_res_level2_rmChrXY_summary.txt")

selected_cond_1 <- "ALS"
selected_cond_2 <- "Control"
selected_cell_type <- "Astro"
# selected_cell_type <- "Micro"
selected_region <- "MCX"
# selected_region <- "mFCX"

# df_shuffle <- read_tsv(
#   str_glue(
#     "./results_shuffle/MAST_res_level2_",
#     "{selected_cell_type}_{selected_region}_",
#     "{selected_cond_1}_vs_{selected_cond_2}",
#     ".txt"
#   ),
#   col_types = "cdddddlldddicc"
# )

df_shuffle <- read_tsv(
  str_glue(
    "./results_shuffle_rmChrXY/MAST_res_level2_",
    "{selected_cell_type}_{selected_region}_",
    "{selected_cond_1}_vs_{selected_cond_2}",
    ".txt"
  ),
  col_types = "cdddddlldddicc"
)

FDR_threshold <- 0.05
FC_threshold <- 1.2

# df_obs_sub <- df_obs %>%
#   filter(
#     cond_1 == selected_cond_1,
#     cond_2 == selected_cond_2,
#     cell_type == selected_cell_type,
#     region == selected_region,
#     FDR < FDR_threshold,
#     abs(model_log2FC) > log2(FC_threshold),
#     conv_C == TRUE,
#     conv_D == TRUE,
#     model_log2FC_ci_hi * model_log2FC_ci_low > 0,
#     abs(model_log2FC - avg_logFC) < 2
#   )

df_obs_sub <- df_obs %>%
  filter(
    cond_1 == selected_cond_1,
    cond_2 == selected_cond_2,
    cell_type == selected_cell_type,
    region == selected_region,
    fdr < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    ci.hi * ci.lo > 0,
    abs(model_log2FC - avg_log2FC) < 2
  )


df_shuffle_sub <- df_shuffle %>%
  filter(
    fdr < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    ci.hi * ci.lo > 0,
    abs(model_log2FC - avg_log2FC) < 2
  )

df_shuffle_sig_count <- df_shuffle_sub %>%
  dplyr::count(seed, permuted_subject_disease_1, permuted_subject_disease_2)

p <- df_shuffle_sig_count %>% ggplot(aes(n))
p +
  geom_histogram() +
  geom_vline(xintercept = nrow(df_obs_sub), linetype = 2, color = "red") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank()
  )


# add sample metadata
sample_metadata <- read_tsv(
  "../../meta_preQC_all_cells_post_cellBender.txt",
  col_types = "cciicccciicdcdl"
) %>%
  distinct(disease, subject, sex, age, PMI)

shuffle_DE_count_by_sex <- df_shuffle_sig_count %>%
  pivot_longer(
    cols = c("permuted_subject_disease_1", "permuted_subject_disease_2"),
    names_to = "group",
    values_to = "subjects"
  ) %>%
  separate(
    subjects,
    into = c("subject_1", "subject_2", "subject_3", "subject_4", "subject_5", "subject_6"),
    sep = ","
  ) %>%
  pivot_longer(
    cols = c("subject_1", "subject_2", "subject_3", "subject_4", "subject_5", "subject_6"),
    names_to = "subject_idx",
    values_to = "subject"
  ) %>%
  left_join(
    sample_metadata,
    by = "subject"
  ) %>%
  dplyr::count(seed, n, group, sex) %>%
  pivot_wider(
    names_from = sex,
    values_from = nn,
    values_fill = 0
  ) %>%
  filter(group == "permuted_subject_disease_1")

p <- shuffle_DE_count_by_sex %>% ggplot(aes(F, n))
p + geom_boxplot(aes(group = F)) +
  geom_hline(yintercept = nrow(df_obs_sub), linetype = 2, color = "red") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank()
  )



p <- df_obs %>%
  filter(
    cond_1 == "ALS",
    cond_2 == "Control",
    cell_type == "Astro",
    region == "MCX"
  ) %>%
  ggplot(aes(model_log2FC, -log10(FDR)))

p + geom_point_rast() +
  geom_point_rast(data = df_obs_sub, color = "red")

p <- df_shuffle %>%
  filter(seed == 18) %>%
  ggplot(aes(model_log2FC, -log10(fdr)))
# p + geom_point_rast()
p + geom_point_rast() +
  geom_point_rast(data = df_shuffle_sub, color = "red")


# check if observed significant DE genes have higher average log2FC than random shuffles
balacned_seed_list <- shuffle_DE_count_by_sex %>%
  filter(F == 3) %>%
  pull(seed)

shuffle_quantile <- df_shuffle %>%
  filter(seed %in% balacned_seed_list) %>%
  group_by(gene) %>%
  summarise(quibble(model_log2FC, avg_log2FC, c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99))) %>%
  ungroup()

################################################################################# 3

quibble <- function(x, y, q = c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)) {
  tibble(
    "{{ x }}_quantile" := quantile(x, q, na.rm = TRUE),
    "{{ y }}_quantile" := quantile(y, q, na.rm = TRUE),
    "q_tile" := q
  )
}


shuffle_quantile <- df_shuffle %>%
  group_by(gene) %>%
  summarise(quibble(model_log2FC, avg_log2FC, c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99))) %>%
  ungroup()

shuffle_quantile_top <- shuffle_quantile %>%
  # filter(q_tile == 0.99) %>%
  filter(q_tile == 0.95) %>%
  select(-q_tile) %>%
  rename(
    model_log2FC_quantile_top = model_log2FC_quantile,
    avg_log2FC_quantile_top = avg_log2FC_quantile
  )
shuffle_quantile_bottom <- shuffle_quantile %>%
  # filter(q_tile == 0.01) %>%
  filter(q_tile == 0.05) %>%
  select(-q_tile) %>%
  rename(
    model_log2FC_quantile_bottom = model_log2FC_quantile,
    avg_log2FC_quantile_bottom = avg_log2FC_quantile
  )

# df_obs_mod <- df_obs %>%
#   left_join(shuffle_quantile_top) %>%
#   left_join(shuffle_quantile_bottom) %>%
#   mutate(
#     shuffle_99quantile_log2FC = if_else(model_log2FC > 0, model_log2FC_quantile_top, model_log2FC_quantile_bottom),
#     over_shuffle = if_else(abs(model_log2FC) > abs(shuffle_99quantile_log2FC), "yes", "no")
#   )
# df_to_plot <- df_obs_mod %>%
#   filter(
#     fdr < 0.05,
#     conv_C == TRUE,
#     conv_D == TRUE,
#     ci.hi * ci.lo > 0,
#     abs(model_log2FC) > log2(1.1)
#   )

# p <- df_to_plot %>%
#   ggplot(aes(shuffle_99quantile_log2FC, model_log2FC))

df_obs_mod <- df_obs_sub %>%
  left_join(shuffle_quantile_top) %>%
  left_join(shuffle_quantile_bottom) %>%
  mutate(
    shuffle_99quantile_log2FC = if_else(model_log2FC > 0, model_log2FC_quantile_top, model_log2FC_quantile_bottom),
    over_shuffle = if_else(abs(model_log2FC) > abs(shuffle_99quantile_log2FC), "yes", "no")
  )

p <- df_obs_mod %>%
  ggplot(aes(shuffle_99quantile_log2FC, model_log2FC))

p + geom_point(aes(color = over_shuffle)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_text_repel(aes(label = gene), data = df_obs_mod %>% filter(over_shuffle == "yes")) +
  xlab("95% quantile of log2FC in shuffle") +
  ylab("Observed log2FC") +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  scale_color_manual(values = c("grey70", "#DE325C"))

# ggsave("test_permutation_res_allrandom.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)
# ggsave("test_permutation_res_halfhalf.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = F)

# distribution of the log2FC in shuffle for individual gene
selected_gene <- "CHI3L1"
selected_gene <- "COLEC12"
selected_gene <- "RANBP3L"
selected_gene <- "MAP2"
selected_gene <- "PCLO"
selected_gene_obs_log2FC <- df_obs_mod %>%
  filter(gene == selected_gene) %>%
  pull(model_log2FC)

selected_gene_shuffle_log2FC <- df_shuffle %>%
  filter(seed %in% balacned_seed_list) %>%
  filter(gene == selected_gene) %>%
  select(model_log2FC, seed)

p <- selected_gene_shuffle_log2FC %>%
  ggplot(aes(model_log2FC))
p +
  geom_histogram() +
  # geom_density() +
  geom_vline(xintercept = selected_gene_obs_log2FC, color = "red") +
  theme_bw(base_size = 10, base_family = "Helvetica")
