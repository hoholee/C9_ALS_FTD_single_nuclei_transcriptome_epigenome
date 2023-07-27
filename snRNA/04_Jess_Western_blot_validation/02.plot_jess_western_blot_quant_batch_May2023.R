# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)

# read in data
df <- read_tsv("./MCX_ALS_vs_Control_JessWesternBlot_protein_quantification_batch_May2023.txt") %>%
  mutate(diagnosis = if_else(grepl("^A", Sample), "ALS", "Control"))

p <- df %>% ggplot(aes(diagnosis, area))
p +
  geom_boxplot(aes(fill = diagnosis)) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("ALS", "Control"))
  ) +
  facet_wrap(~ gene + experiment, scale = "free_y") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
ggsave("batch_May2023_all_boxplot_t_test_stats.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE)
dev.off()

# t-test
get_t_test <- function(test_gene, selected_experiment) {
  to_test <- df %>%
    filter(gene == test_gene, experiment == selected_experiment)
  res <- t.test(area ~ diagnosis, data = to_test)
  tibble(
    gene = test_gene,
    experiment = selected_experiment,
    t = res$statistic,
    df = res$parameter,
    p = res$p.value,
    conf_level = attributes(res$conf.int)$conf.level,
    conf_low = res$conf.int[1],
    conf_hi = res$conf.int[2],
    mean_ALS = res$estimate[1],
    mean_Control = res$estimate[2],
    stderr = res$stderr,
    alternative = res$alternative,
    method = res$method
  )
}

gene_list <- df %>%
  distinct(gene) %>%
  pull(gene)
params <- expand_grid(
  test_gene = gene_list,
  selected_experiment = c("1", "2")
)

t_test_res <- pmap_dfr(params, get_t_test) %>%
  mutate(fdr = p.adjust(p, method = "fdr"))
write_tsv(t_test_res, "./MCX_ALS_vs_Control_JessWesternBlot_protein_quantification_tTest_batch_May2023.txt")

# plot separately for Astro DE genes and Neuronal DE genes
disease_color_panel <- c("ALS" = "#F6416C", "Control" = "#00B8A9")

# Astro
astro_genes <- c("MAOB")
df_astro <- df %>%
  filter(gene %in% astro_genes) %>%
  mutate(gene = factor(gene, levels = astro_genes))

p <- df_astro %>% ggplot(aes(diagnosis, area / 10000))
p +
  geom_boxplot(
    aes(color = diagnosis),
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_beeswarm(
    aes(color = diagnosis),
    size = 1,
    shape = 21,
    priority = "density",
    cex = 6
  ) +
  stat_compare_means(
    aes(label = ..p.signif..),
    method = "t.test",
    comparisons = list(c("ALS", "Control"))
  ) +
  scale_color_manual(values = disease_color_panel) +
  facet_wrap(~ paste(gene, experiment, sep = "_"), scale = "free_y", nrow = 1) +
  ylab("Western blot protein level (A.U.)") +
  xlab("Diagnosis") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(
  "western_protein_quantification_selected_Astro_genes_batch_May2023.pdf",
  device = cairo_pdf(),
  width = 1.6,
  height = 1.5,
  useDingbats = FALSE
)
dev.off()

# Neurons
neu_genes <- c(
  "ATP5A", "CYCS", "COX5A"
)

df_neu <- df %>%
  filter(gene %in% neu_genes) %>%
  mutate(
    gene = factor(gene, levels = neu_genes),
    group = paste(gene, experiment, sep = "_"),
    group = factor(group, levels = c("ATP5A_1", "ATP5A_2", "CYCS_1", "CYCS_2", "COX5A_1", "COX5A_2"))
  )

p <- df_neu %>% ggplot(aes(diagnosis, area / 10000))
p +
  geom_boxplot(
    aes(color = diagnosis),
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_beeswarm(
    aes(color = diagnosis),
    size = 1,
    shape = 21,
    priority = "density",
    cex = 6
  ) +
  stat_compare_means(
    aes(label = ..p.signif..),
    method = "t.test",
    comparisons = list(c("ALS", "Control"))
  ) +
  scale_color_manual(values = disease_color_panel) +
  # facet_wrap(~paste(gene, experiment, sep = "_"), scale = "free_y", nrow = 1) +
  facet_wrap(~group, scale = "free_y", nrow = 1) +
  ylab("Western blot protein level (A.U.)") +
  xlab("Diagnosis") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(
  "western_protein_quantification_selected_Neu_genes_batch_May2023.pdf",
  device = cairo_pdf(),
  width = 6.4,
  height = 1.5,
  useDingbats = FALSE
)
dev.off()

## log sessionInfo
sessionInfo()
