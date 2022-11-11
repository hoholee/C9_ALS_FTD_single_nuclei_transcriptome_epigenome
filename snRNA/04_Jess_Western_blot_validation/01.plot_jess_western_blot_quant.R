# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)

# read in data
df <- read_tsv("./MCX_ALS_vs_Control_JessWesternBlot_protein_quantification.txt")

df_mod <- df %>%
  select(Sample, contains("Area")) %>%
  pivot_longer(-Sample, names_to = "gene", values_to = "area") %>%
  mutate(gene = str_replace(gene, "-Area", "")) %>%
  group_by(Sample, gene) %>%
  summarise(area_avg = mean(area)) %>%
  ungroup() %>%
  mutate(diagnosis = if_else(grepl("C9", Sample), "ALS", "Control"))

p <- df_mod %>% ggplot(aes(diagnosis, area_avg))
p +
  geom_boxplot(aes(fill = diagnosis)) +
  stat_compare_means(
    aes(label = ..p.signif..),
    method = "t.test",
    comparisons = list(c("ALS", "Control"))
  ) +
  facet_wrap(~gene, scale = "free_y") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
ggsave("tmp.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE)


# subset
selected_genes <- c(
  "GFAP", "CD44", "CHI3L1", "RANBP3L", "TGFB2",
  "HSP90", "CLU- Precursor Area", "HSP27", "RAP1GAP", "KCND3/Kv4.3",
  "DNMT3A", "UBB", "HSP70", "NEFL", "SOD1"
)

selected_genes_clean <- selected_genes %>%
  str_replace("- Precursor Area", "") %>%
  str_replace("/Kv4.3", "")

df_mod_sub <- df_mod %>%
  filter(gene %in% selected_genes) %>%
  mutate(
    gene = str_replace(gene, "- Precursor Area", ""),
    gene = str_replace(gene, "/Kv4.3", ""),
    gene = factor(gene, levels = selected_genes_clean)
  )

write_tsv(df_mod_sub, "./MCX_ALS_vs_Control_JessWesternBlot_protein_quantification_selected_genes_signalAvgReps.txt")

# t-test
get_t_test <- function(test_gene) {
  to_test <- df_mod_sub %>%
    filter(gene == test_gene)
  res <- t.test(area_avg ~ diagnosis, data = to_test)
  tibble(
    gene = test_gene,
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

t_test_res <- map_dfr(selected_genes_clean, get_t_test) %>%
  mutate(fdr = p.adjust(p, method = "fdr"))
write_tsv(t_test_res, "./MCX_ALS_vs_Control_JessWesternBlot_protein_quantification_selected_genes_signalAvgReps_tTest.txt")

# plot separately for Astro DE genes and Neuronal DE genes
disease_color_panel <- c("ALS" = "#F6416C", "Control" = "#00B8A9")

# Astro
astro_genes <- c("GFAP", "CD44", "CHI3L1", "RANBP3L", "TGFB2")
df_astro <- df_mod_sub %>%
  filter(gene %in% astro_genes) %>%
  mutate(gene = factor(gene, levels = astro_genes))

p <- df_astro %>% ggplot(aes(diagnosis, area_avg / 10000))
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
  facet_wrap(~gene, scale = "free_y", nrow = 1) +
  ylab("Western blot protein level (A.U.)") +
  xlab("Diagnosis") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(
  "western_protein_quantification_selected_Astro_genes.pdf",
  device = cairo_pdf(),
  width = 4,
  height = 1.5,
  useDingbats = FALSE
)

# Neurons
neu_genes <- c(
  "HSP90", "CLU", "HSP27", "UBB", "HSP70",
  "SOD1", "NEFL", "RAP1GAP", "DNMT3A", "KCND3"
)

df_neu <- df_mod_sub %>%
  filter(gene %in% neu_genes) %>%
  mutate(gene = factor(gene, levels = neu_genes))

p <- df_neu %>% ggplot(aes(diagnosis, area_avg / 10000))
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
  facet_wrap(~gene, scale = "free_y", nrow = 1) +
  ylab("Western blot protein level (A.U.)") +
  xlab("Diagnosis") +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(
  "western_protein_quantification_selected_Neu_genes.pdf",
  device = cairo_pdf(),
  width = 8,
  height = 1.5,
  useDingbats = FALSE
)

## log sessionInfo
sessionInfo()