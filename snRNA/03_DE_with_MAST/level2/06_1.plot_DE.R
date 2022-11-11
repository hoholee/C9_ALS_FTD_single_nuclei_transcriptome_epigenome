library(tidyverse)

DE_level_2 <- read_tsv("MAST_res_level2_summary.txt")

FDR_threshold <- 0.05
FC_threshold <- 1.2

selected_cell_types <- c(
  "Exc_superficial", "Exc_intermediate", "Exc_deep",
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
  "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

DE_level_2_sig <- DE_level_2 %>%
  filter(
    FDR < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    abs(model_log2FC - avg_logFC) < 2
  ) %>%
  group_by(cond_1, cond_2, cell_type, region) %>%
  count() %>%
  ungroup() %>%
  complete(
    fill = list(n = 0),
    cond_1 = unique(cond_1),
    cond_2 = unique(cond_2),
    cell_type = unique(cell_type),
    region = unique(region)
  ) %>%
  filter(cond_2 == "Control", cell_type != "Exc_unknown") %>%
  mutate(cell_type = factor(cell_type, levels = rev(selected_cell_types)))

region_labels <- c("MCX" = "motor cortex", "mFCX" = "frontal cortex")
region_color <- c("MCX" = "#432266", "mFCX" = "#FAA51B")
relabel_region <- function(x) {
  region_labels[x]
}

p <- DE_level_2_sig %>% ggplot(aes(cell_type, n))
p + geom_bar(stat = "identity", aes(fill = region), position = position_dodge(0.8), width = 0.8) +
  facet_wrap(~cond_1, ncol = 2) +
  scale_fill_manual(labels = relabel_region(), values = region_color) +
  ylab(paste0("# significant DE genes (FDR < ", FDR_threshold, ", FC > ", FC_threshold, ")")) +
  xlab("Major cell types") +
  coord_flip() +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
ggsave("MAST_num_sig_DE_fdr0.05_FC1.2_level_2_ordered.pdf",
  device = cairo_pdf(), width = 4, height = 2, useDingbats = F
)

# separate ALS vs. Control and FTD vs. Control
# separate up- and down-regulated DE genes
DE_level_2_sig_up_down <- DE_level_2 %>%
  filter(
    FDR < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    abs(model_log2FC - avg_logFC) < 2
  ) %>%
  mutate(DE_group = if_else(model_log2FC > 0, "up", "down")) %>%
  group_by(cond_1, cond_2, cell_type, region, DE_group) %>%
  count() %>%
  ungroup() %>%
  complete(
    fill = list(n = 0),
    cond_1 = unique(cond_1),
    cond_2 = unique(cond_2),
    cell_type = unique(cell_type),
    region = unique(region),
    DE_group = unique(DE_group),
  ) %>%
  filter(cond_2 == "Control", cell_type != "Exc_unknown") %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(selected_cell_types)),
    DE_group = factor(DE_group, levels = c("up", "down"))
  )

DE_group_labels <- c("up" = "up-regulated", "down" = "down-regulated")
DE_group_color <- c("up" = "#DE315D", "down" = "#14668E")
relabel_DE_group <- function(x) {
  DE_group_labels[x]
}

p <- DE_level_2_sig_up_down %>%
  filter(cond_1 == "ALS") %>%
  ggplot(aes(cell_type, n))
p + geom_bar(
  stat = "identity",
  aes(fill = DE_group),
  position = position_dodge(0.8),
  width = 0.8
) +
  facet_wrap(~ cond_1 + region, ncol = 2) +
  scale_fill_manual(labels = relabel_DE_group(), values = DE_group_color, name = "DE genes") +
  ylab(paste0("# significant DE genes (FDR < ", FDR_threshold, ", FC > ", FC_threshold, ")")) +
  xlab("Major cell types") +
  coord_flip() +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
ggsave("MAST_ALS_vs_Control_num_sig_up_down_DE_fdr0.05_FC1.2_level_2_ordered.pdf",
  device = cairo_pdf(), width = 4, height = 2, useDingbats = F
)

# FTD
p <- DE_level_2_sig_up_down %>%
  filter(cond_1 == "FTD") %>%
  ggplot(aes(cell_type, n))
p + geom_bar(
  stat = "identity",
  aes(fill = DE_group),
  position = position_dodge(0.8),
  width = 0.8
) +
  facet_wrap(~ cond_1 + region, ncol = 2) +
  scale_fill_manual(labels = relabel_DE_group(), values = DE_group_color, name = "DE genes") +
  ylab(paste0("# significant DE genes (FDR < ", FDR_threshold, ", FC > ", FC_threshold, ")")) +
  xlab("Major cell types") +
  coord_flip() +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
ggsave("MAST_FTD_vs_Control_num_sig_up_down_DE_fdr0.05_FC1.2_level_2_ordered.pdf",
  device = cairo_pdf(), width = 4, height = 2, useDingbats = F
)

# FTD, remove all neuron types
p <- DE_level_2_sig_up_down %>%
  filter(cond_1 == "FTD") %>%
  filter(!grepl("Exc_|Inh_", cell_type)) %>%
  ggplot(aes(cell_type, n))
p + geom_bar(
  stat = "identity",
  aes(fill = DE_group),
  position = position_dodge(0.8),
  width = 0.8
) +
  facet_wrap(~ cond_1 + region, ncol = 2) +
  scale_fill_manual(labels = relabel_DE_group(), values = DE_group_color, name = "DE genes") +
  ylab(paste0("# significant DE genes (FDR < ", FDR_threshold, ", FC > ", FC_threshold, ")")) +
  xlab("Major cell types") +
  coord_flip() +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
ggsave("MAST_FTD_vs_Control_num_sig_up_down_DE_fdr0.05_FC1.2_level_2_ordered_rmNeu.pdf",
  device = cairo_pdf(), width = 4, height = 2, useDingbats = F
)


# level 3
DE_level_3 <- read_tsv("../second_round_level_3/MAST_res_level3_summary.txt")

FDR_threshold <- 0.05
FC_threshold <- 1.2

DE_level_3_sig <- DE_level_3 %>%
  filter(
    FDR < FDR_threshold,
    abs(model_log2FC) > log2(FC_threshold),
    conv_C == TRUE,
    conv_D == TRUE,
    model_log2FC_ci_hi * model_log2FC_ci_low > 0,
    abs(model_log2FC - avg_logFC) < 2
  ) %>%
  group_by(cond_1, cond_2, cell_type, region) %>%
  count() %>%
  ungroup()

p <- DE_level_3_sig %>% ggplot(aes(cell_type, n))
p + geom_bar(stat = "identity") +
  facet_wrap(~ cond_1 + cond_2 + region, ncol = 2) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("MAST_num_sig_DE_fdr0.05_FC1.2_level_3.pdf", device = cairo_pdf(), width = 8, height = 8, useDingbats = F)