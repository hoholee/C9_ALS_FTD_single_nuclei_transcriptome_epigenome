# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(viridis)
library(scico)
library(shades)
library(ggbeeswarm)
library(ggpubr)

# read metadata
meta_data <- read_tsv(
    "metadata_all_cells_2nd_round_annotations.txt",
    col_types = "cciiccccddcddliicccccccc"
) %>%
    filter(
        rna_anno_2ndRound_level_3 != "Ambiguous",
        !rna_anno_2ndRound_level_2 %in% c("NK_cell", "Exc_unknown")
    )

color_palette_level_2 <- read_tsv("./color_palette_level_2.txt")
colors_level_2 <- color_palette_level_2$color
names(colors_level_2) <- color_palette_level_2$sub_cluster

cell_count_level2 <- meta_data %>%
    count(orig.ident, region, disease, subject, rna_anno_2ndRound_level_2) %>%
    mutate(rna_anno_2ndRound_level_2 = factor(rna_anno_2ndRound_level_2,
        levels = color_palette_level_2$sub_cluster
    )) %>%
    mutate(
        disease = factor(disease, levels = c("ALS", "FTD", "Control")),
        region = factor(region, levels = c("MCX", "mFCX"))
    )

# cell_count_level3 <- meta_data %>%
#     count(orig.ident, region, disease, subject, rna_anno_2ndRound_level_3)

#
p <- cell_count_level2 %>% ggplot(aes(subject, n))

p + geom_bar(aes(fill = rna_anno_2ndRound_level_2),
    # color = "black",
    position = position_fill(),
    stat = "identity"
) +
    facet_wrap(~ disease + region,
        dir = "v",
        drop = TRUE,
        scales = "free_x",
        nrow = 2
    ) +
    scale_fill_manual(values = colors_level_2, name = "Major cell types") +
    xlab("Individual IDs") +
    ylab("Proportion of types (%)") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        strip.background = element_blank(),
        panel.grid.minor = element_blank()
    )
ggsave("./plots/cell_type_count_proportion_level_2.pdf",
    device = cairo_pdf(),
    width = 6, height = 4, useDingbats = FALSE
)


# group by cell class, separated by individuals
cell_count_class <- meta_data %>%
    count(orig.ident, region, disease, subject, rna_anno_2ndRound_level_1) %>%
    mutate(rna_anno_2ndRound_level_1 = factor(rna_anno_2ndRound_level_1,
        levels = c("Exc_neuron", "Inh_neuron", "Non_neuron")
    ))

cell_class_color <- c(
    "Exc_neuron" = "#96BB45",
    "Inh_neuron" = "#718DC7",
    "Non_neuron" = "#9C482B"
)

#
p <- cell_count_class %>% ggplot(aes(subject, n))

p + geom_bar(aes(fill = rna_anno_2ndRound_level_1),
    position = position_fill(),
    stat = "identity"
) +
    facet_wrap(~ disease + region,
        dir = "v",
        drop = TRUE,
        scales = "free_x",
        nrow = 2
    ) +
    scale_fill_manual(values = cell_class_color, name = "Cell class") +
    xlab("Individual IDs") +
    ylab("Proportion of types (%)") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        strip.background = element_blank(),
        panel.grid.minor = element_blank()
    )
ggsave("./plots/cell_type_count_proportion_level_1.pdf",
    device = cairo_pdf(),
    width = 6, height = 4, useDingbats = FALSE
)


# group by cell class, merged all individuals, show both level 1 and level 2 distributions
cell_count_class_merged <- meta_data %>%
    count(region, disease, rna_anno_2ndRound_level_1) %>%
    mutate(rna_anno_2ndRound_level_1 = factor(rna_anno_2ndRound_level_1,
        levels = c("Exc_neuron", "Inh_neuron", "Non_neuron")
    )) %>%
    rename(group = rna_anno_2ndRound_level_1) %>%
    mutate(disease = paste0(disease, "_level_1"))

cell_count_type_merged <- meta_data %>%
    count(region, disease, rna_anno_2ndRound_level_2) %>%
    mutate(rna_anno_2ndRound_level_2 = factor(rna_anno_2ndRound_level_2,
        levels = color_palette_level_2$sub_cluster
    )) %>%
    rename(group = rna_anno_2ndRound_level_2) %>%
    mutate(disease = paste0(disease, "_level_2"))

cell_count_merged_res <- bind_rows(cell_count_class_merged, cell_count_type_merged) %>%
    mutate(
        disease = factor(disease, levels = c(
            "ALS_level_1", "ALS_level_2",
            "FTD_level_1", "FTD_level_2", "Control_level_1", "Control_level_2"
        )),
        region = factor(region, levels = c("MCX", "mFCX")),
        group = factor(group,
            levels = c(
                c("Exc_neuron", "Inh_neuron", "Non_neuron"),
                color_palette_level_2$sub_cluster
            )
        )
    )

group_colors <- c(cell_class_color, colors_level_2)

#
p <- cell_count_merged_res %>% ggplot(aes(disease, n))

p + geom_bar(aes(fill = group),
    color = "black",
    size = 0.1,
    width = 0.9,
    position = position_fill(),
    stat = "identity"
) +
    facet_wrap(~region,
        dir = "v",
        drop = TRUE,
        scales = "free_x",
        nrow = 1,
        strip.position = "bottom"
    ) +
    scale_fill_manual(
        values = group_colors,
        name = "Cell class/type",
        guide = guide_legend(ncol = 1)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    xlab("Group") +
    ylab("Proportion of types (%)") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key.size = unit(0.5, "lines"),
        legend.text = element_text(size = 6)
    )
ggsave("./plots/cell_type_count_proportion_merged_level_1_2.pdf",
    device = cairo_pdf(),
    width = 4, height = 2.5, useDingbats = FALSE
)

# statistical test of the neuronal cells proportion differences in FTD using individuals

cell_proportion <- cell_count_class %>%
    select(-orig.ident) %>%
    pivot_wider(names_from = rna_anno_2ndRound_level_1, values_from = n) %>%
    mutate(
        total_count = Exc_neuron + Inh_neuron + Non_neuron,
        neuron_proportion = (Exc_neuron + Inh_neuron) / total_count,
        exc_proportion = Exc_neuron / (Exc_neuron + Inh_neuron),
        disease = factor(disease, levels = c("ALS", "FTD", "Control")),
        region = factor(region, levels = c("MCX", "mFCX"))
    )

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")
disease_color_dark <- lightness(disease_color, scalefac(0.8)) %>%
    as.vector()
names(disease_color_dark) <- names(disease_color)

p <- cell_proportion %>% ggplot(aes(disease, 100 * neuron_proportion))
p +
    geom_boxplot(aes(color = disease), width = 0.5, outlier.shape = NA) +
    geom_beeswarm(aes(color = disease), size = 1, cex = 4) +
    facet_wrap(~region, nrow = 1, strip.position = "bottom") +
    stat_compare_means(
        aes(label = ..p.signif..),
        method = "t.test",
        comparisons = list(
            c("ALS", "Control"),
            c("FTD", "Control"),
            c("ALS", "FTD")
        )
    ) +
    xlab("Diagnosis") +
    ylab("Proportion of neurons (%)") +
    coord_cartesian(ylim = c(0, 50)) +
    scale_color_manual(values = disease_color_dark, name = "diagnosis") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none"
    )

ggsave("./plots/compare_neurons_proportion_between_disease_diagnosis.pdf",
    device = cairo_pdf(), width = 2, height = 2, useDingbats = FALSE
)

p <- cell_proportion %>% ggplot(aes(disease, 100 * exc_proportion))
p +
    geom_boxplot(aes(color = disease), width = 0.5, outlier.shape = NA) +
    geom_beeswarm(aes(color = disease), size = 1, cex = 4) +
    facet_wrap(~region, nrow = 1, strip.position = "bottom") +
    stat_compare_means(
        aes(label = ..p.signif..),
        method = "t.test",
        comparisons = list(
            c("ALS", "Control"),
            c("FTD", "Control"),
            c("ALS", "FTD")
        )
    ) +
    xlab("Diagnosis") +
    ylab("Proportion of excitatory neurons in all neurons (%)") +
    coord_cartesian(ylim = c(53, 90)) +
    scale_color_manual(values = disease_color_dark, name = "diagnosis") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none"
    )

ggsave("./plots/compare_ExcNeurons_proportion_between_disease_diagnosis.pdf",
    device = cairo_pdf(), width = 2, height = 2, useDingbats = FALSE
)

# show distibution of cell types in each sample
# mask donor ID
donor_id_masked <- c(
    "110" = "A1",
    "111" = "A2",
    "113" = "A3",
    "332" = "A4",
    "388" = "A5",
    "52" = "A6",
    "902" = "C1",
    "904" = "C2",
    "906" = "C3",
    "1069" = "C4",
    "91" = "C5",
    "945" = "C6",
    "36" = "F1",
    "54" = "F2",
    "55" = "F3",
    "674" = "F4",
    "61" = "F5"
)

p <- meta_data %>%
    mutate(
        rna_anno_2ndRound_level_2 = factor(rna_anno_2ndRound_level_2, levels = names(colors_level_2)),
        disease = factor(disease, levels = c("ALS", "FTD", "Control")),
        donor_id = donor_id_masked[subject]
    ) %>%
    arrange(disease, region, donor_id) %>%
    mutate(orig.ident = factor(orig.ident, levels = unique(orig.ident))) %>%
    ggplot(aes(orig.ident))
p +
    geom_bar(aes(fill = rna_anno_2ndRound_level_2), position = position_fill()) +
    scale_fill_manual(values = colors_level_2) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

# log session info
sessionInfo()
