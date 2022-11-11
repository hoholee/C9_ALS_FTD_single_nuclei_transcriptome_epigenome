# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(viridis)
library(scico)
library(ggrepel)
library(shades)
library(ggrastr)
library(ggbeeswarm)
library(circlize)
library(cowplot)

set.seed(666)
addArchRThreads(threads = 16)
addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep")

celltype_order <- c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep",
    "Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST",
    "Astro", "Micro", "Oligo", "OPC"
)

# add back sex and age information from sample metadata table
sample_meta <- read_tsv("./sample_metadata.txt") %>%
    select(sample_id, sex, age)

# extract metadata, remove "Mixed" clusters in level 2 annotations
metadata <- getCellColData(proj) %>%
    as_tibble() %>%
    filter(!grepl("Mixed", atac_anno_level_2)) %>%
    mutate(
        atac_anno_level_2 = factor(atac_anno_level_2, levels = rev(celltype_order))
    ) %>%
    left_join(sample_meta, by = c("Sample" = "sample_id"))

# get fragment sizes for each sample, return as data frame
fragment_size_df <- plotFragmentSizes(
    ArchRProj = proj,
    groupBy = "Sample",
    chromSizes = getChromSizes(ArchRProj),
    maxSize = 750,
    returnDF = TRUE
)
# save
fragment_size <- fragment_size %>% as_tibble()
write_tsv(fragment_size, "fragment_size_by_sample.txt")

# plot fragment size distribution
individual_color <- c(
    "332" = "#92072a",
    "111" = "#c30938",
    "113" = "#f40b45",
    "52" = "#f63c6b",
    "388" = "#f86d90",
    "110" = "#fa9eb5",
    "36" = "#ffd866",
    "54" = "#ffcb33",
    "55" = "#ffbe00",
    "674" = "#cc9800",
    "61" = "#997200",
    "945" = "#00998d",
    "902" = "#00ccbb",
    "91" = "#00ffea",
    "1069" = "#33ffee",
    "904" = "#66fff3",
    "906" = "#99fff7"
)

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

region_rename <- c("MCX" = "Motor", "mFCX" = "Frontal")

p <- fragment_size %>%
    mutate(sample = group) %>%
    separate(group, c("region", "diagnosis", "individual"), sep = "_") %>%
    mutate(donor_id = donor_id_masked[individual]) %>%
    arrange(diagnosis, region, donor_id) %>%
    mutate(
        region_full = region_rename[region],
        group = paste(donor_id, region_full, sep = ", ")
    ) %>%
    mutate(group = factor(group, levels = unique(group))) %>%
    ggplot(aes(fragmentSize, fragmentPercent))

p + geom_line(aes(color = individual)) +
    facet_wrap(~group, ncol = 6) +
    scale_color_manual(values = individual_color) +
    xlab("snATAC-seq fragment Size (bp)") +
    ylab("Percentage of fragments") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
    )
ggsave(
    "./plots/snATAC_fragment_size_distribution_by_sample.pdf",
    device = cairo_pdf(),
    width = 4.5,
    height = 4,
    useDingbats = FALSE
)

# plot meta data side by side
celltype_color_palette <- read_tsv("./color_palette_level_2.txt")
celltype_color <- celltype_color_palette$color
names(celltype_color) <- celltype_color_palette$sub_cluster

region_labels <- c("MCX" = "motor cortex", "mFCX" = "medial frontal cortex")
region_color <- c("MCX" = "#432266", "mFCX" = "#FAA51B")
relabel_region <- function(x) {
    region_labels[x]
}

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")

sex_labels <- c("F" = "female", "M" = "male")
sex_color <- c("F" = "#FC9E89", "M" = "#3D304C")
relabel_sex <- function(x) {
    sex_labels[x]
}

relabel_donor <- function(x) {
    donor_id_masked[x]
}

p <- metadata %>%
    mutate(subject = factor(subject, levels = names(donor_id_masked))) %>%
    ggplot(aes(atac_anno_level_2))

p_cellcount <- p + geom_bar(stat = "count", aes(fill = atac_anno_level_2)) +
    scale_fill_manual(values = celltype_color) +
    scale_y_log10() +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    coord_flip() +
    annotation_logticks(sides = "b") +
    theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
    )

p_region <- p + geom_bar(aes(fill = region), stat = "count", position = position_fill()) +
    scale_fill_manual(name = "brain region", labels = relabel_region, values = region_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    coord_flip() +
    theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
    )

p_disease <- p + geom_bar(aes(fill = disease), stat = "count", position = position_fill()) +
    scale_fill_manual(values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    coord_flip() +
    theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
    )

p_sex <- p + geom_bar(aes(fill = sex), stat = "count", position = position_fill()) +
    scale_fill_manual(labels = relabel_sex(), values = sex_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    coord_flip() +
    theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
    )

p_individual <- p + geom_bar(aes(fill = subject), stat = "count", position = position_fill()) +
    scale_fill_manual(labels = relabel_donor(), values = individual_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    coord_flip() +
    theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
    )

p_n_unique_fragment <- p + geom_violin(aes(fill = atac_anno_level_2, y = log10(nFrags)), scale = "width") +
    scale_fill_manual(values = celltype_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    coord_flip() +
    theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
    )

p_TSS_enrichment <- p + geom_violin(aes(fill = atac_anno_level_2, y = TSSEnrichment), scale = "width") +
    scale_fill_manual(values = celltype_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    coord_flip() +
    theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
    )

plot_grid(p_cellcount, p_region, p_disease, p_sex, p_individual, p_n_unique_fragment, p_TSS_enrichment,
    align = "h", nrow = 1
)
ggsave(
    "./plots/atac_anno_level2_meta.pdf",
    device = cairo_pdf(),
    width = 4.5,
    height = 3,
    useDingbats = FALSE
)

# plot distribution of TSS enrichment and number of unique fragments
# by major cell types
p <- metadata %>%
    mutate(
        atac_anno_level_2 = factor(atac_anno_level_2, levels = celltype_order),
        disease = factor(disease, levels = c("ALS", "FTD", "Control"))
    ) %>%
    ggplot(aes(atac_anno_level_2, TSSEnrichment))

p +
    geom_hline(yintercept = 4, linetype = 2, color = "darkgrey") +
    geom_boxplot(aes(fill = disease), outlier.shape = NA) +
    scale_fill_manual(name = "Diagnosis", values = disease_color) +
    scale_y_continuous(breaks = seq(4, 20, length.out = 5)) +
    coord_cartesian(ylim = c(4, 20)) +
    xlab("Cell types") +
    ylab("TSS enrichment score") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

ggsave(
    "./plots/atac_TSS_enrichment_by_level2_anno_by_diagnosis.pdf",
    device = cairo_pdf(),
    width = 4,
    height = 2,
    useDingbats = FALSE
)

#
p <- metadata %>%
    mutate(
        atac_anno_level_2 = factor(atac_anno_level_2, levels = celltype_order),
        disease = factor(disease, levels = c("ALS", "FTD", "Control"))
    ) %>%
    ggplot(aes(atac_anno_level_2, log10(nFrags)))

p +
    geom_hline(yintercept = 3, linetype = 2, color = "darkgrey") +
    geom_boxplot(aes(fill = disease), outlier.shape = NA) +
    scale_fill_manual(name = "Diagnosis", values = disease_color) +
    xlab("Cell types") +
    ylab("Number of unique fragments") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

ggsave(
    "./plots/atac_num_unique_fragments_by_level2_anno_by_diagnosis.pdf",
    device = cairo_pdf(),
    width = 4,
    height = 2,
    useDingbats = FALSE
)


# by samples
sample_relabel <- c(
    "A1", "A2", "A3", "A4", "A5", "A6",
    "F1", "F2", "F3", "F4", "F5",
    "C1", "C2", "C3", "C4", "C5", "C6",
    "A1", "A2", "A3", "A4", "A5", "A6",
    "F1", "F2", "F3", "F4", "F5",
    "C1", "C2", "C3", "C4", "C5", "C6"
)
names(sample_relabel) <- c(
    "MCX_ALS_110", "MCX_ALS_111", "MCX_ALS_113", "MCX_ALS_332", "MCX_ALS_388", "MCX_ALS_52",
    "MCX_FTD_36", "MCX_FTD_54", "MCX_FTD_55", "MCX_FTD_674", "MCX_FTD_61",
    "MCX_Control_902", "MCX_Control_904", "MCX_Control_906", "MCX_Control_1069", "MCX_Control_91", "MCX_Control_945",
    "mFCX_ALS_110", "mFCX_ALS_111", "mFCX_ALS_113", "mFCX_ALS_332", "mFCX_ALS_388", "mFCX_ALS_52",
    "mFCX_FTD_36", "mFCX_FTD_54", "mFCX_FTD_55", "mFCX_FTD_674", "mFCX_FTD_61",
    "mFCX_Control_902", "mFCX_Control_904", "mFCX_Control_906", "mFCX_Control_1069", "mFCX_Control_91", "mFCX_Control_945"
)

p <- metadata %>%
    mutate(
        atac_anno_level_2 = factor(atac_anno_level_2, levels = celltype_order),
        disease = factor(disease, levels = c("ALS", "FTD", "Control")),
        subject = factor(subject,
            levels = c(
                "110", "111", "113", "332", "388", "52",
                "36", "54", "55", "674", "61",
                "902", "904", "906", "1069", "91", "945"
            )
        ),
        Sample = factor(Sample, levels = names(sample_relabel))
    ) %>%
    ggplot(aes(Sample, log10(nFrags)))

p +
    geom_hline(yintercept = 3, linetype = 2, color = "darkgrey") +
    geom_boxplot(aes(fill = disease), outlier.shape = NA) +
    scale_x_discrete(labels = sample_relabel) +
    scale_fill_manual(name = "Diagnosis", values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())
coord_cartesian(ylim = c(0))

ggsave(
    "./plots/atac_num_unique_fragments_by_sample.pdf",
    device = cairo_pdf(),
    width = 4.5,
    height = 1.5,
    useDingbats = FALSE
)

#
p <- metadata %>%
    mutate(
        atac_anno_level_2 = factor(atac_anno_level_2, levels = celltype_order),
        disease = factor(disease, levels = c("ALS", "FTD", "Control")),
        subject = factor(subject,
            levels = c(
                "110", "111", "113", "332", "388", "52",
                "36", "54", "55", "674", "61",
                "902", "904", "906", "1069", "91", "945"
            )
        ),
        Sample = factor(Sample, levels = names(sample_relabel))
    ) %>%
    ggplot(aes(Sample, TSSEnrichment))

p +
    geom_hline(yintercept = 4, linetype = 2, color = "darkgrey") +
    geom_boxplot(aes(fill = disease), outlier.shape = NA) +
    scale_x_discrete(labels = sample_relabel) +
    scale_y_continuous(breaks = seq(4, 20, length.out = 5)) +
    coord_cartesian(ylim = c(4, 24)) +
    scale_fill_manual(name = "Diagnosis", values = disease_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())

ggsave(
    "./plots/atac_TSS_enrichment_by_sample.pdf",
    device = cairo_pdf(),
    width = 4.5,
    height = 1.5,
    useDingbats = FALSE
)

# session info
sessionInfo()