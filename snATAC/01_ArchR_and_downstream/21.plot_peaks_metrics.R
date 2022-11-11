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

marker_all <- readRDS("snATAC_marker_peaks_by_annoLevel2_clean.rds")

peak_meta <- getPeakSet(proj) %>% as_tibble()

p <- peak_meta %>% ggplot(aes(factor(Reproducibility)))

p +
    geom_bar(aes(fill = factor(Reproducibility)), width = 0.8) +
    scale_fill_viridis_d(
        option = "B",
        direction = 1,
        begin = 0.1,
        end = 0.9
    ) +
    xlab("Peak reproducibility") +
    ylab("Number of reproduciable peaks") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        legend.position = "none"
    )
ggsave(
    "./plots/num_reproducible_peaks_by_peak_reproducibility.pdf",
    device = cairo_pdf(),
    width = 1.2,
    height = 1.5,
    useDingbats = FALSE
)

p <- peak_meta %>%
    mutate(
        peakType = factor(
            peakType,
            levels = c("Promoter", "Exonic", "Intronic", "Distal")
        )
    ) %>%
    ggplot(aes(peakType))
p +
    geom_bar(aes(fill = peakType), width = 0.8) +
    scale_fill_viridis_d(
        option = "B",
        direction = 1,
        begin = 0.1,
        end = 0.9
    ) +
    xlab("Peak type") +
    ylab("Number of reproduciable peaks") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
ggsave(
    "./plots/num_reproducible_peaks_by_peak_type.pdf",
    device = cairo_pdf(),
    width = 1.2,
    height = 1.5,
    useDingbats = FALSE
)

p <- peak_meta %>%
    mutate(
        peakType = factor(
            peakType,
            levels = c("Promoter", "Exonic", "Intronic", "Distal")
        )
    ) %>%
    # ggplot(aes(log10(distToTSS+1)))
    ggplot(aes(log10(distToTSS)))

p +
    # geom_vline(xintercept = log10(c(2001, 51)), linetype = 2, color = "darkgrey") +
    # geom_vline(xintercept = log10(c(2000, 50)), linetype = 2, color = "darkgrey") +
    geom_density(aes(color = peakType)) +
    scale_color_viridis_d(
        option = "B",
        direction = 1,
        begin = 0.1,
        end = 0.9
    ) +
    xlab("Peak distance to TSS") +
    ylab("Density") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        # legend.position = c(0.85, 0.85)
        legend.position = "none"
    )
ggsave(
    "./plots/reproducible_peak_distanceToTSS_peak_type.pdf",
    device = cairo_pdf(),
    width = 1.5,
    height = 1.25,
    useDingbats = FALSE
)

# plot number of peaks in each group (cell type + diagnosis + region)
num_peaks <- read_tsv("./num_peaks_per_group.txt") %>%
    mutate(
        diagnosis = factor(diagnosis, levels = c("ALS", "FTD", "Control")),
        cell_type = factor(
            cell_type,
            levels = c(
                "Exc_superficial", "Exc_intermediate", "Exc_deep",
                "Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST",
                "Astro", "Micro", "Oligo", "OPC"
            )
        )
    )

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")

p <- num_peaks %>% ggplot(aes(cell_type, num_peaks / 10000))

p +
    geom_bar(
        aes(fill = diagnosis),
        stat = "identity",
        position = position_dodge2(),
        width = 0.75
    ) +
    facet_wrap(~region, nrow = 1) +
    scale_fill_manual(values = disease_color) +
    xlab("Cell type") +
    ylab("Number of peaks (10k)") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
ggsave(
    "./plots/num_reproducible_peaks_by_group.pdf",
    device = cairo_pdf(),
    width = 4,
    height = 2,
    useDingbats = FALSE
)


# session info
sessionInfo()