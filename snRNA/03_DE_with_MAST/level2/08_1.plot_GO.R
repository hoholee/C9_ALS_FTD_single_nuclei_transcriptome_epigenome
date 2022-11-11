# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(scales)
library(viridis)
library(scico)
library(shades)

df <- read_tsv("./GO_webgestalt_results_affinity_rmRedundant_toPlot.tsv") %>%
    mutate(
        types = factor(types,
            levels = c(
                "MCX_Exc-superficial", "mFCX_Exc-superficial",
                "MCX_Exc-intermediate", "mFCX_Exc-intermediate",
                "MCX_Exc-deep", "mFCX_Exc-deep",
                "MCX_Inh-VIP", "mFCX_Inh-VIP",
                "MCX_Inh-LAMP5", "mFCX_Inh-LAMP5",
                "MCX_Inh-PVALB", "mFCX_Inh-PVALB",
                "MCX_Inh-SST", "mFCX_Inh-SST"
            )
        ),
        description = factor(description, levels = rev(c(
            "protein folding",
            "heat shock protein binding",
            "unfolded protein binding",
            "Golgi-associated vesicle",
            "respiratory chain",
            "mitochondrial inner membrane",
            "mitochondrial protein complex",
            "nucleoside triphosphate metabolic process",
            "NADH dehydrogenase complex",
            "NADH dehydrogenase complex assembly",
            "proton transmembrane transport",
            "mitochondrial membrane organization",
            "protein localization to endoplasmic reticulum",
            "polysome",
            "translational elongation",
            "protein localization to nucleus",
            "RNA localization",
            "DNA damage response, detection of DNA damage",
            "nucleotide-excision repair",
            "antibiotic metabolic process",
            "nucleobase-containing compound transport",
            "ribonucleoprotein complex biogenesis",
            "process utilizing autophagic mechanism",
            "microtubule cytoskeleton organization involved in mitosis",
            "tRNA metabolic process",
            "organelle fission",
            "DNA biosynthetic process"
        )))
    )

p <- df %>% ggplot(aes(types, description))

p +
    geom_point(aes(color = -log10(FDR + 1e-15), size = enrichmentRatio)) +
    scale_color_viridis_c(
        name = "-log10(FDR)",
        option = "B",
        direction = 1,
        begin = 0.1,
        end = 0.9,
        breaks = seq(2, 15, length.out = 5),
        label = c("<2", "5.25", "8.5", "11.75", ">15"),
        limits = c(-log10(0.01), -log10(1e-15)),
        oob = squish
    ) +
    scale_size_continuous(
        name = "Enrichment ratio",
        range = c(0, 3),
        limits = c(1, 12),
        breaks = seq(2, 10, length.out = 3),
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_bw(base_size = 6, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key.size = unit(0.5, "line")
    )
ggsave("./GO_dotplot_selected_terms.pdf",
    device = cairo_pdf(),
    width = 4, height = 3, useDingbats = FALSE
)