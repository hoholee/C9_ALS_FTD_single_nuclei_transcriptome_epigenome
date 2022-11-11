# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(pheatmap)
library(viridis)
library(scico)

df <- read_tsv("./snRNA_vs_bulkRNA_FC_correlation_stat_summary.txt")
region_relabel <- c(
    "MCX" = "motor cortex",
    "mFCX" = "frontal cortex"
)

plot_FC_correlation_heatmap <- function(selected_disease,
                                        selected_region,
                                        selected_FC_threshold) {
    df_sub <- df %>%
        filter(
            FC_threshold == selected_FC_threshold,
            snRNA_cond_1 == selected_disease,
            region == selected_region
        )

    rho_mat <- df_sub %>%
        select(snRNA_cell_type, bulkRNA_cell_type, rho) %>%
        pivot_wider(
            names_from = "bulkRNA_cell_type",
            values_from = "rho"
        ) %>%
        mutate(snRNA_cell_type = factor(snRNA_cell_type,
            levels = c(
                "Astro", "Micro", "Oligo", "OPC",
                "Exc_superficial", "Exc_intermediate", "Exc_deep",
                "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
            )
        )) %>%
        arrange(snRNA_cell_type) %>%
        select(snRNA_cell_type, Other_glias, Oligodendrocytes, Neurons) %>%
        as.data.frame() %>%
        column_to_rownames(var = "snRNA_cell_type")

    sig_mat <- df_sub %>%
        mutate(
            FDR = p.adjust(p_value, method = "fdr"),
            sig = case_when(
                FDR < 0.05 ~ "*",
                FDR < 0.01 ~ "**",
                FDR < 0.001 ~ "***",
                TRUE ~ ""
            )
        ) %>%
        select(snRNA_cell_type, bulkRNA_cell_type, sig) %>%
        pivot_wider(
            names_from = "bulkRNA_cell_type",
            values_from = "sig",
            values_fill = ""
        ) %>%
        mutate(snRNA_cell_type = factor(snRNA_cell_type,
            levels = c(
                "Astro", "Micro", "Oligo", "OPC",
                "Exc_superficial", "Exc_intermediate", "Exc_deep",
                "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other", "Inh_PVALB", "Inh_SST"
            )
        )) %>%
        arrange(snRNA_cell_type) %>%
        select(snRNA_cell_type, Other_glias, Oligodendrocytes, Neurons) %>%
        as.data.frame() %>%
        column_to_rownames(var = "snRNA_cell_type")

    pheatmap(rho_mat,
        color = scico(100, palette = "vikO", direction = 1),
        breaks = seq(-1, 1, length.out = 101),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        cellwidth = 8,
        cellheight = 8,
        fontsize = 8,
        show_rownames = TRUE,
        show_colnames = TRUE,
        border_color = NA,
        scale = "none",
        display_numbers = sig_mat,
        number_color = "#000000",
        main = str_glue(
            "Spearman correlation of snRNA and bulkRNA\n",
            "{selected_disease} vs. Control {region_relabel[selected_region]} expression FC\n",
            "in significnat DE genes (FDR < 0.05, FC > {selected_FC_threshold})"
        ),
        filename = str_glue(
            "./plot_vs_bulkRNA/snRNA_vs_bulkRNA_FC_spearman_cor_heatmap_",
            "{selected_disease}_vs_Control_{selected_region}_FC_{selected_FC_threshold}.pdf"
        ),
        width = 6,
        height = 4
    )
}

plot_FC_correlation_heatmap("ALS", "MCX", 2)
plot_FC_correlation_heatmap("ALS", "mFCX", 2)
plot_FC_correlation_heatmap("ALS", "MCX", 1.2)
plot_FC_correlation_heatmap("ALS", "mFCX", 1.2)

plot_FC_correlation_heatmap("FTD", "MCX", 2)
plot_FC_correlation_heatmap("FTD", "mFCX", 2)
plot_FC_correlation_heatmap("FTD", "MCX", 1.2)
plot_FC_correlation_heatmap("FTD", "mFCX", 1.2)

# log session info
sessionInfo()