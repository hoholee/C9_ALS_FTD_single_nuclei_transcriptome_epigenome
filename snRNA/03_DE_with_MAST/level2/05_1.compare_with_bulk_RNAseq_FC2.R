# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)

snRNA <- read_tsv("MAST_res_level2_summary.txt")

FDR_threshold <- 0.05
# FC_threshold <- 1.2
FC_threshold <- 2

snRNA_cell_types <- c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep",
    "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_ADARB2_Other",
    "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC"
)

snRNA_sig <- snRNA %>%
    filter(
        FDR < FDR_threshold,
        abs(model_log2FC) > log2(FC_threshold),
        conv_C == TRUE,
        conv_D == TRUE,
        model_log2FC_ci_hi * model_log2FC_ci_low > 0,
        abs(model_log2FC - avg_logFC) < 2
    )

bulk_RNA <- read_tsv("edgeR_res_glmTest_allExpressedGenes_rmOutlierSamples_simplified.txt")

plot_FC_cor <- function(snRNA_cell_type, bulkRNA_cell_type,
                        snRNA_cond_1, snRNA_cond_2, bulkRNA_comparison,
                        snRNA_region, bulkRNA_region) {
    message(str_glue(
        "Processing snRNA {snRNA_cell_type} vs. bulk RNA {bulkRNA_cell_type} in {snRNA_region} ",
        "with FC > {FC_threshold} and FDR < {FDR_threshold}, {snRNA_cond_1} vs. {snRNA_cond_2}"
    ))
    snRNA_sub <- snRNA_sig %>%
        filter(
            cell_type == snRNA_cell_type,
            cond_1 == snRNA_cond_1,
            cond_2 == snRNA_cond_2,
            region == snRNA_region
        )
    if (nrow(snRNA_sub) == 0) {
        message(str_glue("No significant DE genes in snRNA {snRNA_cell_type}!"))
    } else {
        bulk_RNA_sub <- bulk_RNA %>%
            filter(
                cell_type == bulkRNA_cell_type,
                comparison == bulkRNA_comparison,
                region == bulkRNA_region
            )
        res <- snRNA_sub %>%
            left_join(bulk_RNA_sub, by = c("gene" = "geneName")) %>%
            mutate(bulkRNA_sig = if_else(FDR.y < 0.05, "yes", "no")) %>%
            arrange(bulkRNA_sig)
        cor <- cor.test(res$model_log2FC, res$logFC, method = "spearman")
        rho <- round(cor$estimate, 3)
        p_value <- signif(cor$p.value, 3)
        out_stat_df <- tibble(
            snRNA_cond_1 = snRNA_cond_1,
            snRNA_cond_2 = snRNA_cond_2,
            snRNA_cell_type = snRNA_cell_type,
            bulkRNA_cell_type = bulkRNA_cell_type,
            region = snRNA_region,
            FC_threshold = FC_threshold,
            FDR_threshold = FDR_threshold,
            rho = cor$estimate,
            p_value = cor$p.value
        )
        write_tsv(out_stat_df, str_glue("./plot_vs_bulkRNA/{snRNA_cond_1}_vs_{snRNA_cond_2}_snRNA_{snRNA_cell_type}_expFC_vs_bulkRNA_{bulkRNA_cell_type}_in_{snRNA_region}_FC_{FC_threshold}_FDR_{FDR_threshold}_stat.txt"))
        p <- res %>% ggplot(aes(model_log2FC, logFC))
        p1 <- p +
            geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
            geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
            geom_point(aes(color = bulkRNA_sig), size = 0.1) +
            scale_color_manual(values = c("yes" = "#CD2526", "no" = "#000000")) +
            xlab(str_glue("snRNA gene expresssion log2FC\n({snRNA_cond_1} vs. {snRNA_cond_2})")) +
            ylab(str_glue("bulk RNA gene expression log2FC\n({snRNA_cond_1} vs. {snRNA_cond_2})")) +
            ggtitle(str_glue(
                "snRNA {snRNA_cell_type} vs. bulk RNA {bulkRNA_cell_type} in {snRNA_region}\n",
                "FC > {FC_threshold}, FDR < {FDR_threshold}\n",
                "rho = {rho}, p = {p_value}"
            )) +
            theme_bw(base_size = 8, base_family = "Helvetica") +
            theme(panel.grid.minor = element_blank())
        ggsave(str_glue("./plot_vs_bulkRNA/{snRNA_cond_1}_vs_{snRNA_cond_2}_snRNA_{snRNA_cell_type}_expFC_vs_bulkRNA_{bulkRNA_cell_type}_in_{snRNA_region}_FC_{FC_threshold}_FDR_{FDR_threshold}.pdf"),
            plot = p1, device = cairo_pdf(), width = 4, height = 3, useDingbats = FALSE
        )
        dev.off()
    }
}

# ALS vs. Control
plot_FC_cor("Astro", "Other_glias", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Astro", "Other_glias", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Micro", "Other_glias", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Micro", "Other_glias", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Oligo", "Oligodendrocytes", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Oligo", "Oligodendrocytes", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("OPC", "Oligodendrocytes", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("OPC", "Oligodendrocytes", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Exc_superficial", "Neurons", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Exc_superficial", "Neurons", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Exc_intermediate", "Neurons", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Exc_intermediate", "Neurons", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Exc_deep", "Neurons", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Exc_deep", "Neurons", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Inh_PVALB", "Neurons", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Inh_PVALB", "Neurons", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Inh_SST", "Neurons", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Inh_SST", "Neurons", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Inh_VIP", "Neurons", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Inh_VIP", "Neurons", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Inh_LAMP5", "Neurons", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Inh_LAMP5", "Neurons", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Inh_ADARB2_Other", "Neurons", "ALS", "Control", "ALS_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Inh_ADARB2_Other", "Neurons", "ALS", "Control", "ALS_vs_Control", "mFCX", "mid_frontal_cortex")

# FTD vs. Control
plot_FC_cor("Astro", "Other_glias", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Astro", "Other_glias", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Micro", "Other_glias", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Micro", "Other_glias", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Oligo", "Oligodendrocytes", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Oligo", "Oligodendrocytes", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("OPC", "Oligodendrocytes", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("OPC", "Oligodendrocytes", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Exc_superficial", "Neurons", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Exc_superficial", "Neurons", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Exc_intermediate", "Neurons", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Exc_intermediate", "Neurons", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Exc_deep", "Neurons", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Exc_deep", "Neurons", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Inh_PVALB", "Neurons", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Inh_PVALB", "Neurons", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Inh_SST", "Neurons", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Inh_SST", "Neurons", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

# plot_FC_cor("Inh_VIP", "Neurons", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex") # not enough finite observations
plot_FC_cor("Inh_VIP", "Neurons", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

plot_FC_cor("Inh_LAMP5", "Neurons", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
# plot_FC_cor("Inh_LAMP5", "Neurons", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex") # not enough finite observations

plot_FC_cor("Inh_ADARB2_Other", "Neurons", "FTD", "Control", "FTD_vs_Control", "MCX", "motor_cortex")
plot_FC_cor("Inh_ADARB2_Other", "Neurons", "FTD", "Control", "FTD_vs_Control", "mFCX", "mid_frontal_cortex")

# log session info
sessionInfo()