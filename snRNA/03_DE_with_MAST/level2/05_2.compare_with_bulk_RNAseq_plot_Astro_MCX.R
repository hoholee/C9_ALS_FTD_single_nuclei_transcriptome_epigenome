# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ggrastr)
library(ggrepel)

snRNA <- read_tsv("MAST_res_level2_summary.txt")

FDR_threshold <- 0.05
FC_threshold <- 1.2
# FC_threshold <- 2

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


# ALS vs. Control
snRNA_cell_type <- "Astro"
bulkRNA_cell_type <- "Other_glias"
snRNA_cond_1 <- "ALS"
snRNA_cond_2 <- "Control"
bulkRNA_comparison <- "ALS_vs_Control"
snRNA_region <- "MCX"
bulkRNA_region <- "motor_cortex"

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

selected_genes <- c("GFAP", "CD44", "CHI3L1", "RANBP3L", "GOLGB1")
res_selected_genes <- res %>% filter(gene %in% selected_genes)

p <- res %>% ggplot(aes(model_log2FC, logFC))
p1 <- p +
    geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
    geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
    geom_point_rast(color = "#EF6075", raster.dpi = 600, scale = 0.5, alpha = 0.6) +
    geom_point(color = "black", size = 1, data = res_selected_genes, shape = 1) +
    geom_text_repel(data = res_selected_genes, aes(label = gene)) +
    xlab(str_glue("snRNA gene expresssion log2FC\n({snRNA_cond_1} vs. {snRNA_cond_2})")) +
    ylab(str_glue("bulk RNA gene expression log2FC\n({snRNA_cond_1} vs. {snRNA_cond_2})")) +
    ggtitle(str_glue(
        "snRNA {snRNA_cell_type} vs. bulk RNA {bulkRNA_cell_type} in {snRNA_region}\n",
        "FC > {FC_threshold}, FDR < {FDR_threshold}\n",
        "rho = {rho}, p = {p_value}"
    )) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())
ggsave(str_glue("./plot_vs_bulkRNA/{snRNA_cond_1}_vs_{snRNA_cond_2}_snRNA_{snRNA_cell_type}_expFC_vs_bulkRNA_{bulkRNA_cell_type}_in_{snRNA_region}_FC_{FC_threshold}_FDR_{FDR_threshold}_rasterized.pdf"),
    plot = p1, device = cairo_pdf(), width = 2.3, height = 2.3, useDingbats = FALSE
)

# log session info
sessionInfo()
