# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(viridis)
library(scico)
library(ggrepel)
library(broom)
library(future.apply)
library(ggridges)
library(shades)
library(ggrastr)
library(ggbeeswarm)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

set.seed(666)
addArchRThreads(threads = 16)
addArchRGenome("hg38")

proj <- loadArchRProject(path = "Dracheva_ALS_FTD_latest_atac_anno_addPseudobulkRep")

# extract the chromVAR deviation and Z score as a `SummarizedExperiment` object
# two assays: `deviations` and `z`
chromVAR_se <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "vierstra_motif_archetypeMatrix"
)
saveRDS(chromVAR_se, "chromVAR_vierstra_motif_archetype_deviations_and_z_in_summarizedExperiment.rds")
# chromVAR_se <- readRDS("chromVAR_vierstra_motif_archetype_deviations_and_z_in_summarizedExperiment.rds")

# test whether the chromVAR deviation or Z score is correlated with read depth
mat_z <- assays(chromVAR_se)$z
mat_deviations <- assays(chromVAR_se)$deviations

# test_df <- tibble(
#     z = mat_z[2,],
#     reads_in_TSS = chromVAR_se$ReadsInTSS
# )
# p <- test_df %>% ggplot(aes(reads_in_TSS, z))
# p + geom_hex(bins = 100) + theme_bw()

# test_df2 <- tibble(
#     deviation = mat_deviations[2,],
#     reads_in_TSS = chromVAR_se$ReadsInTSS
# )
# p <- test_df2 %>% ggplot(aes(reads_in_TSS, deviation))
# p + geom_hex(bins = 100) + theme_bw()

# read the variability (of z scores)
var_dev <- read_tsv("./vierstra_archetype_motif_chromVAR_variability_ranks.txt")

p <- var_dev %>% ggplot(aes(rank, combinedVars))
p +
    geom_point_rast(
        aes(color = combinedVars),
        size = 0.5,
        raster.dpi = 600,
        dev = "cairo",
        scale = 1
    ) +
    geom_text_repel(
        aes(label = name_clean),
        data = var_dev %>%
            filter(rank <= 10) %>%
            separate(name, c("archetype_id", "tf_name", "family_name"), sep = "\\|") %>%
            mutate(
                tf_name = str_replace_all(tf_name, "/", ", "),
                name_clean = str_glue("{tf_name} ({family_name})")
            ),
        nudge_x = 2,
        size = 1
    ) +
    scale_color_scico(palette = "lajolla", name = "Variability") +
    xlab("Rank sorted motif archetypes") +
    ylab("Variability of chromVAR Z scores") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "none"
    )
ggsave("./plots/chromVAR_rank_sorted_vierstra_motif_archetype_variability_cleanName.pdf",
    device = cairo_pdf(), width = 2, height = 2, useDingbats = FALSE
)

# add back sex and age information from sample metadata table
sample_meta <- read_tsv("./sample_metadata.txt") %>%
    select(sample_id, sex, age)
to_add_meta <- tibble(
    sample_id = chromVAR_se$Sample %>% as.character()
) %>%
    left_join(sample_meta, by = "sample_id")

# eta-square to estimate variance explained
# define a function to perform anova test
test_anova <- function(x) {
    df <- tibble(
        z = x,
        cluster = chromVAR_se$atac_anno_level_2 %>% as.character(),
        disease = chromVAR_se$disease %>% as.character(),
        region = chromVAR_se$region %>% as.character(),
        subject = chromVAR_se$subject %>% as.character(),
        sex = to_add_meta$sex
    )

    aov_res <- aov(z ~ cluster + disease + region + sex + cluster:disease, data = df)

    stats_summary <- glance(aov_res) %>%
        dplyr::rename(
            r_squared = r.squared
        )

    terms_summary <- tidy(aov_res) %>% mutate(eta_squared = sumsq / sum(sumsq))

    eta_sq <- terms_summary %>%
        select(term, eta_squared) %>%
        spread(term, eta_squared) %>%
        dplyr::select(cluster, disease, region, sex, `cluster:disease`, Residuals) %>%
        dplyr::rename(cluster_disease_interaction = `cluster:disease`) %>%
        rename_all(funs(paste0("eta_squared_", .)))

    terms_p <- terms_summary %>%
        select(term, p.value) %>%
        spread(term, p.value) %>%
        dplyr::select(cluster, disease, region, sex, `cluster:disease`) %>%
        dplyr::rename(cluster_disease_interaction = `cluster:disease`) %>%
        rename_all(funs(paste0("p_", .)))

    bind_cols(stats_summary, eta_sq, terms_p)
}

# use future version of `apply` to run in parallel
# availableCores()
# plan(multisession)
options(future.globals.maxSize = 2 * 1024^3)
plan(multicore, workers = 32)
anova_res <- future_apply(mat_z, 1, test_anova)

anova_res_df <- anova_res %>%
    bind_rows() %>%
    mutate(
        motif_id = names(anova_res),
        fdr_cluster = p.adjust(p_cluster, method = "fdr"),
        fdr_disease = p.adjust(p_disease, method = "fdr"),
        fdr_region = p.adjust(p_region, method = "fdr"),
        fdr_cluster_disease_interaction = p.adjust(p_cluster_disease_interaction, method = "fdr")
    ) %>%
    select(motif_id, everything())
write_tsv(anova_res_df, "chromVAR_z_vierstra_motif_archetype_anova_eta_square.txt")
# anova_res_df <- read_tsv("./chromVAR_z_vierstra_motif_archetype_anova_eta_square.txt")

p <- anova_res_df %>%
    select(motif_id, starts_with("eta_")) %>%
    pivot_longer(starts_with("eta_"), names_to = "term", values_to = "eta_sq") %>%
    mutate(term = str_replace(term, "eta_squared_", "")) %>%
    ggplot(aes(eta_sq, factor(term)))

p +
    geom_density_ridges_gradient(
        aes(fill = stat(x)),
        rel_min_height = 0.001,
        size = 0.3,
        scale = 1
    ) +
    scale_fill_viridis_c(name = "Eta squared", option = "C") +
    xlab("Eta squared") +
    ylab("Term") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )
ggsave("./plots/chromVAR_z_vierstra_motif_archetype_anova_eta_square.pdf",
    device = cairo_pdf(), width = 4, height = 3, useDingbats = FALSE
)

# hide the residuals, sort term and zoom in
p <- anova_res_df %>%
    select(motif_id, starts_with("eta_")) %>%
    select(-eta_squared_Residuals) %>%
    pivot_longer(starts_with("eta_"), names_to = "term", values_to = "eta_sq") %>%
    mutate(
        term = str_replace(term, "eta_squared_", ""),
        term_mod = case_when(
            term == "cluster" ~ "Cell type",
            term == "disease" ~ "Diagnosis",
            term == "region" ~ "Brain region",
            term == "sex" ~ "Sex",
            term == "cluster_disease_interaction" ~ "Interaction (Cell type and Diagnosis)",
        )
    ) %>%
    mutate(
        term_mod = factor(
            term_mod,
            levels = rev(c(
                "Cell type",
                "Diagnosis",
                "Interaction (Cell type and Diagnosis)",
                "Brain region",
                "Sex"
            ))
        )
    ) %>%
    ggplot(aes(eta_sq, term_mod))

p +
    geom_density_ridges_gradient(
        aes(fill = stat(x)),
        rel_min_height = 0.001,
        size = 0.3,
        scale = 1
    ) +
    # scale_fill_viridis_c(name = "Eta squared", option = "C") +
    scale_fill_scico(palette = "lajolla", name = "Eta squared") +
    xlab("Eta squared") +
    ylab("Term") +
    xlim(NA, 0.7) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )
ggsave("./plots/chromVAR_z_vierstra_motif_archetype_anova_eta_square_zoomed.pdf",
    device = cairo_pdf(), width = 3, height = 2, useDingbats = FALSE
)

## plot motifs with highest variability
# motif_id <- c("AC0259", "AC0227", "AC0622", "AC0610", "AC0606", "AC0261", "AC0489", "AC0405", "AC0024", "AC0525")
# marker_motifs <- getFeatures(proj, select = paste(motif_id, collapse = "|"), useMatrix = "vierstra_motif_archetypeMatrix")
# marker_motifs

# marker_motifs_z <- grep("z:", marker_motifs, value = TRUE)
# marker_motifs_z

# # density plot, grouped by `atac_anno_level_2`
# p <- plotGroups(
#     ArchRProj = proj,
#     groupBy = "atac_anno_level_2",
#     colorBy = "vierstra_motif_archetypeMatrix",
#     name = marker_motifs_z,
#     imputeWeights = getImputeWeights(proj),
#     # imputeWeights = NULL,
#     quantCut = c(0, 1),
#     maxCells = 100000
# )

# p2 <- lapply(seq_along(p), function(x) {
#     if (x != 1) {
#         p[[x]] + guides(color = FALSE, fill = FALSE) +
#             theme_ArchR(baseSize = 8) +
#             theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
#             theme(
#                 axis.text.y = element_blank(),
#                 axis.ticks.y = element_blank(),
#                 axis.title.y = element_blank()
#             ) + ylab("")
#     } else {
#         p[[x]] + guides(color = FALSE, fill = FALSE) +
#             theme_ArchR(baseSize = 8) +
#             theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
#             theme(
#                 axis.ticks.y = element_blank(),
#                 axis.title.y = element_blank()
#             ) + ylab("")
#     }
# })
# do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))), p2))

# plotPDF(p,
#     name = "Plot-Groups-Deviations-w-Imputation-by-atac-anno-level-2",
#     width = 6, height = 4.5, ArchRProj = proj, addDOC = FALSE
# )

# manually make the ridge plot to provide more flexibility on plotting styles
# todo: check the implementation of `maxCells` and `quantCut`; ignore them for now
impute_weight <- getImputeWeights(proj)
mat_z_imputed <- imputeMatrix(
    mat = as.matrix(mat_z),
    imputeWeights = impute_weight
)

# to_plot_motif_id <- marker_motifs_z %>% str_replace("z:", "")
# to_plot_motif_id <- var_dev$name[1:10]
to_plot_motif_id <- var_dev$name[c(1, 2, 5)]
# to_plot_value <- mat_z[to_plot_motif_id, , drop = F] %>%
#     t() %>%
#     as.data.frame() %>%
#     rownames_to_column("cell_id") %>%
#     as_tibble()
to_plot_value <- mat_z_imputed[to_plot_motif_id, , drop = F] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()
to_plot_meta <- getCellColData(proj) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()

to_plot_res <- to_plot_meta %>%
    left_join(to_plot_value, by = "cell_id") %>%
    select(cell_id, region, disease, subject, atac_anno_level_2, starts_with("AC")) %>%
    filter(!grepl("Mixed__", atac_anno_level_2)) %>%
    pivot_longer(starts_with("AC"), names_to = "motif", values_to = "chromVAR_z") %>%
    mutate(
        atac_anno_level_2 = factor(atac_anno_level_2,
            levels = rev(c(
                "Exc_superficial", "Exc_intermediate", "Exc_deep",
                "Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST",
                "Astro", "Micro", "Oligo", "OPC"
            ))
        ),
        motif = factor(motif, levels = to_plot_motif_id)
    )

level2_colors_palette <- read_tsv("./color_palette_level_2.txt")
level2_colors <- level2_colors_palette$color
names(level2_colors) <- level2_colors_palette$sub_cluster

p <- to_plot_res %>% ggplot(aes(y = atac_anno_level_2, x = chromVAR_z))
p +
    geom_density_ridges(
        aes(fill = atac_anno_level_2),
        rel_min_height = 0.001,
        size = 0.3,
        scale = 2,
    ) +
    facet_wrap(~motif, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = level2_colors) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )
ggsave("./plots/chromVAR_z_vierstra_motif_archetype_selected_motif_Z_distribution_in_level2_anno.pdf",
    device = cairo_pdf(), width = 4, height = 2, useDingbats = FALSE
)

# make a heatmap of average chromVAR Z score in level 2 cell types, grouped by diagnosis
mat_z_imputed_sp <- as(mat_z_imputed, "dgCMatrix")
col_order <- colnames(mat_z_imputed_sp)

ordered_meta <- tibble(
    cell_id = col_order
) %>%
    left_join(to_plot_meta, by = "cell_id")

all.equal(ordered_meta$cell_id, colnames(mat_z_imputed_sp))

group <- ordered_meta$group %>% as_factor()
names(group) <- ordered_meta$cell_id

mm <- model.matrix(~ 0 + group)
colnames(mm) <- levels(group)

mat_z_imputed_sum <- mat_z_imputed_sp %*% mm

group_count <- ordered_meta %>% dplyr::count(group)
count_vec <- group_count$n
names(count_vec) <- group_count$group
count_vec <- count_vec[colnames(mat_z_imputed_sum)]

all.equal(colnames(mat_z_imputed_sum), names(count_vec))

mat_z_imputed_avg <- mat_z_imputed_sum / count_vec[col(mat_z_imputed_sum)]

to_plot_motif_id <- var_dev$name[1:100]
# to_plot_motif_id <- var_dev$name[1:10]
# to_plot_motif_id <- var_dev$name[1:20]
# to_plot_motif_id <- var_dev$name

# to_plot_motif_id <- c(
#     "AC0239|BACH/NFE|bZIP",
#     "AC0244|MAFG/MAFF|bZIP",
#     "AC0242|FOSL/JUND|bZIP"
#     )

# to_plot_motif_id <- diff_deviation_cluster_t %>%
#     filter(fdr < 0.05) %>%
#     group_by(cluster) %>%
#     arrange(-mean_x) %>%
#     top_n(2, wt = mean_x) %>%
#     pull(motif_id)

to_plot_group <- ordered_meta %>%
    distinct(group) %>%
    filter(!grepl("Mixed__", group)) %>%
    arrange(group) %>%
    pull(group)

# to_plot_group <- c(
#     "Astro_MCX_ALS", "Astro_MCX_Control", "Astro_MCX_FTD",
#     "Astro_mFCX_ALS", "Astro_mFCX_Control", "Astro_mFCX_FTD"
# )

to_plot_mat <- mat_z_imputed_avg[to_plot_motif_id, to_plot_group, drop = FALSE]

# pheatmap(
#     to_plot_mat,
#     color = scico(20, palette = "vik"),
#     # color = scico(20, palette = "lajolla"),
#     breaks = seq(from = -5, to = 5, length.out = 21),
#     # breaks = seq(from = -10, to = 10, length.out = 21),
#     cluster_rows = TRUE,
#     cluster_cols = FALSE,
#     scale = "none",
#     show_rownames = TRUE
#     # show_rownames = FALSE
# )

to_plot_mat_full <- as.matrix(to_plot_mat)

color_fun <- colorRamp2(seq(-5, 5, length.out = 21), scico(n = 21, palette = "vik"))
mat_column_order <- paste(
    rep(c(
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST",
        "Astro", "Micro", "Oligo", "OPC"
    ), each = 6),
    rep(rep(c("MCX", "mFCX"), each = 3), 11),
    rep(c("ALS", "FTD", "Control"), 22),
    sep = "_"
)

selected_cell_types <- c(
    "Exc_superficial", "Exc_intermediate", "Exc_deep",
    "Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST",
    "Astro", "Micro", "Oligo", "OPC"
)

region_color <- c("MCX" = "#432266", "mFCX" = "#FAA51B")
cell_type_color_df <- read_tsv("./color_palette_level_2.txt") %>%
    filter(sub_cluster %in% selected_cell_types) %>%
    mutate(sub_cluster = factor(sub_cluster, levels = selected_cell_types))
cell_type_color <- cell_type_color_df$color
names(cell_type_color) <- cell_type_color_df$sub_cluster

cell_class_color <- c(
    "Non_neuron" = "#9C482B",
    "Exc_neuron" = "#96BB45",
    "Inh_neuron" = "#718DC7"
)

disease_color <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")

col_meta <- tibble(
    group = colnames(to_plot_mat),
    region = if_else(grepl("_MCX_", group), "MCX", "mFCX"),
    cell_type = str_extract(colnames(to_plot_mat), pattern = "\\w+(?=_(MCX|mFCX)_)"),
    diagnosis = case_when(
        grepl("ALS", group) ~ "ALS",
        grepl("FTD", group) ~ "FTD",
        grepl("Control", group) ~ "Control"
    ),
    cell_class = case_when(
        grepl("Exc_", cell_type) ~ "Exc_neuron",
        grepl("Inh_", cell_type) ~ "Inh_neuron",
        TRUE ~ "Non_neuron"
    )
)

mat_col_anno <- HeatmapAnnotation(
    region = col_meta$region,
    cell_class = col_meta$cell_class,
    cell_type = col_meta$cell_type,
    diagnosis = col_meta$diagnosis,
    col = list(
        region = region_color,
        cell_class = cell_class_color,
        cell_type = cell_type_color,
        diagnosis = disease_color
    )
)

top_motifs_each_cluster <- read_tsv("./chromVAR_diff_deviations_by_atac_anno_level2_t_vierstra_motif_archetype_top3_each_cluster.txt") %>%
    mutate(motif_id2 = motif_id) %>%
    separate(motif_id2, c("archetype_id", "tf_name", "family_name"), sep = "\\|") %>%
    mutate(
        tf_name = str_replace_all(tf_name, "/", ", "),
        name_clean = str_glue("{tf_name} ({family_name})")
    )
mat_row_anno <- rowAnnotation(
    motif = anno_mark(
        at = match(top_motifs_each_cluster$motif_id, rownames(to_plot_mat)),
        labels = top_motifs_each_cluster$name_clean,
        lines_gp = gpar(),
        labels_gp = gpar(fontsize = 5)
    )
)

pdf("./plots/chromVAR_diff_Z_by_disease_vierstra_motif_archetype_top_heatmap.pdf", width = 8, height = 8)
ht <- Heatmap(
    matrix = to_plot_mat_full,
    col = color_fun,
    name = "chromVAR Z",
    na_col = "grey",
    color_space = "LAB",
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    # row_order = mat_row_order,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    row_dend_reorder = TRUE,
    show_row_dend = TRUE,
    cluster_columns = FALSE,
    column_order = mat_column_order,
    top_annotation = mat_col_anno,
    right_annotation = mat_row_anno,
    use_raster = TRUE,
    raster_by_magick = TRUE,
    raster_magick_filter = "Lanczos",
    raster_device = "CairoPNG",
    width = unit(4, "inches"),
    height = unit(3, "inches"),
)
draw(ht)
dev.off()


# differential test on deviations or z-score
# adapted from `chromVAR`:
# https://github.com/GreenleafLab/chromVAR/blob/0f27fcc8463d537d770f164a55f949701beb6add/R/differential_tests.R#L18

# define helper funtions
# for deviations, t-test (parametric) or wilcoxon test (non-parametric)
t_helper_mod <- function(x, groups) {
    splitx <- split(x, groups)
    t_res <- t.test(splitx[[1]], splitx[[2]],
        alternative = "two.sided",
        paired = FALSE,
        var.equal = FALSE
    )
    tibble(
        p = t_res$p.value,
        mean_x = t_res$estimate[1],
        mean_y = t_res$estimate[2]
    )
}

wilcoxon_helper_mod <- function(x, groups) {
    splitx <- split(x, groups)
    wilcoxon_res <- wilcox.test(splitx[[1]], splitx[[2]],
        alternative = "two.sided",
        mu = 0,
        paired = FALSE,
        conf.int = TRUE
    )
    tibble(
        p = wilcoxon_res$p.value,
        diff_estimate = wilcoxon_res$estimate,
        diff_conf_int_low = wilcoxon_res$conf.int[1],
        diff_conf_int_hi = wilcoxon_res$conf.int[2]
    )
}

# for z-scores: Brown-Forsythe test
bf_var_test_mod <- function(x, groups) {
    medians <- aggregate(x, list(groups), median, na.rm = TRUE)$x
    median_diff <- abs(x - unsplit(medians, groups))
    anova_res <- anova(lm(median_diff ~ groups))
    splitx <- split(x, groups)
    mean_target <- mean(splitx[[1]])
    mean_ref <- mean(splitx[[2]])
    tibble(
        p = anova_res[1, 5],
        mean_x = mean_target,
        mean_y = mean_ref
    )
}


# testing between two groups:
# one cell type vs. the rest
# deviations diff
differentialDeviations_cluster <- function(selected_cluster, parametric = TRUE) {
    df <- tibble(
        cluster = chromVAR_se$atac_anno_level_2 %>% as.character(),
    ) %>%
        mutate(group = if_else(cluster == selected_cluster, selected_cluster, "other"))
    groups <- factor(df$group, levels = c(selected_cluster, "other"))
    if (parametric) {
        test_res <- apply(mat_deviations, 1, t_helper_mod, groups)
        bind_rows(test_res) %>%
            mutate(
                motif_id = names(test_res),
                mean_diff = mean_x - mean_y,
                fdr = p.adjust(p, method = "fdr"),
                cluster = selected_cluster
            ) %>%
            select(cluster, motif_id, mean_x, mean_y, mean_diff, p, fdr) %>%
            arrange(fdr)
    } else {
        test_res <- apply(mat_deviations, 1, wilcoxon_helper_mod, groups)
        bind_rows(test_res) %>%
            mutate(
                motif_id = names(test_res),
                fdr = p.adjust(p, method = "fdr"),
                cluster = selected_cluster
            ) %>%
            select(cluster, motif_id, diff_estimate, diff_conf_int_low, diff_conf_int_hi, p, fdr) %>%
            arrange(fdr)
    }
}

cluster_to_test <- chromVAR_se$atac_anno_level_2 %>%
    sort() %>%
    unique()

diff_deviation_cluster_t <- map_dfr(
    cluster_to_test,
    differentialDeviations_cluster,
    parametric = TRUE
)
write_tsv(diff_deviation_cluster_t, "chromVAR_diff_deviations_by_atac_anno_level2_t_vierstra_motif_archetype.txt")

# super slow... consider refactor the wilcoxon test part with `future_apply` later
# diff_deviation_cluster_wilcoxon <- map_dfr(
#     cluster_to_test,
#     differentialDeviations_cluster,
#     parametric = FALSE
# )

# z-score diff
differentialZ_cluster <- function(selected_cluster) {
    df <- tibble(
        cluster = chromVAR_se$atac_anno_level_2 %>% as.character(),
    ) %>%
        mutate(group = if_else(cluster == selected_cluster, selected_cluster, "other"))
    groups <- factor(df$group, levels = c(selected_cluster, "other"))
    test_res <- apply(mat_z, 1, bf_var_test_mod, groups)
    bind_rows(test_res) %>%
        mutate(
            motif_id = names(test_res),
            mean_diff = mean_x - mean_y,
            fdr = p.adjust(p, method = "fdr"),
            cluster = selected_cluster
        ) %>%
        select(cluster, motif_id, mean_x, mean_y, mean_diff, p, fdr) %>%
        arrange(fdr)
}

diff_z_cluster <- map_dfr(
    cluster_to_test,
    differentialZ_cluster
)
write_tsv(diff_z_cluster, "chromVAR_diff_z_by_atac_anno_level2_vierstra_motif_archetype.txt")


# testing between two groups:
# disease (ALS or FTD) vs. control, separated by brain region
# deviations diff
differentialDeviations_disease <- function(selected_cluster, selected_region, selected_disease) {
    idx <- which(
        chromVAR_se$atac_anno_level_2 == selected_cluster &
            chromVAR_se$region == selected_region &
            chromVAR_se$disease %in% c(selected_disease, "Control")
    )
    df <- tibble(
        disease = chromVAR_se$disease[idx]
    )
    groups <- factor(df$disease, levels = c(selected_disease, "Control"))
    inputs <- mat_deviations[, idx]
    test_res <- apply(inputs, 1, t_helper_mod, groups)
    bind_rows(test_res) %>%
        mutate(
            motif_id = names(test_res),
            mean_diff = mean_x - mean_y,
            fdr = p.adjust(p, method = "fdr"),
            cluster = selected_cluster,
            disease = selected_disease,
            region = selected_region
        ) %>%
        select(cluster, region, disease, motif_id, mean_x, mean_y, mean_diff, p, fdr) %>%
        arrange(fdr)
}

params <- expand_grid(
    selected_cluster = cluster_to_test,
    selected_region = c("MCX", "mFCX"),
    selected_disease = c("ALS", "FTD")
)
diff_dev_disease <- pmap_dfr(params, differentialDeviations_disease)

write_tsv(diff_dev_disease, "chromVAR_diff_deviations_by_disease_vierstra_motif_archetype.txt")
# diff_dev_disease <- read_tsv("./chromVAR_diff_deviations_by_disease_vierstra_motif_archetype.txt")

# z-score diff
differentialZ_disease <- function(selected_cluster, selected_region, selected_disease) {
    idx <- which(
        chromVAR_se$atac_anno_level_2 == selected_cluster &
            chromVAR_se$region == selected_region &
            chromVAR_se$disease %in% c(selected_disease, "Control")
    )
    df <- tibble(
        disease = chromVAR_se$disease[idx]
    )
    groups <- factor(df$disease, levels = c(selected_disease, "Control"))
    inputs <- mat_z[, idx]
    test_res <- apply(inputs, 1, bf_var_test_mod, groups)
    bind_rows(test_res) %>%
        mutate(
            motif_id = names(test_res),
            mean_diff = mean_x - mean_y,
            fdr = p.adjust(p, method = "fdr"),
            cluster = selected_cluster,
            disease = selected_disease,
            region = selected_region
        ) %>%
        select(cluster, region, disease, motif_id, mean_x, mean_y, mean_diff, p, fdr) %>%
        arrange(fdr)
}

diff_z_disease <- pmap_dfr(params, differentialZ_disease)

write_tsv(diff_z_disease, "chromVAR_diff_z_by_disease_vierstra_motif_archetype.txt")

# plot chromVAR differential results
# cluster-specific
top_each_cluster <- diff_deviation_cluster_t %>%
    filter(!grepl("Mixed", cluster)) %>%
    filter(mean_x > 0, mean_diff > 0, fdr < 0.05) %>%
    arrange(-mean_diff) %>%
    group_by(cluster) %>%
    dplyr::slice(1:3)

p <- diff_deviation_cluster_t %>%
    filter(!grepl("Mixed", cluster)) %>%
    ggplot(aes(mean_diff, -log10(fdr + 1e-50)))
p +
    geom_point(color = "grey50", size = 0.5) +
    geom_point(
        aes(color = cluster),
        size = 0.5,
        data = diff_deviation_cluster_t %>%
            filter(!grepl("Mixed", cluster)) %>%
            filter(fdr < 0.05, mean_x > 0, mean_diff > 0)
    ) +
    # geom_text_repel(aes(label = motif_id),
    #     data = top_each_cluster,
    #     size = 2
    # ) +
    facet_wrap(~cluster, scales = "free_x") +
    xlab("Mean difference in accessibility\n(chromVAR deviations, vs. remaining cells)") +
    ylab("-log10(FDR)") +
    #   scale_color_manual(values = cluster_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )

ggsave("./plots/chromVAR_diff_deviations_by_cluster_vierstra_motif_archetype_top_TF_volcano.pdf",
    device = cairo_pdf(), width = 6, height = 6, useDingbats = FALSE
)

top_each_cluster_z <- diff_z_cluster %>%
    filter(!grepl("Mixed", cluster)) %>%
    filter(mean_x > 0, mean_diff > 0, fdr < 0.05) %>%
    arrange(-mean_diff) %>%
    group_by(cluster) %>%
    dplyr::slice(1:3)

p <- diff_z_cluster %>%
    filter(!grepl("Mixed", cluster)) %>%
    ggplot(aes(mean_diff, -log10(fdr + 1e-50)))
p +
    geom_point(color = "grey50", size = 0.5) +
    geom_point(
        aes(color = cluster),
        size = 0.5,
        data = diff_z_cluster %>%
            filter(!grepl("Mixed", cluster)) %>%
            filter(mean_x > 0, mean_diff > 0, fdr < 0.05)
    ) +
    # geom_text_repel(
    #     aes(label = motif_id),
    #     data = top_each_cluster_z,
    #     size = 2
    # ) +
    facet_wrap(~cluster, scales = "free_x") +
    xlab("Mean difference in variability\n(chromVAR Z, vs. remaining cells)") +
    ylab("-log10(FDR)") +
    # scale_color_manual(values = cluster_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )
ggsave("./plots/chromVAR_diff_z_by_cluster_vierstra_motif_archetype_top_TF_volcano.pdf",
    device = cairo_pdf(), width = 6, height = 6, useDingbats = FALSE
)

# disease
top_each_disease <- diff_dev_disease %>%
    filter(!grepl("Mixed", cluster)) %>%
    filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05)) %>%
    arrange(-abs(mean_diff)) %>%
    group_by(cluster, region, disease) %>%
    dplyr::slice(1:3)

p <- diff_dev_disease %>%
    filter(!grepl("Mixed", cluster)) %>%
    ggplot(aes(mean_diff, -log10(fdr + 1e-50)))
p +
    geom_point(color = "grey50", size = 0.5) +
    geom_point(aes(color = cluster),
        size = 0.5,
        data = diff_dev_disease %>%
            filter(!grepl("Mixed", cluster)) %>%
            filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05))
    ) +
    # geom_text_repel(
    #     aes(label = motif_id),
    #     data = top_each_disease,
    #     size = 1.8
    # ) +
    facet_grid(disease + region ~ cluster, scales = "free_x") +
    xlab("Mean difference in accessibility\n(chromVAR deviations, disease vs. control)") +
    ylab("-log10(FDR)") +
    # scale_color_manual(values = cluster_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )
ggsave("./plots/chromVAR_diff_deviations_by_disease_vierstra_motif_archetype_top_TF_volcano.pdf",
    device = cairo_pdf(), width = 10, height = 6, useDingbats = FALSE
)

###

top_each_disease_z <- diff_z_disease %>%
    filter(!grepl("Mixed", cluster)) %>%
    filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05)) %>%
    arrange(-abs(mean_diff)) %>%
    group_by(cluster, region, disease) %>%
    dplyr::slice(1:3)

p <- diff_z_disease %>%
    filter(!grepl("Mixed", cluster)) %>%
    ggplot(aes(mean_diff, -log10(fdr + 1e-50)))
p +
    geom_point(color = "grey50", size = 0.5) +
    geom_point(aes(color = cluster),
        size = 0.5,
        data = diff_z_disease %>%
            filter(!grepl("Mixed", cluster)) %>%
            filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05))
    ) +
    # geom_text_repel(
    #     aes(label = motif_id),
    #     data = top_each_disease_z,
    #     size = 1.8
    # ) +
    facet_grid(disease + region ~ cluster, scales = "free") +
    xlab("Mean difference in variability\n(chromVAR Z, disease vs. control)") +
    ylab("-log10(FDR)") +
    #   scale_color_manual(values = cluster_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )
ggsave("./plots/chromVAR_diff_Z_by_disease_vierstra_motif_archetype_top_TF_volcano.pdf",
    device = cairo_pdf(), width = 10, height = 6, useDingbats = FALSE
)

# zoomed in Astro, MCX, ALS vs. Control
p <- diff_dev_disease %>%
    filter(cluster == "Astro", region == "MCX", disease == "ALS") %>%
    ggplot(aes(mean_diff, -log10(fdr + 1e-50)))
p +
    geom_point(color = "grey50", size = 0.5) +
    geom_point(aes(color = cluster),
        size = 0.5,
        data = diff_dev_disease %>%
            filter(cluster == "Astro", region == "MCX", disease == "ALS") %>%
            filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05))
    ) +
    geom_text_repel(
        aes(label = motif_id),
        data = top_each_disease %>%
            filter(cluster == "Astro", region == "MCX", disease == "ALS"),
        size = 1.8
    ) +
    facet_wrap(~cluster, scales = "free_x") +
    xlab("Mean difference in accessibility\n(chromVAR deviations, disease vs. control)") +
    ylab("-log10(FDR)") +
    # scale_color_manual(values = cluster_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )
ggsave("./plots/chromVAR_diff_deviations_by_disease_vierstra_motif_archetype_top_TF_volcano_Astro_MCX_ALS_vs_Control.pdf",
    device = cairo_pdf(), width = 4, height = 3, useDingbats = FALSE
)

p <- diff_z_disease %>%
    filter(cluster == "Astro", region == "MCX", disease == "ALS") %>%
    ggplot(aes(mean_diff, -log10(fdr + 1e-50)))
p +
    geom_point(color = "grey50", size = 0.5) +
    geom_point(aes(color = cluster),
        size = 0.5,
        data = diff_z_disease %>%
            filter(cluster == "Astro", region == "MCX", disease == "ALS") %>%
            filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05))
    ) +
    geom_text_repel(
        aes(label = motif_id),
        data = top_each_disease_z %>%
            filter(cluster == "Astro", region == "MCX", disease == "ALS"),
        size = 1.8
    ) +
    facet_wrap(~cluster, scales = "free_x") +
    xlab("Mean difference in variability\n(chromVAR Z, disease vs. control)") +
    ylab("-log10(FDR)") +
    #   scale_color_manual(values = cluster_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )
ggsave("./plots/chromVAR_diff_Z_by_disease_vierstra_motif_archetype_top_TF_volcano_Astro_MCX_ALS_vs_Control.pdf",
    device = cairo_pdf(), width = 4, height = 3, useDingbats = FALSE
)


# link snRNA expression FC between disease and control
# read snRNA MAST DE results (level 2)
snRNA_DE <- read_tsv("./differential_gene_scores_by_snATAC_anno/MAST_res_level2_summary.txt")
# read Vierstra motif archetype and TF match table
motif_TF_lookup <- read_tsv("./Vierstra_archetype_motif/metadata.tsv")

# test on Astro, MCX, ALS vs. Control
snRNA_DE_sub <- snRNA_DE %>%
    filter(
        cond_1 == "ALS",
        cond_2 == "Control",
        cell_type == "Astro",
        region == "MCX"
    ) %>%
    mutate(snRNA_sig = if_else(
        FDR < 0.05 &
            abs(model_log2FC) > log2(1.2) &
            conv_C == TRUE &
            conv_D == TRUE &
            model_log2FC_ci_hi * model_log2FC_ci_low > 0 &
            abs(model_log2FC - avg_logFC) < 2,
        "yes",
        "no"
    )) %>%
    select(gene, model_log2FC, FDR, avg_log2CPM_cond_1, avg_log2CPM_cond_2, snRNA_sig)

diff_dev_disease_sub <- diff_dev_disease %>%
    filter(
        cluster == "Astro",
        region == "MCX",
        disease == "ALS"
    ) %>%
    mutate(archetype_id = str_extract(motif_id, "AC\\d+")) %>%
    dplyr::rename(archetype_name = motif_id) %>%
    select(archetype_id, archetype_name, fdr, mean_x, mean_y, mean_diff)

diff_z_disease_sub <- diff_z_disease %>%
    filter(
        cluster == "Astro",
        region == "MCX",
        disease == "ALS"
    ) %>%
    mutate(archetype_id = str_extract(motif_id, "AC\\d+")) %>%
    dplyr::rename(archetype_name = motif_id) %>%
    select(archetype_id, archetype_name, fdr, mean_x, mean_y, mean_diff)

merged_df_Astro <- motif_TF_lookup %>%
    left_join(diff_dev_disease_sub, by = c("cluster" = "archetype_id")) %>%
    left_join(snRNA_DE_sub, by = c("tf_name" = "gene")) %>%
    drop_na() %>%
    mutate(sig_both = if_else(snRNA_sig == "yes" & fdr < 0.05, "yes", "no"))

p <- merged_df_Astro %>%
    ggplot(aes(model_log2FC, mean_diff))

p +
    geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
    geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
    geom_point(color = "grey", size = 0.5, data = merged_df_Astro %>% filter(sig_both == "no")) +
    geom_point(color = "#CD2526", size = 1, data = merged_df_Astro %>% filter(sig_both == "yes")) +
    geom_text_repel(
        aes(label = tf_name),
        data = merged_df_Astro %>%
            filter(sig_both == "yes") %>%
            select(tf_name, model_log2FC, mean_diff) %>%
            distinct(),
        size = 1.8
    ) +
    xlab("snRNA TF expresison log2FC\n(ALS vs. Control)") +
    ylab("snATAC Motif archetype chromVAR deviation difference\n(ALS vs. Control)") +
    ggtitle("Astro, motor cortex") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())
ggsave("./plots/chromVAR_diff_deviations_vs_snRNA_expFC_Astro_MCX_ALS_vs_Control.pdf",
    device = cairo_pdf(), width = 6, height = 4.5, useDingbats = FALSE
)


merged_df_z_Astro <- motif_TF_lookup %>%
    left_join(diff_z_disease_sub, by = c("cluster" = "archetype_id")) %>%
    left_join(snRNA_DE_sub, by = c("tf_name" = "gene")) %>%
    drop_na() %>%
    mutate(sig_both = if_else(snRNA_sig == "yes" & fdr < 0.05, "yes", "no"))

p <- merged_df_z_Astro %>%
    ggplot(aes(model_log2FC, mean_diff))

p +
    geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
    geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
    geom_point(color = "grey", size = 0.5, data = merged_df_z_Astro %>% filter(sig_both == "no")) +
    geom_point(color = "#CD2526", size = 1, data = merged_df_z_Astro %>% filter(sig_both == "yes")) +
    geom_text_repel(
        aes(label = tf_name),
        data = merged_df_z_Astro %>%
            filter(sig_both == "yes") %>%
            select(tf_name, model_log2FC, mean_diff) %>%
            distinct(),
        size = 1.8
    ) +
    xlab("snRNA TF expresison log2FC\n(ALS vs. Control)") +
    ylab("snATAC Motif archetype chromVAR Z-score difference\n(ALS vs. Control)") +
    ggtitle("Astro, motor cortex") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())
ggsave("./plots/chromVAR_diff_Z_vs_snRNA_expFC_Astro_MCX_ALS_vs_Control.pdf",
    device = cairo_pdf(), width = 6, height = 4.5, useDingbats = FALSE
)

# test on Exc superficial, MCX, ALS vs. Control
snRNA_DE_sub <- snRNA_DE %>%
    filter(
        cond_1 == "ALS",
        cond_2 == "Control",
        cell_type == "Exc_superficial",
        region == "MCX"
    ) %>%
    mutate(snRNA_sig = if_else(
        FDR < 0.05 &
            abs(model_log2FC) > log2(1.2) &
            conv_C == TRUE &
            conv_D == TRUE &
            model_log2FC_ci_hi * model_log2FC_ci_low > 0 &
            abs(model_log2FC - avg_logFC) < 2,
        "yes",
        "no"
    )) %>%
    select(gene, model_log2FC, FDR, avg_log2CPM_cond_1, avg_log2CPM_cond_2, snRNA_sig)

diff_dev_disease_sub <- diff_dev_disease %>%
    filter(
        cluster == "Exc_superficial",
        region == "MCX",
        disease == "ALS"
    ) %>%
    mutate(archetype_id = str_extract(motif_id, "AC\\d+")) %>%
    dplyr::rename(archetype_name = motif_id) %>%
    select(archetype_id, archetype_name, fdr, mean_x, mean_y, mean_diff)

diff_z_disease_sub <- diff_z_disease %>%
    filter(
        cluster == "Exc_superficial",
        region == "MCX",
        disease == "ALS"
    ) %>%
    mutate(archetype_id = str_extract(motif_id, "AC\\d+")) %>%
    dplyr::rename(archetype_name = motif_id) %>%
    select(archetype_id, archetype_name, fdr, mean_x, mean_y, mean_diff)

merged_df_ExcSuperficial <- motif_TF_lookup %>%
    left_join(diff_dev_disease_sub, by = c("cluster" = "archetype_id")) %>%
    left_join(snRNA_DE_sub, by = c("tf_name" = "gene")) %>%
    drop_na() %>%
    mutate(sig_both = if_else(snRNA_sig == "yes" & fdr < 0.05, "yes", "no"))

p <- merged_df_ExcSuperficial %>%
    ggplot(aes(model_log2FC, mean_diff))

p +
    geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
    geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
    geom_point(color = "grey", size = 0.5, data = merged_df_ExcSuperficial %>% filter(sig_both == "no")) +
    geom_point(color = "#CD2526", size = 1, data = merged_df_ExcSuperficial %>% filter(sig_both == "yes")) +
    geom_text_repel(
        aes(label = tf_name),
        data = merged_df_ExcSuperficial %>%
            filter(sig_both == "yes") %>%
            select(tf_name, model_log2FC, mean_diff) %>%
            distinct(),
        size = 1.8
    ) +
    xlab("snRNA TF expresison log2FC\n(ALS vs. Control)") +
    ylab("snATAC Motif archetype chromVAR deviation difference\n(ALS vs. Control)") +
    ggtitle("Exc superficial, motor cortex") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())
ggsave("./plots/chromVAR_diff_deviations_vs_snRNA_expFC_ExcSuperficial_MCX_ALS_vs_Control.pdf",
    device = cairo_pdf(), width = 6, height = 4.5, useDingbats = FALSE
)

merged_df_z_ExcSuperficial <- motif_TF_lookup %>%
    left_join(diff_z_disease_sub, by = c("cluster" = "archetype_id")) %>%
    left_join(snRNA_DE_sub, by = c("tf_name" = "gene")) %>%
    drop_na() %>%
    mutate(sig_both = if_else(snRNA_sig == "yes" & fdr < 0.05, "yes", "no"))

p <- merged_df_z_ExcSuperficial %>%
    ggplot(aes(model_log2FC, mean_diff))

p +
    geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
    geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
    geom_point(color = "grey", size = 0.5, data = merged_df_z_ExcSuperficial %>% filter(sig_both == "no")) +
    geom_point(color = "#CD2526", size = 1, data = merged_df_z_ExcSuperficial %>% filter(sig_both == "yes")) +
    geom_text_repel(
        aes(label = tf_name),
        data = merged_df_z_ExcSuperficial %>%
            filter(sig_both == "yes") %>%
            select(tf_name, model_log2FC, mean_diff) %>%
            distinct(),
        size = 1.8
    ) +
    xlab("snRNA TF expresison log2FC\n(ALS vs. Control)") +
    ylab("snATAC Motif archetype chromVAR Z-score difference\n(ALS vs. Control)") +
    ggtitle("Exc superficial, motor cortex") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())
ggsave("./plots/chromVAR_diff_Z_vs_snRNA_expFC_ExcSuperficial_MCX_ALS_vs_Control.pdf",
    device = cairo_pdf(), width = 6, height = 4.5, useDingbats = FALSE
)

# density plot to show disease effects of selected motif archetypes
disease_color_panel <- c("ALS" = "#F6416C", "FTD" = "#FFDE7D", "Control" = "#00B8A9")

plot_ridge_chromVAR_dev_by_disease <- function(selected_motifs, selected_cluster, selected_region, selected_disease) {
    idx <- which(
        chromVAR_se$atac_anno_level_2 == selected_cluster &
            chromVAR_se$region == selected_region &
            chromVAR_se$disease %in% c(selected_disease, "Control")
    )

    df <- tibble(
        disease = factor(chromVAR_se$disease[idx], levels = c(selected_disease, "Control"))
    )

    inputs <- mat_deviations[selected_motifs, idx]

    motif_data_df <- inputs %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("cell_id") %>%
        as_tibble() %>%
        mutate(
            disease = df$disease,
            region = selected_region,
            cluster = selected_cluster
        ) %>%
        pivot_longer(cols = -one_of(c("cell_id", "disease", "cluster", "region")), names_to = "motif_id", values_to = "dev_score")

    p <- motif_data_df %>% ggplot(aes(dev_score, disease))
    p +
        geom_density_ridges(
            aes(fill = disease, color = disease),
            rel_min_height = 0.001,
            size = 0.3,
            scale = 1,
            quantile_lines = TRUE,
            quantiles = 2
        ) +
        # geom_vline(xintercept = 0, linetype = 2, color = "darkgrey", alpha = 0.8) +
        facet_wrap(~motif_id, ncol = 4, scales = "free_x") +
        xlab("chromVAR deviation score") +
        ylab("disease") +
        ggtitle(paste(selected_cluster, selected_region, sep = ", ")) +
        scale_fill_manual(values = disease_color_panel[c(selected_disease, "Control")]) +
        scale_color_manual(values = lightness(disease_color_panel[c(selected_disease, "Control")], scalefac(0.3))) +
        theme_bw(base_size = 8, base_family = "Helvetica") +
        theme(
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            legend.position = "none"
        )
}

plot_ridge_chromVAR_z_by_disease <- function(selected_motifs, selected_cluster, selected_region, selected_disease) {
    idx <- which(
        chromVAR_se$atac_anno_level_2 == selected_cluster &
            chromVAR_se$region == selected_region &
            chromVAR_se$disease %in% c(selected_disease, "Control")
    )

    df <- tibble(
        disease = factor(chromVAR_se$disease[idx], levels = c(selected_disease, "Control"))
    )

    inputs <- mat_z[selected_motifs, idx]

    motif_data_df <- inputs %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("cell_id") %>%
        as_tibble() %>%
        mutate(
            disease = df$disease,
            region = selected_region,
            cluster = selected_cluster
        ) %>%
        pivot_longer(cols = -one_of(c("cell_id", "disease", "cluster", "region")), names_to = "motif_id", values_to = "dev_score")

    p <- motif_data_df %>% ggplot(aes(dev_score, disease))
    p +
        geom_density_ridges(
            aes(fill = disease, color = disease),
            rel_min_height = 0.001,
            size = 0.3,
            scale = 1,
            quantile_lines = TRUE,
            quantiles = 2
        ) +
        # geom_vline(xintercept = 0, linetype = 2, color = "darkgrey", alpha = 0.8) +
        facet_wrap(~motif_id, ncol = 4, scales = "free_x") +
        xlab("chromVAR Z") +
        ylab("disease") +
        ggtitle(paste(selected_cluster, selected_region, sep = ", ")) +
        scale_fill_manual(values = disease_color_panel[c(selected_disease, "Control")]) +
        scale_color_manual(values = lightness(disease_color_panel[c(selected_disease, "Control")], scalefac(0.3))) +
        theme_bw(base_size = 8, base_family = "Helvetica") +
        theme(
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            legend.position = "none"
        )
}

plot_ridge_chromVAR_dev_by_disease(top_each_disease$motif_id[1:3], "Astro", "MCX", "ALS")
plot_ridge_chromVAR_z_by_disease(top_each_disease$motif_id[1:3], "Astro", "MCX", "ALS")


sig_Astro_motif <- merged_df_z_Astro %>%
    filter(sig_both == "yes") %>%
    distinct(archetype_name) %>%
    pull(archetype_name)
plot_ridge_chromVAR_z_by_disease(sig_Astro_motif, "Astro", "MCX", "ALS")
ggsave("./plots/chromVAR_diff_Z_vs_snRNA_expFC_Astro_MCX_ALS_vs_Control_sig_motif_ridges.pdf",
    device = cairo_pdf(), width = 6, height = 4.5, useDingbats = FALSE
)

sig_ExcSuperficial_motif <- merged_df_z_ExcSuperficial %>%
    filter(sig_both == "yes") %>%
    distinct(archetype_name) %>%
    pull(archetype_name)
plot_ridge_chromVAR_z_by_disease(sig_ExcSuperficial_motif, "Exc_superficial", "MCX", "ALS")
ggsave("./plots/chromVAR_diff_Z_vs_snRNA_expFC_ExcSuperficial_MCX_ALS_vs_Control_sig_motif_ridges.pdf",
    device = cairo_pdf(), width = 6, height = 4.5, useDingbats = FALSE
)

###
# testing between two groups, using the average of each individual as observations
# compute the matrix of average chromVAR deviations by group and individual
meta_data <- getCellColData(proj) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()

mat_deviations_col_order <- colnames(mat_deviations)

meta_data_ordered <- tibble(cell_id = mat_deviations_col_order) %>%
    left_join(meta_data, by = "cell_id") %>%
    mutate(group_subject = paste(group, subject, sep = "_"))

all.equal(meta_data_ordered$cell_id, colnames(mat_deviations))

group_subject <- meta_data_ordered$group_subject %>% as_factor()
names(group_subject) <- meta_data_ordered$cell_id

mm_subject <- model.matrix(~ 0 + group_subject)
colnames(mm_subject) <- levels(group_subject)

mat_deviations_subject_sum <- mat_deviations %*% mm_subject

group_subject_count <- meta_data_ordered %>% dplyr::count(group_subject)
count_vec2 <- group_subject_count$n
names(count_vec2) <- group_subject_count$group_subject
count_vec2 <- count_vec2[colnames(mat_deviations_subject_sum)]

all.equal(colnames(mat_deviations_subject_sum), names(count_vec2))

mat_deviations_subject_avg <- mat_deviations_subject_sum / count_vec2[col(mat_deviations_subject_sum)]

# disease (ALS or FTD) vs. control, separated by brain region
# deviations diff
differentialDeviations_disease_by_individual <- function(selected_cluster, selected_region, selected_disease) {
    group_id <- colnames(mat_deviations_subject_avg)
    idx <- grepl(selected_cluster, group_id) &
        grepl(selected_region, group_id) &
        (grepl(selected_disease, group_id) | grepl("Control", group_id)) &
        !(grepl("Mixed__", group_id))
    inputs <- mat_deviations_subject_avg[, idx]
    df <- tibble(group_name = colnames(inputs)) %>%
        mutate(
            disease = if_else(grepl("Control", group_name), "Control", selected_disease)
        )
    groups <- factor(df$disease, levels = c(selected_disease, "Control"))
    test_res <- apply(inputs, 1, t_helper_mod, groups)
    bind_rows(test_res) %>%
        mutate(
            motif_id = names(test_res),
            mean_diff = mean_x - mean_y,
            fdr = p.adjust(p, method = "fdr"),
            cluster = selected_cluster,
            disease = selected_disease,
            region = selected_region
        ) %>%
        select(cluster, region, disease, motif_id, mean_x, mean_y, mean_diff, p, fdr) %>%
        arrange(fdr)
}

params <- expand_grid(
    selected_cluster = c(
        "Astro", "Micro", "Oligo", "OPC",
        "Exc_superficial", "Exc_intermediate", "Exc_deep",
        "Inh_VIP", "Inh_LAMP5", "Inh_PVALB", "Inh_SST"
    ),
    selected_region = c("MCX", "mFCX"),
    selected_disease = c("ALS", "FTD")
)
diff_dev_disease_by_individual <- pmap_dfr(params, differentialDeviations_disease_by_individual)

write_tsv(diff_dev_disease_by_individual, "chromVAR_diff_deviations_by_disease_vierstra_motif_archetype_individualMeanAsObs.txt")
# diff_dev_disease_by_individual <- read_tsv("./chromVAR_diff_deviations_by_disease_vierstra_motif_archetype_individualMeanAsObs.txt")

cell_type_color_palette <- read_tsv("./color_palette_level_2.txt")
cluster_color <- cell_type_color_palette$color
names(cluster_color) <- cell_type_color_palette$sub_cluster

to_plot_diff_dev_disease_by_individual <- diff_dev_disease_by_individual %>%
    filter(!grepl("Mixed", cluster)) %>%
    mutate(cluster = factor(cluster, levels = cell_type_color_palette$sub_cluster)) %>%
    separate(motif_id, c("archetype_id", "tf_name", "family_name"), sep = "\\|") %>%
    mutate(
        tf_name = str_replace_all(tf_name, "/", ", "),
        # name_clean = str_glue("{tf_name} ({family_name})")
        name_clean = str_glue("{tf_name}")
    )

top_each_disease2 <- to_plot_diff_dev_disease_by_individual %>%
    filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05)) %>%
    arrange(-abs(mean_diff)) %>%
    group_by(cluster, region, disease) %>%
    dplyr::slice(1:3)

p <- to_plot_diff_dev_disease_by_individual %>%
    ggplot(aes(mean_diff, -log10(fdr)))
p +
    geom_point_rast(color = "grey50", size = 0.5, scale = 0.8) +
    geom_point_rast(aes(color = cluster),
        size = 0.5,
        data = to_plot_diff_dev_disease_by_individual %>%
            # filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05))
            filter(fdr < 0.05, abs(mean_diff) > 0.01),
        scale = 0.8
    ) +
    geom_point(
        aes(color = "black"),
        size = 0.5,
        data = top_each_disease2,
        shape = 1
    ) +
    geom_text_repel(
        aes(label = name_clean),
        data = top_each_disease2,
        nudge_x = 2,
        size = 1
    ) +
    facet_grid(
        disease + region ~ cluster,
        scales = "fixed",
        switch = "y"
    ) +
    xlab("Mean difference in accessibility\n(chromVAR deviations, disease vs. control)") +
    ylab("-log10(FDR)") +
    scale_color_manual(values = cluster_color) +
    xlim(-0.08, 0.08) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside"
    )
ggsave("./plots/chromVAR_diff_deviations_by_disease_vierstra_motif_archetype_top_TF_volcano_individualMeanAsObs.pdf",
    device = cairo_pdf(), width = 7, height = 3, useDingbats = FALSE
)

motif_diff_mat <- diff_dev_disease_by_individual %>%
    filter(!grepl("Mixed", cluster)) %>%
    filter(disease == "ALS") %>%
    filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05)) %>%
    # filter((mean_x > 0 & mean_diff > 0 & fdr < 0.01) | (mean_y > 0 & mean_diff < 0 & fdr < 0.01)) %>%
    # filter((mean_x > 0 & mean_diff > 0.02 & fdr < 0.05) | (mean_y > 0 & mean_diff < -0.02 & fdr < 0.05)) %>%
    arrange(-abs(mean_diff)) %>%
    mutate(group = paste(cluster, region, disease, sep = "__")) %>%
    select(group, motif_id, mean_diff) %>%
    pivot_wider(names_from = "motif_id", values_from = "mean_diff", values_fill = 0) %>%
    arrange(group) %>%
    as.data.frame() %>%
    column_to_rownames("group") %>%
    as.matrix() %>%
    t()

pheatmap(
    motif_diff_mat,
    color = scico(20, palette = "vik"),
    # color = scico(20, palette = "lajolla"),
    breaks = seq(from = -0.05, to = 0.05, length.out = 21),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "none",
    # show_rownames = TRUE
    show_rownames = FALSE
)

###
# plot differential results as ridge plot by individual
mat_deviations_imputed <- imputeMatrix(
    mat = as.matrix(mat_deviations),
    imputeWeights = impute_weight
)

# to_plot_motif_id <- var_dev$name[c(1, 2, 5)]
to_plot_motif_id <- c(
    "AC0242|FOSL/JUND|bZIP",
    "AC0239|BACH/NFE|bZIP",
    "AC0244|MAFG/MAFF|bZIP",
    "AC0078|CTCF/CTCFL|C2H2_ZF",
    "AC0038|THAP|['C2CH_THAP-type_zinc_finger_factors']",
    "AC0083|ZBTB|C2H2_ZF",
    "AC0505|ZNF|C2H2_ZF",
    "AC0042|ZSCAN|C2H2_ZF",
    "AC0626|ELF/ETV|Ets",
    "AC0525|NFYA/NFYB|CBF/NF-Y"
)

to_plot_motif_id <- c(
    "AC0015|RUNX|Runt",
    "AC0231|PAX|Homeodomain,Paired_box",
    "AC0243|ZSCAN5C|C2H2_ZF",
    "AC0265|FOXJ|Forkhead",
    "AC0269|CUX|CUT,Homeodomain",
    "AC0270|BBX|Sox",
    "AC0283|ZSCAN|C2H2_ZF,MADF",
    "AC0295|ZNF|C2H2_ZF",
    "AC0367|ZNF|C2H2_ZF",
    "AC0505|ZNF|C2H2_ZF",
    "AC0604|ZNF|C2H2_ZF"
)

to_plot_motif_id <- c(
    # "AC0376|ZBED|BED_ZF"
    "AC0398|POU6F|Homeodomain"
)

to_plot_motif_id <- c(
    "AC0012|SATB|['Homeo_domain_factors']",
    "AC0019|AIRE|SAND",
    "AC0033|PRDM|C2H2_ZF",
    "AC0042|ZSCAN|C2H2_ZF",
    "AC0045|ZNF|C2H2_ZF"
)

# to_plot_value <- mat_deviations_imputed[to_plot_motif_id, , drop = F] %>%
to_plot_value <- mat_deviations[to_plot_motif_id, , drop = F] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()
to_plot_meta <- getCellColData(proj) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    as_tibble()

to_plot_res <- to_plot_meta %>%
    left_join(to_plot_value, by = "cell_id") %>%
    select(cell_id, region, disease, subject, atac_anno_level_2, starts_with("AC")) %>%
    filter(!grepl("Mixed__", atac_anno_level_2)) %>%
    pivot_longer(starts_with("AC"), names_to = "motif", values_to = "chromVAR_dev") %>%
    mutate(
        atac_anno_level_2 = factor(atac_anno_level_2,
            levels = rev(c(
                "Exc_superficial", "Exc_intermediate", "Exc_deep",
                "Inh_LAMP5", "Inh_VIP", "Inh_PVALB", "Inh_SST",
                "Astro", "Micro", "Oligo", "OPC"
            ))
        ),
        motif = factor(motif, levels = to_plot_motif_id)
    )

p <- to_plot_res %>%
    filter(region == "MCX", disease != "FTD", atac_anno_level_2 == "Astro") %>%
    # filter(region == "MCX", disease != "FTD", atac_anno_level_2 == "Exc_deep") %>%
    arrange(disease) %>%
    mutate(subject = factor(subject,
        levels = c("945", "91", "906", "904", "902", "1069", "52", "388", "332", "113", "111", "110")
    )) %>%
    ggplot(aes(y = subject, x = chromVAR_dev))

p +
    # geom_density_ridges(
    #     aes(fill = disease),
    #     rel_min_height = 0.001,
    #     size = 0.3,
    #     scale = 2,
    # ) +
    stat_density_ridges(
        aes(fill = disease),
        rel_min_height = 0.001,
        size = 0.3,
        scale = 2,
        quantile_lines = TRUE,
        quantiles = 2
    ) +
    facet_wrap(~motif, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = c("ALS" = "red", "Control" = "grey")) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )

df_deviations_subject_avg <- mat_deviations_subject_avg %>%
    as.data.frame() %>%
    rownames_to_column("motif_id") %>%
    as_tibble() %>%
    pivot_longer(-motif_id, names_to = "group", values_to = "chromVAR_dev") %>%
    mutate(
        cell_type = str_extract(group, ".*(?=_(mFCX|MCX))"),
        region = if_else(grepl("_MCX_", group), "MCX", "mFCX"),
        disease = case_when(
            grepl("ALS", group) ~ "ALS",
            grepl("Control", group) ~ "Control",
            grepl("FTD", group) ~ "FTD",
        ),
        subject = str_extract(group, "(?<=_(Control|ALS|FTD)_).*")
    ) %>%
    filter(!grepl("Mixed", cell_type))

p <- df_deviations_subject_avg %>%
    filter(
        motif_id %in% to_plot_motif_id,
        # cell_type == "Astro",
        cell_type == "Oligo",
        # region == "MCX",
        region == "mFCX",
        # disease %in% c("ALS", "Control")
        disease %in% c("FTD", "Control")
    ) %>%
    mutate(motif_id = factor(motif_id, levels = to_plot_motif_id)) %>%
    ggplot(aes(disease, chromVAR_dev))

p +
    geom_boxplot(aes(color = disease), width = 0.4, outlier.shape = NA) +
    geom_beeswarm(aes(color = disease)) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    facet_wrap(~motif_id, scales = "free_y", nrow = 2)

## Astro, diff results using individuals as observations
p <- diff_dev_disease_by_individual %>%
    filter(cluster == "Astro", region == "MCX", disease == "ALS") %>%
    ggplot(aes(mean_diff, -log10(fdr)))
p +
    geom_point(color = "grey50", size = 0.5) +
    geom_point(aes(color = cluster),
        size = 0.5,
        data = diff_dev_disease_by_individual %>%
            filter(cluster == "Astro", region == "MCX", disease == "ALS") %>%
            filter((mean_x > 0 & mean_diff > 0 & fdr < 0.05) | (mean_y > 0 & mean_diff < 0 & fdr < 0.05))
    ) +
    geom_text_repel(
        aes(label = motif_id),
        data = top_each_disease2 %>%
            filter(cluster == "Astro", region == "MCX", disease == "ALS"),
        size = 1.8
    ) +
    facet_wrap(~cluster, scales = "free_x") +
    xlab("Mean difference in accessibility\n(chromVAR deviations, disease vs. control)") +
    ylab("-log10(FDR)") +
    # scale_color_manual(values = cluster_color) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )

diff_dev_disease_sub_by_individual <- diff_dev_disease_by_individual %>%
    filter(
        cluster == "Astro",
        region == "MCX",
        disease == "ALS"
    ) %>%
    mutate(archetype_id = str_extract(motif_id, "AC\\d+")) %>%
    dplyr::rename(archetype_name = motif_id) %>%
    select(archetype_id, archetype_name, fdr, mean_x, mean_y, mean_diff)

merged_df_Astro_by_individual <- motif_TF_lookup %>%
    left_join(diff_dev_disease_sub_by_individual, by = c("cluster" = "archetype_id")) %>%
    left_join(snRNA_DE_sub, by = c("tf_name" = "gene")) %>%
    drop_na() %>%
    mutate(sig_both = if_else(snRNA_sig == "yes" & fdr < 0.05, "yes", "no"))

p <- merged_df_Astro_by_individual %>%
    ggplot(aes(model_log2FC, mean_diff))

p +
    geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
    geom_hline(yintercept = 0, color = "darkgrey", linetype = 2) +
    geom_point(color = "grey", size = 0.5, data = merged_df_Astro_by_individual %>% filter(sig_both == "no")) +
    geom_point(color = "#CD2526", size = 1, data = merged_df_Astro_by_individual %>% filter(sig_both == "yes")) +
    geom_text_repel(
        aes(label = tf_name),
        data = merged_df_Astro_by_individual %>%
            filter(sig_both == "yes") %>%
            select(tf_name, model_log2FC, mean_diff) %>%
            distinct(),
        size = 1.8
    ) +
    xlab("snRNA TF expresison log2FC\n(ALS vs. Control)") +
    ylab("snATAC Motif archetype chromVAR deviation difference\n(ALS vs. Control)") +
    ggtitle("Astro, motor cortex") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())

###

# session info
sessionInfo()