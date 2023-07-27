# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)

# difine samples to check
# selected_region <- "MCX"
selected_region <- "mFCX"
selected_disease_1 <- "ALS"
selected_disease_2 <- "Control"

# read observed scDist results
df <- read_tsv(
  str_glue(
    "./scDist_results_SCT_scaled_",
    "{selected_region}_",
    "{selected_disease_1}_vs_{selected_disease_2}.tsv"
  )
)

# add shuffle results
shuffle_files <- list.files(
  "./scDist_shuffle/",
  pattern = str_glue(".*{selected_region}_{selected_disease_1}_vs_{selected_disease_2}_permuted.*tsv"),
  full.names = TRUE
)

read_shuffle_files <- function(file) {
  read_tsv(file, col_types = "cdddcccc")
}

shuffle_df <- map_dfr(shuffle_files, read_shuffle_files)

# define summary functions
sem <- function(x) {
  N <- sum(!is.na(x))
  sd(x, na.rm = T) / sqrt(N)
}

ci <- function(x, conf.interval = 0.95) {
  N <- sum(!is.na(x))
  sem <- sd(x, na.rm = T) / sqrt(N)
  qt(1 - (1 - conf.interval) / 2, N - 1) * sem
}

shuffle_summary <- shuffle_df %>%
  group_by(cell_type) %>%
  summarize(
    dist_mean = mean(Dist.),
    dist_sd = sd(Dist.),
    dist_median = median(Dist.),
    dist_sem = sem(Dist.),
    dist_ci = ci(Dist.)
  ) %>%
  rename(Dist. = dist_mean)


# make plot
p <- df %>%
  ggplot(aes(cell_type, Dist.))
p +
  geom_point(stat = "identity", color = "red") +
  geom_point(aes(cell_type, Dist.), data = shuffle_summary, color = "darkgrey") +
  # geom_errorbar(aes(cell_type, ymin = Dist. - dist_sem, ymax = Dist. + dist_sem),
  #   data = shuffle_summary, width = 0.2, color = "darkgrey"
  # ) +
  # geom_errorbar(aes(cell_type, ymin = Dist. - dist_sd, ymax = Dist. + dist_sd),
  #   data = shuffle_summary, width = 0.2, color = "darkgrey"
  # ) +
  geom_errorbar(aes(cell_type, ymin = Dist. - dist_ci, ymax = Dist. + dist_ci),
    data = shuffle_summary, width = 0.2, color = "darkgrey"
  ) +
  xlab("Cell type level 2") +
  ylab(
    str_glue(
      "{selected_disease_1} vs {selected_disease_2} ",
      "differences in {selected_region} (scDist distance)"
    )
  ) +
  coord_flip() +
  theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# save plot
ggsave(
  str_glue(
    "./scDist_results_SCT_scaled_",
    "{selected_region}_",
    "{selected_disease_1}_vs_{selected_disease_2}_add_shuffles.pdf"
  ),
  device = cairo_pdf(),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
dev.off()


## log sessionInfo
sessionInfo()
