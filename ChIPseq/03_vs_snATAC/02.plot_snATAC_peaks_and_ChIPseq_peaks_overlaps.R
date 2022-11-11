library(tidyverse)

df <- read_tsv("snATAC_peaks_overlaps_with_ChIPseq_consensus_peaks.txt") %>%
	mutate(num_non_overlaps = num_snATAC_merged_peaks - num_overlaps_with_ChIPseq_consensus_peaks) %>%
	select(-num_snATAC_merged_peaks) %>%
	pivot_longer(-cell_type, names_to = "group", values_to = "count")

p <- df %>% ggplot(aes(cell_type, count))
p + geom_bar(aes(fill=group), stat = "identity", position = position_stack()) +
	theme_bw(base_size = 10, base_family = "Helvetica")
ggsave("snATAC_peaks_overlaps.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE)

df <- read_tsv("ChIPseq_consensus_peaks_overlaps_with_snATAC_peaks.txt") %>%
	mutate(num_non_overlaps = num_ChIPseq_consensus_peaks - num_overlaps_with_snATAC_merged_peaks) %>%
	select(-num_ChIPseq_consensus_peaks) %>%
	pivot_longer(-cell_type, names_to = "group", values_to = "count")

p <- df %>% ggplot(aes(cell_type, count))
p + geom_bar(aes(fill=group), stat = "identity", position = position_stack()) +
	theme_bw(base_size = 10, base_family = "Helvetica")
ggsave("ChIPseq_peaks_overlaps.pdf", device = cairo_pdf(), width = 6, height = 4, useDingbats = FALSE)
