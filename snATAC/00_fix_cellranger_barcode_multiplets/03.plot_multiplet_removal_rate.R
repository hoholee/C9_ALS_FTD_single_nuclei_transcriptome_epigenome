library(tidyverse)

df <- read_tsv("count_cell_barcodes_pre_post_multiplets_removal.txt", col_names = c("sample", "original_count", "filtered_count")) %>% 
  mutate(rate = 100*(1 - (filtered_count/original_count)),
         rate = round(rate, 2))

p <- df %>% ggplot(aes(sample, original_count))

p + geom_bar(fill = "#DE315C", stat = "identity") +
  geom_bar(fill = "#629CFF", stat = "identity", data = df, aes(sample, filtered_count)) +
  geom_text(aes(sample, original_count + 200, label = paste0(rate, "%"))) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  ylab("# Cells") +
  coord_flip()
  
