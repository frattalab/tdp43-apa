library(tidyverse)
library(tidytext)

# read in kmer distribution tables for all experiments
dbrn_paths <- list.files(path = "data/peka_qapa",
                         pattern = "_6mer_distribution_genome.tsv$",
                         recursive = T,
                         full.names = T) %>%
  set_names(str_remove(basename(.), "_6mer_distribution_genome.tsv$"))

l_dbrn_tbls <- map(dbrn_paths,
                   ~ read_tsv(.x, show_col_types = F) %>% rename(kmer = `...1`)
)

#combine into single table, 
dbrn_tbl <- bind_rows(l_dbrn_tbls, .id = "experiment_comparison_name")
# split experiment_comparison_name into sub strings for easier categorisation
dbrn_tbl <- separate(dbrn_tbl,
                     experiment_comparison_name,
                     into = c("experiment_name", "contrast", "pas", NA, "direction"),
                     sep = "\\.",
                     remove = F
)

simp_dbrn_tbl <- dbrn_tbl %>%
  select(kmer,
         experiment_comparison_name,
         experiment_name,
         contrast,
         pas,
         direction,
         kmer,
         mtxn,
         artxn,
         aroxn,
         etxn,
         PEKA_score = `PEKA-score`,
         pvalue = `p-value`)

# Rank each kmer within experiment, contrast, group by descending PEKA-score
kmer_ranks <- simp_dbrn_tbl %>%
  group_by(experiment_comparison_name) %>%
  arrange(desc(PEKA_score), .by_group = T) %>%
  mutate(peka_rank = row_number(),
         peka_rank_frac = peka_rank / n()) %>%
  ungroup() %>%
  select(kmer, experiment_comparison_name, experiment_name, contrast, pas, direction, peka_rank, peka_rank_frac)



# kmer peka score ranks overall across comparisons (irrespective of direction)
ave_all_kmer_ranks <- kmer_ranks %>%
  group_by(kmer) %>%
  summarise(ranksum = sum(peka_rank),
            mean_rank = mean(peka_rank),
            median_rank = median(peka_rank)
  ) %>%
  arrange(mean_rank)

ave_all_kmer_ranks

# kmer ranks by site type and direction
ave_split_kmer_ranks <- kmer_ranks %>%
  group_by(kmer, pas, direction) %>%
  summarise(ranksum = sum(peka_rank),
            mean_rank = mean(peka_rank),
            median_rank = median(peka_rank)
            ) %>%
  ungroup() %>%
  arrange(mean_rank)

ave_split_kmer_ranks


# plot median ranks of top for each comparison
median_top20_bar <- ave_split_kmer_ranks %>%
  group_by(pas, direction) %>%
  slice_min(mean_rank, n = 20) %>%
  ungroup() %>%
  mutate(kmer = tidytext::reorder_within(kmer, median_rank ,list(pas, direction))) %>%
  mutate(pas = factor(pas, levels = c("proximal", "distal")),
         direction = factor(direction, levels = c("up", "down"))
         ) %>%
  ggplot(aes(x = kmer, y = median_rank)) +
  facet_wrap("pas ~ direction", scales = "free") +
  geom_col() +
  scale_x_reordered() + 
  scale_y_continuous(breaks = seq(0,500,100),
                     minor_breaks = seq(0,500,50)) +
  theme_bw(base_size = 14) +
  labs(title = "Top 20 most enriched 6mers around regulated PAS across datasets & comparisons",
       x = "kmer",
       y = "Median rank of PEKA score") +
  theme(axis.text.x = element_text(angle = 90))
    

median_top20_bar

ggsave(plot = median_top20_bar,
       filename = "2023-06-06_median_peka_score_split_top20_bar.png",
       path = "processed",
       device = "png",
       units = "in",
       height = 10,
       width = 12,
       dpi = "retina")


