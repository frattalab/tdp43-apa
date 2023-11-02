library(tidyverse)

dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-02_papa_cryptics_kmer6_window_200_distal_window_250.cleaned_6mer_distribution_genome.tsv")
simp_dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-02_papa_cryptics_kmer6_window_200_distal_window_250.cleaned_6mer_distribution_genome_simple.tsv")


# adjust p-values across all comparisons
simp_dbrn_tbl <- simp_dbrn_tbl %>%
  mutate(padj = p.adjust(pvalue, method = "BH"))


# Rank each kmer within experiment, contrast, group by descending PEKA-score
kmer_ranks <- simp_dbrn_tbl %>%
  group_by(comparison_name) %>%
  arrange(desc(PEKA_score), .by_group = T) %>%
  mutate(peka_rank = row_number(),
         peka_rank_frac = peka_rank / n()) %>%
  ungroup() %>%
  select(kmer, comparison_name, PEKA_score, padj, peka_rank, peka_rank_frac)



# kmer peka score ranks overall across comparisons (irrespective of direction)
ave_all_kmer_ranks <- kmer_ranks %>%
  group_by(kmer) %>%
  summarise(ranksum = sum(peka_rank),
            mean_rank = mean(peka_rank),
            median_rank = median(peka_rank)
  ) %>%
  arrange(median_rank)

ave_all_kmer_ranks

# if stratify by relative genomic position, what are the most consistently enriched kmers
kmer_ranks %>%
  mutate(rel_pos = if_else(str_ends(comparison_name, "exonstart|proximal"), "proximal", "distal")) %>%
  group_by(kmer, rel_pos) %>%
  summarise(ranksum = sum(peka_rank),
            mean_rank = mean(peka_rank),
            median_rank = median(peka_rank)
  ) %>%
  arrange(desc(rel_pos), median_rank)

# if just look at splice sites, what are the most consistently enriched kmers
kmer_ranks %>%
  filter(str_ends(comparison_name, "exonstart")) %>%
  group_by(kmer) %>%
  summarise(ranksum = sum(peka_rank),
            mean_rank = mean(peka_rank),
            median_rank = median(peka_rank)
  ) %>%
  arrange(median_rank)
  

# if just look at pas of internal, what are the most consistently enriched kmers
kmer_ranks %>%
  filter(str_ends(comparison_name, "pas")) %>%
  group_by(kmer) %>%
  summarise(ranksum = sum(peka_rank),
            mean_rank = mean(peka_rank),
            median_rank = median(peka_rank)
  ) %>%
  arrange(median_rank)


# plot median ranks of top n motifs across all experiments
ave_all_kmer_ranks %>%
  slice_min(median_rank, n = 20) %>%
  ungroup() %>%
  mutate(kmer = fct_reorder(kmer, median_rank)) %>%
  ggplot(aes(x = kmer, y = median_rank)) +
  geom_col() +
  scale_y_continuous(breaks = seq(0,500,100),
                     minor_breaks = seq(0,500,50)) +
  theme_bw(base_size = 20) +
  labs(title = "Top 20 most enriched 6mers across comparisons",
       x = "kmer",
       y = "Median rank of PEKA score") +
  theme(axis.text.x = element_text(angle = 90))

# plot ranks of top 20 in each comparison
kmer_ranks %>%
  group_by(comparison_name) %>%
  slice_min(peka_rank, n = 25) %>%
  ungroup() %>%
  mutate(comparison_name = factor(comparison_name, levels = names(dbrn_paths)),
         kmer = tidytext::reorder_within(kmer, peka_rank, comparison_name)) %>%
  ggplot(aes(x = kmer, y = peka_rank)) +
  facet_wrap("comparison_name", scales = "free", ncol = 2) +
  geom_col() +
  scale_x_reordered() +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90))
  

