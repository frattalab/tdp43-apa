library(tidyverse)
library(tidytext)

dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-03_papa_cryptics_kmer6_window_250_distal_window_500.cleaned_6mer_distribution_genome.tsv")
simp_dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-03_papa_cryptics_kmer6_window_250_distal_window_500.cleaned_6mer_distribution_genome_simple.tsv")

# map comparison names to cleaned event types & region types
plot_clean_names <- tibble(comparison_name = c("bleedthrough_exonstart", "bleedthrough_pas", "spliced_exonstart", "spliced_pas", "d3utr_pas_proximal", "d3utr_pas_distal"),
                           plot_event_type = c("Bleedthrough-ALE", "Bleedthrough-ALE", "AS-ALE", "AS-ALE", "3'UTR-ALE", "3'UTR-ALE"),
                           plot_region_type = factor(c("Exon Start", "PAS", "Exon Start", "PAS", "Proximal", "Distal"),
                                                     levels = c("Exon Start", "PAS", "Proximal", "Distal")
                           )
)

simp_dbrn_tbl <- simp_dbrn_tbl %>%
  left_join(plot_clean_names,by = "comparison_name") %>%
  # construct a cleaned combined name
  mutate(plot_comparison_name = paste(plot_event_type, plot_region_type, sep = " "))



comparison_names <- unique(simp_dbrn_tbl$comparison_name)
plot_comparison_names <- unique(simp_dbrn_tbl$plot_comparison_name)

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
  select(kmer, comparison_name, plot_comparison_name, plot_event_type, PEKA_score, padj, peka_rank, peka_rank_frac)

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

# plot PEKA score distributions across comparisons
# Essentially a Z-score, so good to check if looks normally distributed and pvalues are valid
simp_dbrn_tbl %>%
  ggplot(aes(x = PEKA_score)) +
  facet_wrap("~ plot_comparison_name") +
  geom_density() +
  theme_bw() +
  labs(title = "relpos = auto")

# They all have peak around 0, but distribution has a really long extreme right-tail

simp_dbrn_tbl %>%
  ggplot(aes(x = pvalue)) +
  facet_wrap("~ plot_comparison_name", scales = "free_y") +
  geom_histogram(bins = 200) +
  theme_bw()

# some peak at low p-values, but generally have a weird peak a ~ 0.75

comparison_colours <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
# plot median ranks of top n motifs across all experiments
median_all_kmer_ranks_bar <- ave_all_kmer_ranks %>%
  slice_min(median_rank, n = 25) %>%
  ungroup() %>%
  mutate(kmer = fct_reorder(kmer, median_rank)) %>%
  ggplot(aes(x = kmer, y = median_rank)) +
  geom_col() +
  scale_y_continuous(breaks = seq(0,500,100),
                     minor_breaks = seq(0,500,50)) +
  theme_bw(base_size = 20) +
  labs(title = "Top 25 most enriched 6mers across comparisons",
       subtitle = "relpos = auto",
       x = "kmer",
       y = "Median rank of PEKA score") +
  theme(axis.text.x = element_text(angle = 90))

median_all_kmer_ranks_bar

median_all_kmer_ranks_box_colour <- ave_all_kmer_ranks %>%
  # pick top 25 enriched by median rank across comparisons ('generally enriched motifs')
  slice_min(median_rank, n = 25) %>%
  left_join(kmer_ranks, by = "kmer") %>%
  # sort axis from min ranks to max ranks
  mutate(kmer = fct_reorder(kmer, median_rank)) %>%
  ggplot(aes(x = kmer, y = log2(peka_rank))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(x = kmer, y = log2(peka_rank), colour = plot_comparison_name, shape = plot_comparison_name),
              width = 0.3) +
  scale_colour_manual(name = "",
                      values = comparison_colours) +
  theme_bw(base_size = 20) +
  labs(title = "Top 25 most enriched 6mers across comparisons",
       x = "kmer",
       y = "log2 Median rank of PEKA score",
       shape = ""
       ) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")

median_all_kmer_ranks_box_colour

# plot ranks of top 25 in each comparison
top_comp_kmer_ranks_bar <- kmer_ranks %>%
  group_by(plot_comparison_name) %>%
  slice_min(peka_rank, n = 25) %>%
  ungroup() %>%
  mutate(#plot_comparison_name = factor(plot_comparison_name, levels = plot_comparison_names),
         kmer = tidytext::reorder_within(kmer, peka_rank, plot_comparison_name)) %>%
  ggplot(aes(x = kmer, y = peka_rank)) +
  facet_wrap("plot_comparison_name", scales = "free", ncol = 2) +
  geom_col() +
  scale_x_reordered() +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "",
       y = "Rank of PEKA score",
       )

top_comp_kmer_ranks_bar

# save plots
if (!dir.exists("processed/peka/papa/plots/")) {dir.create("processed/peka/papa/plots", recursive = T)}

ggsave(filename = "2023-11-16_median_peka_score_rank_bar_top25.png",
       plot = median_all_kmer_ranks_bar,
       device = "png",
       path = "processed/peka/papa/plots/",
       width = 10,
       height = 10,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-11-16_median_peka_score_rank_bar_top25.svg",
       plot = median_all_kmer_ranks_bar,
       device = svg,
       path = "processed/peka/papa/plots/",
       width = 10,
       height = 10,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-11-16_median_peka_score_rank_box_colourcomp_top25.png",
       plot = median_all_kmer_ranks_box_colour,
       device = "png",
       path = "processed/peka/papa/plots/",
       width = 10,
       height = 10,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-11-16_median_peka_score_rank_box_colourcomp_top25.svg",
       plot = median_all_kmer_ranks_box_colour,
       device = svg,
       path = "processed/peka/papa/plots/",
       width = 10,
       height = 10,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-11-16_percomparison_peka_score_rank_bar_top25.png",
       plot = top_comp_kmer_ranks_bar,
       device = "png",
       path = "processed/peka/papa/plots/",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-11-10_percomparison_peka_score_rank_bar_top25.svg",
       plot = top_comp_kmer_ranks_bar,
       device = svg,
       path = "processed/peka/papa/plots/",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")
