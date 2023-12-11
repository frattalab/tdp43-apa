library(tidyverse)
source("scripts/fncs_plot_peka.R")


dbrn_tbl <- read_tsv("processed/peka/papa/2023-11_27_papa_cryptics_kmer6_window_250_distal_window_500_relpos_0.cleaned_6mer_distribution_genome.tsv")

# map comparison names to cleaned event types & region types
plot_clean_names <- tibble(comparison_name = c("bleedthrough_exonstart", "bleedthrough_pas", "spliced_exonstart", "spliced_pas", "d3utr_pas_proximal", "d3utr_pas_distal"),
                           plot_event_type = c("Bleedthrough-ALE", "Bleedthrough-ALE", "AS-ALE", "AS-ALE", "3'UTR-ALE", "3'UTR-ALE"),
                           plot_region_type = factor(c("Exon Start", "PAS", "Exon Start", "PAS", "Proximal", "Distal"),
                                                     levels = c("Exon Start", "PAS", "Proximal", "Distal")
                           )
)

# Set the A4 paper size dimensions in inches
a4_width <- 8.27
a4_height <- 11.69

## preprocess PEKA dbrn tables for ploting

dbrn_tbl_long_gu <- peka_wide_to_long(dbrn_tbl, c("GUGUGU"), first_posn_idx = 10) %>% 
  separate(comparison_name, into = c("comparison_name", "background_type"), "\\.") %>%
  left_join(plot_clean_names, by = "comparison_name") %>%
  mutate(plot_background_type = if_else(background_type == "background_all", "Any dataset", "SH-SY5Y only"))

dbrn_tbl_long_ug <- peka_wide_to_long(dbrn_tbl, c("UGUGUG"), first_posn_idx = 10) %>% 
  separate(comparison_name, into = c("comparison_name", "background_type"), "\\.") %>%
  left_join(plot_clean_names, by = "comparison_name") %>%
  mutate(plot_background_type = if_else(background_type == "background_all", "Any dataset", "SH-SY5Y only"))

dbrn_tbl_long_both <- peka_wide_to_long(dbrn_tbl, c("UGUGUG", "GUGUGU"),
                                        first_posn_idx = 10, sum_occur = T, 
                                        sum_group_cols = c( "comparison_name", "rel_posn")) %>% 
  separate(comparison_name, into = c("comparison_name", "background_type"), "\\.") %>%
  left_join(plot_clean_names, by = "comparison_name") %>%
  mutate(plot_background_type = if_else(background_type == "background_all", "Any dataset", "SH-SY5Y only"))

dbrn_tbl_long_both_nosum <- peka_wide_to_long(dbrn_tbl, c("UGUGUG", "GUGUGU"), first_posn_idx = 10, sum_occur = F) %>% 
  separate(comparison_name, into = c("comparison_name", "background_type"), "\\.") %>%
  left_join(plot_clean_names, by = "comparison_name") %>%
  mutate(plot_background_type = if_else(background_type == "background_all", "Any dataset", "SH-SY5Y only"))


# plot just GU 6mer
all_map_gu <- plot_kmer_dbrn(mutate(dbrn_tbl_long_gu, direction = plot_background_type),
               rolling_mean = T,
               facet_w = "plot_event_type ~ plot_region_type",
               n_col = 2,
               n_row = 3,
               breaks = seq(-500,500,100),
               rolling_k = 10
               ) +
  labs(title = "GUGUGU only",
       colour = "Background")

# plot just UG 6mer
all_map_ug <- plot_kmer_dbrn(mutate(dbrn_tbl_long_ug, direction = background_type),
               rolling_mean = T,
               facet_w = "plot_event_type ~ plot_region_type",
               n_col = 2,
               n_row = 3,
               breaks = seq(-500,500,100),
               rolling_k = 10
) +
  labs(title = "UGUGUG only",
       colour = "Background")

# plot GU and UG dbrns with rel exprn combined
all_map_both <- plot_kmer_dbrn(mutate(dbrn_tbl_long_both, direction = background_type),
               rolling_mean = T,
               facet_w = "plot_event_type ~ plot_region_type",
               n_col = 2,
               n_row = 3,
               breaks = seq(-500,500,100),
               rolling_k = 10
) +
  labs(title = "GUGUGU & UGUGUG summed",
       colour = "Background")

# # plot GU and UG dbrns as separate traces
# all_map_both_nosum <- plot_kmer_dbrn(mutate(dbrn_tbl_long_both_nosum, direction = kmer),
#                rolling_mean = T,
#                facet_w = "plot_event_type ~ plot_region_type",
#                n_col = 2,
#                n_row = 3,
#                breaks = seq(-500,500,100),
#                rolling_k = 10
#                ) +
#   labs(title = "GUGUGU & UGUGUG",
#        direction = "kmer")
               



### Since somewhat arbitrarily picked out the GU & UG motifs, let's check PEKA's top n 'enriched' kmers and see if they have positional enrichment resembling binding profile by iCLIP


# Extract top 5-10 kmers for each comparison
tmp_top10_grpd <- dbrn_tbl %>%
  select(kmer, comparison_name, PEKA_score) %>%
  group_by(comparison_name) %>%
  arrange(desc(PEKA_score), .by_group = T) %>%
  slice_max(PEKA_score, n = 10)

# get a named list of vectors of top 10 enriched kmers for each comparioson
 comparison_top10_kmer_list <- tmp_top10_grpd %>%
  group_split() %>%
  set_names(pull(group_keys(tmp_top10_grpd))) %>%
  map(~ pull(.x, kmer))
 
# for each comparison, get a plot ready df

dbrn_tbl_top10_comparison_list <- map2(names(comparison_top10_kmer_list),
     comparison_top10_kmer_list,
     ~ dbrn_tbl %>%
       filter(comparison_name == .x) %>%
       peka_wide_to_long(kmers = .y, first_posn_idx = 10, sum_occur = FALSE) %>%
       separate(comparison_name, into = c("comparison_name", "background_type"), "\\.") %>%
       left_join(plot_clean_names, by = "comparison_name") %>%
       mutate(direction = kmer)
     ) %>%
   set_names(names(comparison_top10_kmer_list))

# make plots for top 10 kmers 
plots_top10_comparison <- map2(.x = dbrn_tbl_top10_comparison_list,
                               .y = names(dbrn_tbl_top10_comparison_list),
    ~ plot_kmer_dbrn(.x,
                    rolling_mean = T,
                    facet_w = "plot_event_type ~ plot_region_type",
                    n_col = 2,
                    n_row = 3,
                    breaks = seq(-500,500,100),
                    rolling_k = 10,
                    plot_colours = str_split_1("#a6cee3\n#1f78b4\n#b2df8a\n#33a02c\n#fb9a99\n#e31a1c\n#fdbf6f\n#ff7f00\n#cab2d6\n#6a3d9a", "\\n")
                    ) +
      labs(title = str_split_i(.y, "\\.", 2),
           colour = "kmer")
    )
 



# Create a multi-page PDF file for each group of plots
if (!dir.exists("processed/peka/papa/plots/")) { dir.create("processed/peka/papa/plots/", recursive = T)}

pdf("processed/peka/papa/plots/2023-11-27_papa_all_comparisons_gu_ug_map.pdf", width = a4_height, height = a4_width )

# Loop through each ggplot object and print it to the PDF file
plots_ug_list <- list(all_map_ug, all_map_gu, all_map_both)
for (i in seq_along(plots_ug_list)) {
  print(plots_ug_list[[i]])
}

# Close the PDF file
dev.off()

# top 10 kmers plots
pdf("processed/peka/papa/plots/2023-11-27_papa_peka_top10_kmers_by_comparison_map.pdf", width = a4_height, height = a4_width )

# Loop through each ggplot object and print it to the PDF file
for (i in seq_along(plots_top10_comparison)) {
  print(plots_top10_comparison[[i]])
}

# Close the PDF file
dev.off()
