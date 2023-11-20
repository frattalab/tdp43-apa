library(tidyverse)
source("scripts/fncs_plot_peka.R")

#' Save a list of ggplots to pdf file, with one plot per page
plot_list_to_pdf <- function(plot_list, path, width, height) {
  
  pdf(path, width = width, height = height)
  
  # Loop through each ggplot object and print it to the PDF file
  for (i in seq_along(plot_list)) {
    print(plot_list[[i]])
  }
  # Close the PDF file
  dev.off()
  
}

dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-16_papa_cryptics_cvcoverage_window_500.gugugu_ugugug.distribution_genome.tsv")

unique(dbrn_tbl$comparison_id)
# map comparison names to cleaned event types & region types
plot_clean_names <- tibble(comparison_name = c("bleedthrough_exonstart", "bleedthrough_pas", "spliced_exonstart", "spliced_pas", "d3utr_pas_proximal", "d3utr_pas_distal"),
                           plot_event_type = c("Bleedthrough-ALE", "Bleedthrough-ALE", "AS-ALE", "AS-ALE", "3'UTR-ALE", "3'UTR-ALE"),
                           plot_region_type = factor(c("Exon Start", "PAS", "Exon Start", "PAS", "Proximal", "Distal"),
                                                     levels = c("Exon Start", "PAS", "Proximal", "Distal")
                                                     )
                           )

# add plot labels, convert to % age
dbrn_tbl <- dbrn_tbl %>%
  left_join(plot_clean_names, by = "comparison_name") %>%
  mutate(plot_group = if_else(group == "foreground", "Cryptic", "Background"),
         rel_occur = rel_occur * 100)

# map with free scales across y
map_freey <- plot_kmer_dbrn(mutate(dbrn_tbl, direction = plot_group),
               rolling_mean = T,
               facet_w = "plot_event_type ~ plot_region_type",
               n_col = 2,
               n_row = 3,
               breaks = seq(-500,500,100),
               rolling_k = 10,
               plot_colours = c("#000000", "#d95f02")
               ) + 
  labs(y = "% Coverage",
       colour = "") +
  theme(legend.position = "top")

# map with fixed scales across events
map_fixed <- plot_kmer_dbrn(mutate(dbrn_tbl, direction = plot_group),
               rolling_mean = T,
               facet_w = "plot_event_type ~ plot_region_type",
               n_col = 2,
               n_row = 3,
               breaks = seq(-500,500,100),
               rolling_k = 10,
               plot_colours = c("#000000", "#d95f02"),
               facet_scales = "fixed") + 
  labs(y = "% Coverage",
       colour = "") +
  theme(legend.position = "top")

map_freey
map_fixed

# Set the A4 paper size dimensions in inches
a4_width <- 8.27
a4_height <- 11.69

# save each plot as 1 page per a4
# Create a multi-page PDF file for each group of plots
if (!dir.exists("processed/peka/papa/plots/")) { dir.create("processed/peka/papa/plots/", recursive = T)}

plots_ug_list <- list(map_freey, map_fixed)
pdf("processed/peka/papa/plots/2023-11-16_papa_all_comparisons_cvcoverage_gu_ug_map.pdf",
    width = a4_height, height = a4_width )

# Loop through each ggplot object and print it to the PDF file

for (i in seq_along(plots_ug_list)) {
  print(plots_ug_list[[i]])
}

# Close the PDF file
dev.off()

#  repeat for Halleger et al motif groups
yg_dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-19_papa_cryptics_cvcoverage_window_500.yg_6mers.per_kmer_distribution_genome.tsv")
ya_dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-19_papa_cryptics_cvcoverage_window_500.ya_6mers.per_kmer_distribution_genome.tsv")
yg_6mers <- c("UGUGUG", "GUGUGU","UGUGCG", "UGCGUG","CGUGUG","GUGUGC")
ya_6mers <- c("AUGUGU", "GUAUGU", "GUGUAU", "UGUGUA", "UGUAUG", "UGCAUG")

# Sum occurrence of kmers across the yg groups
yg_sum_dbrn_tbl <-  yg_dbrn_tbl %>%
  group_by(region_type, group, rel_posn) %>%
  summarise(rel_occur = sum(rel_occur)) %>%
  ungroup()

ya_sum_dbrn_tbl <-  ya_dbrn_tbl %>%
  group_by(region_type, group, rel_posn) %>%
  summarise(rel_occur = sum(rel_occur)) %>%
  ungroup()

# prepare for ploitting - column names, cryptic/not etc. 

yg_sum_dbrn_tbl <- yg_sum_dbrn_tbl %>%
  rename(comparison_name = region_type) %>%
  left_join(plot_clean_names, by = "comparison_name") %>%
  mutate(plot_group = if_else(group == "foreground", "Cryptic", "Background"),
         rel_occur = rel_occur * 100)

ya_sum_dbrn_tbl <- ya_sum_dbrn_tbl %>%
  rename(comparison_name = region_type) %>%
  left_join(plot_clean_names, by = "comparison_name") %>%
  mutate(plot_group = if_else(group == "foreground", "Cryptic", "Background"),
         rel_occur = rel_occur * 100)

# maps with free scales across y
yg_ya_sum_map_freey <- map2(.x = list(YG_group = yg_sum_dbrn_tbl, YA_group = ya_sum_dbrn_tbl),
                            .y = list(YG_group = yg_6mers, YA_group = ya_6mers),
                            ~ plot_kmer_dbrn(mutate(.x, direction = plot_group),
                                             rolling_mean = T,
                                             facet_w = "plot_event_type ~ plot_region_type",
                                             n_col = 2,
                                             n_row = 3,
                                             breaks = seq(-500,500,100),
                                             rolling_k = 10,
                                             plot_colours = c("#000000", "#d95f02")
                            ) + 
                              labs(title = paste(.y, collapse = ","),
                                   y = "% Coverage",
                                   colour = "") +
                              theme(legend.position = "top")
)

# maps with fixed scales across y
yg_ya_sum_map_fixed <- map2(.x = list(YG_group = yg_sum_dbrn_tbl, YA_group = ya_sum_dbrn_tbl),
                            .y = list(YG_group = yg_6mers, YA_group = ya_6mers),
                            ~ plot_kmer_dbrn(mutate(.x, direction = plot_group),
                                             rolling_mean = T,
                                             facet_w = "plot_event_type ~ plot_region_type",
                                             n_col = 2,
                                             n_row = 3,
                                             breaks = seq(-500,500,100),
                                             rolling_k = 10,
                                             plot_colours = c("#000000", "#d95f02"),
                                             facet_scales = "fixed"
                            ) + 
                              labs(title = paste(.y, collapse = ","),
                                   y = "% Coverage",
                                   colour = "") +
                              theme(legend.position = "top")
)

# save to PDF, one per page in landscape
plot_list_to_pdf(yg_ya_sum_map_freey,
                 "processed/peka/papa/plots/2023-11-20_papa_all_comparisons_cvcoverage_halleger_yg_map_freey.pdf",
                 width = a4_height,
                 height = a4_width)

plot_list_to_pdf(yg_ya_sum_map_fixed,
                 "processed/peka/papa/plots/2023-11-20_papa_all_comparisons_cvcoverage_halleger_yg_map_fixedy.pdf",
                 width = a4_height,
                 height = a4_width)
