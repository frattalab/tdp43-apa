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

# Set the A4 paper size dimensions in inches
a4_width <- 8.27
a4_height <- 11.69

# paths to different processed kmer distribution tbls
dbrn_tbl_paths <- list("YG-containing-motifs" = "processed/peka/papa/2023-12-15_papa_cryptics_cvcoverage_window_500.yg_6mers.per_kmer_distribution_genome.tsv",
                   "YA-containing-motifs" = "processed/peka/papa/2023-12-15_papa_cryptics_cvcoverage_window_500.ya_6mers.per_kmer_distribution_genome.tsv",
                   "AA-containing-motifs" = "processed/peka/papa/2023-12-15_papa_cryptics_cvcoverage_window_500.aa_6mers.per_kmer_distribution_genome.tsv")

# dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-27_papa_cryptics_cvcoverage_window_500.yg_6mers.per_kmer_distribution_genome.tsv")

# map comparison names to cleaned event types & region types
plot_clean_names <- tibble(comparison_name = c("bleedthrough_exonstart", "bleedthrough_pas", "spliced_exonstart", "spliced_pas", "d3utr_pas_proximal", "d3utr_pas_distal", "proximald3utr_pas_proximal", "proximald3utr_pas_distal"),
                           plot_event_type = c("IPA", "IPA", "ALE", "ALE", "3'Ext", "3'Ext", "3'Ext Proximal", "3'Ext Proximal"),
                           plot_region_type = factor(c("Intron Start", "PAS", "Exon Start", "PAS", "Proximal", "Distal", "Proximal", "Distal"),
                                                     levels = c("Intron Start","Exon Start", "PAS", "Proximal", "Distal")
                           )
)

# 6mer groups enriched in Halleger et al. TDP-43 + CTD MUT iCLIP
# Create list ready for use with fgsea
tdp_motif_groups <- list("YG-containing-motifs" = c("UGUGUG", "GUGUGU","UGUGCG", "UGCGUG","CGUGUG","GUGUGC"),
                         "YA-containing-motifs" = c("AUGUGU", "GUAUGU", "GUGUAU", "UGUGUA", "UGUAUG", "UGCAUG"),
                         "AA-containing-motifs" = c("GUGUGA", "AAUGAA", "GAAUGA", "UGAAUG", "AUGAAU", "GUGAAU", "GAAUGU", "UUGAAU")
)

# read in dfs, add plot labels & covert coverage to % values, remove bleedthroughs
dbrn_tbls <- dbrn_tbl_paths %>%
  map(~ read_tsv(.x, show_col_types = F) %>%
        separate(region_type, into = c("comparison_name", "background_type"), sep = "\\.", remove = F) %>%
        left_join(plot_clean_names, by = "comparison_name") %>%
        mutate(plot_group = if_else(group == "foreground", "Cryptic", "Background"),
               rel_occur = rel_occur * 100) # %>%
        # filter(str_starts(region_type, "bleedthrough", negate = T))
      )

# make dfs with summed occurence within each kmer group
dbrn_tbls_sum <- map(dbrn_tbls,
    ~ .x %>%
      group_by(comparison_name, background_type, kmer_group, group, rel_posn) %>%
      summarise(rel_occur = sum(rel_occur)) %>%
      ungroup() %>%
      # add back in the plot name variables
      left_join(plot_clean_names, by = "comparison_name") %>%
      mutate(plot_group = if_else(group == "foreground", "Cryptic", "Background"))
  )


# make maps with sum of kmers, separately for each kmer group
separate_sum_maps_freey <-  map2(.x = dbrn_tbls_sum,
     .y = names(dbrn_tbls_sum),
     ~ plot_kmer_dbrn(mutate(.x, direction = plot_group),
                      rolling_mean = T,
                      facet_w = "plot_event_type ~ plot_region_type",
                      n_col = 2,
                      n_row = 4,
                      breaks = seq(-500,500,100),
                      rolling_k = 12,
                      plot_colours = c("#000000", "#d95f02"),
                      base_size = 20
                      ) + 
       labs(title = .y,
            y = "% Coverage",
            colour = "") +
       theme(legend.position = "top")
       )

separate_sum_maps_fixed <-  map2(.x = dbrn_tbls_sum,
                                 .y = names(dbrn_tbls_sum),
                                 ~ plot_kmer_dbrn(mutate(.x, direction = plot_group),
                                                  rolling_mean = T,
                                                  facet_w = "plot_event_type ~ plot_region_type",
                                                  n_col = 2,
                                                  n_row = 4,
                                                  breaks = seq(-500,500,100),
                                                  rolling_k = 12,
                                                  plot_colours = c("#000000", "#d95f02"),
                                                  facet_scales = "fixed",
                                                  base_size = 20
                                 ) + 
                                   labs(title = .y,
                                        y = "% Coverage",
                                        colour = "") +
                                   theme(legend.position = "top")
)

separate_sum_maps_freey

# save to PDF, one per page in landscape
plot_list_to_pdf(separate_sum_maps_freey,
                 "processed/peka/papa/plots/2023-12-17_papa_cvcoverage_halleger_motif_groups_summed_map_freey.pdf",
                 width = a4_height,
                 height = a4_width)

plot_list_to_pdf(separate_sum_maps_fixed,
                 "processed/peka/papa/plots/2023-12-17_papa_cvcoverage_halleger_motif_groups_summed_map_fixedy.pdf",
                 width = a4_height,
                 height = a4_width)


# Also save individual plots to SVG
walk2(.x = separate_sum_maps_freey,
      .y = names(separate_sum_maps_freey),
      ~ ggsave(filename = paste("2023-12-17_papa_cvcoverage.",
                                .y,
                                ".summed_map_freey.svg",
                                sep = ""),
               plot = .x + labs(title = "", x = "Position"),
               path = "processed/peka/papa/plots/",
               device = svg,
               width = 19.5,
               height = 6.5*2,
               units = "in",
               dpi = "retina")
      )

walk2(.x = separate_sum_maps_fixed,
      .y = names(separate_sum_maps_fixed),
      ~ ggsave(filename = paste("2023-12-17_papa_cvcoverage.",
                                .y,
                                ".summed_map_fixedy.svg",
                                sep = ""),
               plot = .x + labs(title = "", x = "Position"),
               path = "processed/peka/papa/plots/",
               device = svg,
               width = 19.5,
               height = 6.5*2,
               units = "in",
               dpi = "retina")
)


# # map with free scales across y
# map_freey <- plot_kmer_dbrn(mutate(dbrn_tbl, direction = plot_group),
#                rolling_mean = T,
#                facet_w = "plot_event_type ~ plot_region_type",
#                n_col = 2,
#                n_row = 3,
#                breaks = seq(-500,500,100),
#                rolling_k = 12,
#                plot_colours = c("#000000", "#d95f02")
#                ) + 
#   labs(y = "% Coverage",
#        colour = "") +
#   theme(legend.position = "top")
# 
# # map with fixed scales across events
# map_fixed <- plot_kmer_dbrn(mutate(dbrn_tbl, direction = plot_group),
#                rolling_mean = T,
#                facet_w = "plot_event_type ~ plot_region_type",
#                n_col = 2,
#                n_row = 3,
#                breaks = seq(-500,500,100),
#                rolling_k = 12,
#                plot_colours = c("#000000", "#d95f02"),
#                facet_scales = "fixed") + 
#   labs(y = "% Coverage",
#        colour = "") +
#   theme(legend.position = "top")

# map_freey
# map_fixed
# 

# 
# # save each plot as 1 page per a4
# # Create a multi-page PDF file for each group of plots
# if (!dir.exists("processed/peka/papa/plots/")) { dir.create("processed/peka/papa/plots/", recursive = T)}
# 
# plots_ug_list <- list(map_freey, map_fixed)
# pdf("processed/peka/papa/plots/2023-11-16_papa_all_comparisons_cvcoverage_gu_ug_map.pdf",
#     width = a4_height, height = a4_width )
# 
# # Loop through each ggplot object and print it to the PDF file
# 
# for (i in seq_along(plots_ug_list)) {
#   print(plots_ug_list[[i]])
# }
# 
# # Close the PDF file
# dev.off()

# #  repeat for Halleger et al motif groups
# yg_dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-19_papa_cryptics_cvcoverage_window_500.yg_6mers.per_kmer_distribution_genome.tsv")
# ya_dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-19_papa_cryptics_cvcoverage_window_500.ya_6mers.per_kmer_distribution_genome.tsv")
# yg_6mers <- c("UGUGUG", "GUGUGU","UGUGCG", "UGCGUG","CGUGUG","GUGUGC")
# ya_6mers <- c("AUGUGU", "GUAUGU", "GUGUAU", "UGUGUA", "UGUAUG", "UGCAUG")
# 
# # Sum occurrence of kmers across the yg groups
# yg_sum_dbrn_tbl <-  yg_dbrn_tbl %>%
#   group_by(region_type, group, rel_posn) %>%
#   summarise(rel_occur = sum(rel_occur)) %>%
#   ungroup()
# 
# ya_sum_dbrn_tbl <-  ya_dbrn_tbl %>%
#   group_by(region_type, group, rel_posn) %>%
#   summarise(rel_occur = sum(rel_occur)) %>%
#   ungroup()
# 
# # prepare for ploitting - column names, cryptic/not etc. 
# 
# yg_sum_dbrn_tbl <- yg_sum_dbrn_tbl %>%
#   rename(comparison_name = region_type) %>%
#   left_join(plot_clean_names, by = "comparison_name") %>%
#   mutate(plot_group = if_else(group == "foreground", "Cryptic", "Background"),
#          rel_occur = rel_occur * 100)
# 
# ya_sum_dbrn_tbl <- ya_sum_dbrn_tbl %>%
#   rename(comparison_name = region_type) %>%
#   left_join(plot_clean_names, by = "comparison_name") %>%
#   mutate(plot_group = if_else(group == "foreground", "Cryptic", "Background"),
#          rel_occur = rel_occur * 100)
# 
# # maps with free scales across y
# yg_ya_sum_map_freey <- map2(.x = list(YG_group = yg_sum_dbrn_tbl, YA_group = ya_sum_dbrn_tbl),
#                             .y = list(YG_group = yg_6mers, YA_group = ya_6mers),
#                             ~ plot_kmer_dbrn(mutate(.x, direction = plot_group),
#                                              rolling_mean = T,
#                                              facet_w = "plot_event_type ~ plot_region_type",
#                                              n_col = 2,
#                                              n_row = 3,
#                                              breaks = seq(-500,500,100),
#                                              rolling_k = 12,
#                                              plot_colours = c("#000000", "#d95f02")
#                             ) + 
#                               labs(title = paste(.y, collapse = ","),
#                                    y = "% Coverage",
#                                    colour = "") +
#                               theme(legend.position = "top")
# )
# 
# # maps with fixed scales across y
# yg_ya_sum_map_fixed <- map2(.x = list(YG_group = yg_sum_dbrn_tbl, YA_group = ya_sum_dbrn_tbl),
#                             .y = list(YG_group = yg_6mers, YA_group = ya_6mers),
#                             ~ plot_kmer_dbrn(mutate(.x, direction = plot_group),
#                                              rolling_mean = T,
#                                              facet_w = "plot_event_type ~ plot_region_type",
#                                              n_col = 2,
#                                              n_row = 3,
#                                              breaks = seq(-500,500,100),
#                                              rolling_k = 12,
#                                              plot_colours = c("#000000", "#d95f02"),
#                                              facet_scales = "fixed"
#                             ) + 
#                               labs(title = paste(.y, collapse = ","),
#                                    y = "% Coverage",
#                                    colour = "") +
#                               theme(legend.position = "top")
# )
# 
# # save to PDF, one per page in landscape
# plot_list_to_pdf(yg_ya_sum_map_freey,
#                  "processed/peka/papa/plots/2023-11-20_papa_all_comparisons_cvcoverage_halleger_yg_map_freey.pdf",
#                  width = a4_height,
#                  height = a4_width)
# 
# plot_list_to_pdf(yg_ya_sum_map_fixed,
#                  "processed/peka/papa/plots/2023-11-20_papa_all_comparisons_cvcoverage_halleger_yg_map_fixedy.pdf",
#                  width = a4_height,
#                  height = a4_width)


# What about all groups on one plot, separated 



### repeat for per-motif (attempted, but cancelling as want split by background/cryptic still)
# line-type to split cryptic vs background?
# 

# clean up dfs ready for plotting
# yg_dbrn_tbl <- yg_dbrn_tbl %>%
# rename(comparison_name = region_type) %>%
#   left_join(plot_clean_names, by = "comparison_name") %>%
#   mutate(plot_group = if_else(group == "foreground", "Cryptic", "Background"),
#          rel_occur = rel_occur * 100,
#          direction = kmer)
# 
# ya_dbrn_tbl <- ya_dbrn_tbl %>%
#   rename(comparison_name = region_type) %>%
#   left_join(plot_clean_names, by = "comparison_name") %>%
#   mutate(plot_group = if_else(group == "foreground", "Cryptic", "Background"),
#          rel_occur = rel_occur * 100,
#          direction = kmer)
# 
# 
# 
# # maps with fixed scales across y
# yg_ya_perk_map_fixed <- map2(.x = list(YG_group = yg_dbrn_tbl, YA_group = ya_dbrn_tbl),
#                             .y = list(YG_group = yg_6mers, YA_group = ya_6mers),
#                             ~ plot_kmer_dbrn(.x,
#                                              rolling_mean = T,
#                                              facet_w = "plot_event_type ~ plot_region_type",
#                                              n_col = 2,
#                                              n_row = 3,
#                                              breaks = seq(-500,500,100),
#                                              rolling_k = 12,
#                                              plot_colours = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'),
#                             ) + 
#                               labs(title = paste(.y, collapse = ","),
#                                    y = "% Coverage",
#                                    colour = "") +
#                               theme(legend.position = "top")
# )
# 
# yg_ya_perk_map_fixed$YG_group
# yg_ya_perk_map_fixed$YA_group
