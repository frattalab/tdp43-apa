library(tidyverse)
source("scripts/fncs_plot_peka.R")
source("scripts/fncs_plot_iclip.R")

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


#' go from cv coverage df with % occurrence for each event to a mean count + se confidence interval
#' Assume right now that comparison_name corresponds to the region type, group cryptic/bg
#' Multiple comparison_names/groups appear to work, but no checks made to ensure counts_df and coverage df match properly 
cv_coverage_to_cis <- function(coverage_df, counts_df, join_cols = c("comparison_name", "group"),
                               ci_se_mult = 1.96) {
  
  # add number of events as columns for each position
  joined <- left_join(coverage_df, counts_df, by = join_cols)
  
  # Calculate number of events that have binding from % occurrence
  
  
  # so can generate the distribution of 1s and 0s for each group
  #then can calculate sd and confidence intervals properly
  groupcols <- c(join_cols, "rel_posn")
  joined_summ <- joined %>%
    group_by(across(all_of(groupcols))) %>%
    summarise(occur = (rel_occur / 100) * count,
              rel_occur_est = occur / count,
              rel_occur_orig = unique(rel_occur),
              count_orig = unique(count)) %>%
    ungroup()
  
  # create a df containing each position replicated for the number of events in each group
  # get n rows for number of events in each group
  posns <- joined_summ %>% 
    select(all_of(c(groupcols, "count_orig"))) %>%
    # duplicate rows according to value in count_orig column (number of events)
    uncount(count_orig) %>%
    group_by(across(all_of(groupcols))) %>%
    mutate(row_n = row_number()) %>%
    ungroup()
  
  # create a df containing each position replicated for the number of events with matches in each group
  occur_dbrn <- joined_summ %>%
    select(all_of(c(groupcols, "occur"))) %>%
    mutate(occur = round(occur)) %>%
    # remove rows with no matches
    filter(occur != 0) %>%
    # duplicate the rows/positions based on number of events with match 
    uncount(occur) %>%
    # set dummy value for a match
    mutate(occur = 1) %>%
    group_by(across(all_of(groupcols))) %>%
    # assign ID for matching with total number of events
    mutate(row_n = row_number()) %>%
    ungroup()
  
  
  # combine two, so join and set occur to 1 where matches, otherwise fill with 0 (full dbrn)
  join_cols_dbrn <- c(groupcols, "row_n")
  posns_occur <- left_join(posns, occur_dbrn, by = join_cols_dbrn) %>%
    replace_na(list(occur = 0))
  
  # now calculate confidence intervals based on SEs
  posns_occur %>%
    group_by(across(all_of(groupcols))) %>%
    summarize(avg_coverage = mean(occur),
              se = sd(occur) / sqrt(n()),
              # 95 % confidence intervals
              lwr = avg_coverage - 1.96*se,
              upr = avg_coverage + 1.96*se,
              n_overlaps = sum(occur),
              n_events = n(),
              frac_overlaps = n_overlaps / n_events) %>%
    ungroup()
  
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

# add cleaned comparison name for simpler iCLIP style maps
plot_clean_names <- plot_clean_names %>%
  mutate(plot_comparison_name = paste(plot_event_type, plot_region_type, sep = " - "),
         plot_comparison_name = factor(plot_comparison_name,
                                       levels = c("ALE - Exon Start", "ALE - PAS",
                                                  "IPA - Intron Start", "IPA - PAS",
                                                  "3'Ext - Proximal", "3'Ext - Distal",
                                                  "3'Ext Proximal - Proximal", "3'Ext Proximal - Distal")
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


#### PLOTTING IN ICLIP STYLE FORMAT WITH CONFIDENCE INTERVALS


event_counts_ale <- read_tsv("processed/iclip_regions/2023-12-14_papa_bleedthrough_spliced.event_counts.tsv")
event_counts_d3utr <- read_tsv("processed/iclip_regions/2023-12-15_papa_d3utr.event_counts.tsv")

event_counts_ale %>%
  filter(str_ends(event_type, "all_background$", negate = T))

# subset to bleedtrhoughs, annotating for each bleedthrough category
event_counts_bld <- c("bleedthrough_exonstart", "bleedthrough_pas") %>%
  set_names() %>%
  map(~  event_counts_ale %>%
        filter(event_type == "bleed_out_bed_shsy5y_all") %>%
        select(group = reg_status, count)) %>%
  bind_rows(.id = "comparison_name")

# subset to spliced, annotating for each category
event_counts_spliced <- c("spliced_exonstart", "spliced_pas") %>%
  set_names() %>%
  map(~  event_counts_ale %>%
        filter(event_type == "spliced_shsy5y_background") %>%
        select(group = reg_status, count)) %>%
  bind_rows(.id = "comparison_name")

# subset to 3'exts, annotating for each comparison 
event_counts_d3utr <- c("d3utr_pas_proximal", "d3utr_pas_distal") %>%
  set_names() %>%
  map(~ event_counts_d3utr %>%
        select(group = reg_status, count) %>%
        distinct(.keep_all = T)
      ) %>%
  bind_rows(.id = "comparison_name")

# TODO: proximal 3'UTR extensions are currently missing, need to extract event counts from scropt
# BUT as know identical background and cryptic n = 20, can construct for now

event_counts_d3utrprox <- tibble("comparison_name" = c(rep("proximald3utr_pas_proximal", 2), rep("proximald3utr_pas_distal",2)),
                                 "group" = rep(c("background", "cryptic"), 2),
                                 "count" = rep(c(798, 20), 2)
                                 )

# combine all event counts dfs
event_counts_all <- bind_rows(event_counts_spliced, event_counts_bld, event_counts_d3utr, event_counts_d3utrprox)

# porcess cv coverage tables to format compatible with iCLIP maps (i.e. mean coverage + calculate standard errors)
dbrn_tbls_sum_ci <- map(dbrn_tbls_sum,
    ~ .x %>%
      mutate(group = if_else(group == "foreground", "cryptic", group)) %>%
      cv_coverage_to_cis(., event_counts_all, ci_se_mult = 1) %>%
      left_join(plot_clean_names, by = "comparison_name")
    )


# Make combined map of all comparisons on a single plot, separately for each motif group
separate_sum_maps_freey_ci <- map2(.x = dbrn_tbls_sum_ci,
     .y = names(dbrn_tbls_sum_ci),
    ~ .x %>%
      rename(position = rel_posn) %>%
      mutate(plot_cryptic = if_else(group == "cryptic", "Cryptic", "Background"),
             plot_type = plot_comparison_name,
             position = position + 500
             ) %>%
      plot_coverage(ci_se_mult = 1,
                    loess_span = 0.1,
                    y_scales = scale_y_continuous(),
                    title_lab = .y,
                    facet_scales = "free_y")
    )


# save to PDF, one per page in landscape
plot_list_to_pdf(separate_sum_maps_freey_ci,
                 "processed/peka/papa/plots/2023-12-20_papa_cvcoverage_halleger_motif_groups_summed_map_freey_ci.pdf",
                 width = a4_height,
                 height = a4_width)

# generate single panel plots per region for YG group
per_region_yg_group_sum_maps_freey_ci <- unique(dbrn_tbls_sum_ci$`YG-containing-motifs`$comparison_name) %>%
  set_names() %>%
  map(  ~ dbrn_tbls_sum_ci$`YG-containing-motifs` %>%
  filter(comparison_name == .x) %>%
  rename(position = rel_posn) %>%
  mutate(plot_cryptic = if_else(group == "cryptic", "Cryptic", "Background"),
         plot_type = plot_region_type,
         position = position + 500
  ) %>%
  plot_coverage(ci_se_mult = 1,
                loess_span = 0.1,
                y_scales = scale_y_continuous(),
                facet_scales = "free_y")
  )

# save to SVG
walk2(.x = per_region_yg_group_sum_maps_freey_ci,
      .y = names(per_region_yg_group_sum_maps_freey_ci),
      ~ ggsave(filename = paste("2023-12-20_papa_cvcoverage.yg_group_summed_map.",
                                .y,
                                ".ci_freey.svg",
                                sep = ""),
               plot = .x,
               path = "processed/peka/papa/plots/",
                              device = svg,
                              height = 6*0.8,
                              width = 18*0.8,
                              units = "in",
                              dpi = "retina")
      )


# Finally, make horizontally faceted plots of each event type for YG group
per_event_yg_group_sum_maps_fixedy_ci <- unique(dbrn_tbls_sum_ci$`YG-containing-motifs`$plot_event_type) %>%
  set_names() %>%
  map(~ dbrn_tbls_sum_ci$`YG-containing-motifs` %>%
        filter(plot_event_type == .x) %>%
        rename(position = rel_posn) %>%
        mutate(plot_cryptic = if_else(group == "cryptic", "Cryptic", "Background"),
               plot_type = plot_region_type,
               position = position + 500
        ) %>%
        plot_coverage(ci_se_mult = 1,
                      loess_span = 0.1,
                      y_scales = scale_y_continuous(),
                      title_lab = .x,
                      facet_scales = "fixed")
      )
        

# save to SVG
walk2(.x = per_event_yg_group_sum_maps_fixedy_ci,
      .y = names(per_event_yg_group_sum_maps_fixedy_ci),
      ~ ggsave(filename = paste("2023-12-20_papa_cvcoverage.yg_group_summed_map_horiz_stack.",
                                str_replace_all(.y, "'|-|[[:space:]]","_"),
                                ".ci_fixed.svg",
                                sep = ""),
               plot = .x,
               path = "processed/peka/papa/plots/",
               device = svg,
               height = 6*0.8,
               width = 18*0.8,
               units = "in",
               dpi = "retina")
               )



