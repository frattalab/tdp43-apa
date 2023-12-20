library(tidyverse)
source("scripts/fncs_plot_iclip.R")

parse_coverage <- function(file, flank_interval) {
  f <- read_tsv(file, col_names = c("chr", "start", "end", "ID", ".", "strand", "position", "coverage"))
  # make sure positions are strand-aware
  f$position[f$strand == "-"] <- (2 + flank_interval*2) - f$position[f$strand == "-"]
  
  average_coverage <- f %>%
    group_by(position) %>%
    summarize(avg_coverage = mean(coverage),
              se = sd(coverage) / sqrt(n()),
              # 95 % confidence intervals
              lwr = avg_coverage - 1.96*se, 
              upr = avg_coverage + 1.96*se,
              n_overlaps = sum(coverage),
              n_events = n(),
              frac_overlaps = sum(coverage) / n()
    )
  average_coverage
  
}

#' read in coverages files - assumes a named vector
get_combined_coverages <- function(files, flank_interval) {
  
  files %>%
    map(~ parse_coverage(.x, 500)) %>%
    bind_rows(.id = "origin") %>%
    # pull out site type & status
    separate(origin, into = c("background_type", "type", "region_type", "cryptic"), sep = "\\.", remove = T)
  
}



d3utr_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/d3utr", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 
spliced_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/spliced", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 
bleedthrough_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/bleedthrough_uniq", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 
d3utrprox_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/d3utr_proximal", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 

# Extract info from file path to construct IDs
# background_type.event_type.region_type.cryptic/bg

# processed/iclip_maps/coverage/background_shsy5y/d3utr/distal/regions.flank_500.coverage.bg.txt.gz -> background_shsy5y.d3utr.distal.bg 
d3utr_nm <- paste(str_remove(str_replace_all(dirname(d3utr_paths), "\\/", "."), "processed.iclip_maps.coverage."),
                  str_extract(basename(d3utr_paths), pattern = "bg|cr"),
                  sep = "."
                  )

# processed/iclip_maps/coverage/background_shsy5y/d3utr/distal/regions.flank_500.coverage.bg.txt.gz -> background_shsy5y.d3utr.distal.bg 
d3utrprox_nm <- paste(str_remove(str_replace_all(dirname(d3utrprox_paths), "\\/", "."),
                                 "processed.iclip_maps.coverage."),
                      str_extract(basename(d3utrprox_paths), pattern = "bg|cr"),
                      sep = "."
                      )

# background_shsy5y/spliced/pas/regions.flank_500.coverage.bg.txt.gz" -> background_shsy5y.spliced.pas.bg
spliced_nm <- paste(str_remove(str_replace_all(dirname(spliced_paths), "\\/", "."), "processed.iclip_maps.coverage."),
                    str_extract(basename(spliced_paths), pattern = "bg|cr"),
                    sep = "."
                    )

bleedthrough_nm <- paste(str_remove(str_replace_all(dirname(bleedthrough_paths), "\\/", "."), "processed.iclip_maps.coverage."),
                    str_extract(basename(bleedthrough_paths), pattern = "bg|cr"),
                    sep = "."
                    )



# read in coverages to combined dfs

spliced_average_coverage <- get_combined_coverages(set_names(spliced_paths, spliced_nm), 500)
bleedthrough_average_coverage <- get_combined_coverages(set_names(bleedthrough_paths, bleedthrough_nm), 500)
d3utr_average_coverage <- get_combined_coverages(set_names(d3utr_paths, d3utr_nm), 500)
d3utrprox_average_coverage <- get_combined_coverages(set_names(d3utrprox_paths, d3utrprox_nm), 500)
# Group specific tidying for plotting

d3utr_average_coverage <- d3utr_average_coverage %>%
  mutate(plot_type = if_else(region_type == "proximal",
                             "Proximal",
                             "Distal"),
         plot_type = factor(plot_type, levels = c("Proximal", "Distal")),
         plot_cryptic = if_else(cryptic == "cr", "Cryptic", "Background"),
         plot_cryptic = factor(plot_cryptic, levels = c("Background", "Cryptic"))
         )

d3utrprox_average_coverage <- d3utrprox_average_coverage %>%
  mutate(plot_type = if_else(region_type == "proximal",
                             "Proximal",
                             "Distal"),
         plot_type = factor(plot_type, levels = c("Proximal", "Distal")),
         plot_cryptic = if_else(cryptic == "cr", "Cryptic", "Background"),
         plot_cryptic = factor(plot_cryptic, levels = c("Background", "Cryptic"))
  )

spliced_average_coverage <- spliced_average_coverage %>%
  mutate(plot_type = if_else(region_type == "le_start", "Exon Start", "PAS"),
         plot_type = factor(plot_type, levels = c("Exon Start", "PAS")),
         plot_cryptic = if_else(cryptic == "cr", "Cryptic", "Background"),
         plot_cryptic = factor(plot_cryptic, levels = c("Background", "Cryptic"))
         )

bleedthrough_average_coverage <- bleedthrough_average_coverage %>%
  mutate(plot_type = if_else(region_type == "le_start", "Intron Start", "PAS"),
         plot_type = factor(plot_type, levels = c("Intron Start", "PAS")),
         plot_cryptic = if_else(cryptic == "cr", "Cryptic", "Background"),
         plot_cryptic = factor(plot_cryptic, levels = c("Background", "Cryptic"))
         )

# construct list of events to facilitate generating common plots
event_lists <- list("3'Ext" = d3utr_average_coverage,
                    "3'Ext Proximal" = d3utrprox_average_coverage,
                    "ALE" = spliced_average_coverage,
                    "IPA" = bleedthrough_average_coverage)


# More liberal confidence intervals (1*se) - faceted/side-by-side
iclip_maps_1se <- map2(.x = event_lists, .y = names(event_lists),
     ~ plot_coverage(.x, ci_se_mult = 1, title_lab = .y, loess_span = 0.1,
                     y_scales = scale_y_continuous(limits = c(NA, 0.12),
                                                   breaks = seq(0, 0.12, 0.02))
                     )
)

# get underlying dfs
iclip_dfs_maps_1se <- map2(.x = event_lists, .y = names(event_lists),
                           ~ plot_coverage_df(.x, ci_se_mult = 1)
)


# generate single event type maps (for spliced & bleedthrough)
iclip_maps_single_le_start_1se <- map2(.x = event_lists[3:length(event_lists)],
                                       .y = names(event_lists)[3:length(event_lists)],
                                       ~ plot_coverage(filter(.x, region_type == "le_start"),
                                                       ci_se_mult = 1,
                                                       title_lab = .y, loess_span = 0.1,
                                                       y_scales = scale_y_continuous(limits = c(NA, 0.12),
                                                                                     breaks = seq(0, 0.16, 0.02)))
                                       )

iclip_maps_single_pas_1se <- map2(.x = event_lists[3:length(event_lists)],
                                  .y = names(event_lists)[3:length(event_lists)],
                                  ~ plot_coverage(filter(.x, region_type == "pas"), ci_se_mult = 1, title_lab = .y, loess_span = 0.1,
                                                  y_scales = scale_y_continuous(limits = c(NA, 0.12),
                                                                                breaks = seq(0, 0.16, 0.02)))
                                  )

# # now make single panel plots for 3'UTR-APA
iclip_maps_single_d3utr_1se <- c("proximal", "distal") %>%
  set_names() %>%
  map(~ plot_coverage(filter(d3utr_average_coverage, region_type == .x),
                      ci_se_mult = 1,
                      title_lab = "3'Ext", loess_span = 0.1,
                      y_scales = scale_y_continuous(limits = c(NA, 0.12),
                                                    breaks = seq(0, 0.16, 0.02)))
      )

iclip_maps_single_d3utrprox_1se <- c("proximal", "distal") %>%
  set_names() %>%
  map(~ plot_coverage(filter(d3utrprox_average_coverage, region_type == .x),
                      ci_se_mult = 1,
                      title_lab = "3'Ext Proximal", loess_span = 0.1,
                      y_scales = scale_y_continuous(limits = c(NA, 0.12),
                                                    breaks = seq(0, 0.16, 0.02)))
      )


if (!dir.exists("processed/iclip_maps/plots")) { dir.create("processed/iclip_maps/plots", recursive = T)}

# write to file (PNG and SVG)
walk2(.x = iclip_maps_1se,
      .y = names(iclip_maps_1se),
      ~ ggsave(filename = paste("2023-12-20_background_shsy5y_papa_cryptic_iclip_map.horiz_stack.fixed_ylim.1_se_ci.",
                                str_replace_all(.y, "'|-", "_"),
                                ".png",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = "png",
               height = 6 * 0.8,
               width = 18 * 0.8,
               units = "in",
               dpi = "retina")
)

walk2(.x = iclip_maps_1se,
      .y = names(iclip_maps_1se),
      ~ ggsave(filename = paste("2023-12-20_background_shsy5y_papa_cryptic_iclip_map.horiz_stack.fixed_ylim.1_se_ci.",
                                str_replace_all(.y, "'|-", "_"),
                                ".svg",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = svg,
               height = 6 * 0.8,
               width = 18 * 0.8,
               units = "in",
               dpi = "retina")
)

# write plot dfs to file
walk2(.x = iclip_dfs_maps_1se,
      .y = names(iclip_dfs_maps_1se),
      ~ write_tsv(.x, 
                  paste("processed/iclip_maps/plots/2023-12-20_background_shsy5y_papa_cryptic_iclip_df.horiz_stack.fixed_ylim.1_se_ci.",
                            str_replace_all(.y, "'|-", "_"),
                            ".tsv",
                            sep = ""),
                  col_names = T
                  )
      )


# write single panel d3'utrs to file
walk2(.x = iclip_maps_single_d3utr_1se,
      .y = names(iclip_maps_single_d3utr_1se),
      ~ ggsave(filename = paste("2023-12-20_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.",
                                .y,
                                ".svg",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = svg,
               height = 6 * 0.8,
               width = 18 * 0.8,
               units = "in",
               dpi = "retina")
      )

walk2(.x = iclip_maps_single_d3utr_1se,
      .y = names(iclip_maps_single_d3utr_1se),
      ~ ggsave(filename = paste("2023-12-20_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.",
                                .y,
                                ".png",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = "png",
               height = 6 * 0.8,
               width = 18 * 0.8,
               units = "in",
               dpi = "retina")
)

# proximal cryptics
walk2(.x = iclip_maps_single_d3utrprox_1se,
      .y = names(iclip_maps_single_d3utrprox_1se),
      ~ ggsave(filename = paste("2023-12-20_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.",
                                .y,
                                ".svg",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = svg,
               height = 6 * 0.8,
               width = 18 * 0.8,
               units = "in",
               dpi = "retina")
)

walk2(.x = iclip_maps_single_d3utrprox_1se,
      .y = names(iclip_maps_single_d3utrprox_1se),
      ~ ggsave(filename = paste("2023-12-20_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.",
                                .y,
                                ".png",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = "png",
               height = 6 * 0.8,
               width = 18 * 0.8,
               units = "in",
               dpi = "retina")
)

# write single panel le starts to file
walk2(.x = iclip_maps_single_le_start_1se,
      .y = names(iclip_maps_single_le_start_1se),
      ~ ggsave(filename = paste("2023-12-20_background_shsy5y_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.le_start.",
                                str_replace_all(.y, "'|-", "_"),
                                ".svg",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = svg,
               height = 6 * 0.8,
               width = 18 * 0.8,
               units = "in",
               dpi = "retina")
        )

# write single panel PAS to file
walk2(.x = iclip_maps_single_pas_1se,
      .y = names(iclip_maps_single_pas_1se),
      ~ ggsave(filename = paste("2023-12-20_background_shsy5y_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.pas.",
                                str_replace_all(.y, "'|-", "_"),
                                ".svg",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = svg,
               height = 6 * 0.8,
               width = 18 * 0.8,
               units = "in",
               dpi = "retina")
)
