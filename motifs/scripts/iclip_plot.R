library(tidyverse)

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

plot_coverage_df <- function(df, ci_se_mult = 1.96, event_col = "plot_type", group_col = "plot_cryptic") {
  
  group_cols <- c(event_col, group_col)
  
  # generate confidence interval values
  df %>%
    mutate(plot_ymin = avg_coverage - (ci_se_mult*se),
           plot_ymax = avg_coverage + (ci_se_mult*se)) %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(ymin_smooth = stats::predict(loess(plot_ymin~position, span=0.1)),
           ymax_smooth = stats::predict(loess(plot_ymax~position, span=0.1))) %>%
    ungroup()
  
}

plot_coverage <- function(df, ci_se_mult = 1.96, event_col = "plot_type", group_col = "plot_cryptic", facet_ncol = 2, fill_colours = c("#000000", "#d95f02"), line_colours = c("#000000", "#d95f02"), fill_lab = "", colour_lab = "", title_lab = "") {
  
  # generate confidence interval values
  plot_df <- plot_coverage_df(df, ci_se_mult, event_col, group_col)
  
  plot_df %>%
    ggplot(aes(x = position, y = avg_coverage, color=!!sym(group_col), fill=!!sym(group_col))) +
    geom_smooth(method = "loess", span=0.2, se = F) +
    geom_ribbon(aes(ymin = ymin_smooth,  ymax = ymax_smooth, fill = !!sym(group_col)),
                alpha = 0.31) +
    xlab("Position") +
    ylab("Average Coverage") +
    geom_vline(xintercept = 500, linetype = "dashed", alpha = 0.5) +
    facet_wrap(event_col, ncol = facet_ncol, scales = "fixed") +
    scale_x_continuous(
      limits = c(0, 1001),
      breaks = seq(0,1000,100),
      labels = as.character(seq(-500,500,100))
    ) +
    scale_y_continuous(limits = c(NA, 0.1),
                       breaks = seq(0, 0.1, 0.02)) +
    scale_fill_manual(values = fill_colours) +
    scale_color_manual(values = line_colours) +
    theme_bw(base_size = 20) +
    theme(legend.position = "top",
          axis.text = element_text(size = rel(1.25)),
          axis.title = element_text(size = rel(1.25)),
          strip.text = element_text(size = rel(1.25))
          ) +
    labs(fill = fill_lab, colour = colour_lab, title = title_lab)
  
  
}


#' read in coverages files - assumes a named vector
get_combined_coverages <- function(files, flank_interval) {
  
  files %>%
    map(~ parse_coverage(.x, 500)) %>%
    bind_rows(.id = "origin") %>%
    # pull out site type & status
    separate(origin, into = c("background_type", "type", "region_type", "cryptic"), sep = "\\.", remove = T)
  
}



# d3utr_paths <- list.files(path = "processed/iclip_maps/d3utr", full.names = T) 
spliced_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/spliced", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 
bleedthrough_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/bleedthrough_uniq", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 

# extract coord type & group from filenames
# split by '.', extract the two elements before '.txt.gz'
# d3utr_nm <- str_split(basename(d3utr_paths), "\\.",simplify = T) %>%
#   apply(MARGIN = 1,
#         function(x) {paste(c(x[length(x)-3], x[length(x)-2]),
#                            collapse = ".")
#           })



# background_shsy5y/spliced/pas/regions.flank_500.coverage.bg.txt.gz" -> background_shsy5y.spliced.pas.bg
spliced_nm <- paste(str_remove(str_replace_all(dirname(spliced_paths), "\\/", "."), "processed.iclip_maps.coverage."),
                    str_extract(basename(spliced_paths), pattern = "bg|cr"),
                    sep = "."
                    )

bleedthrough_nm <- paste(str_remove(str_replace_all(dirname(bleedthrough_paths), "\\/", "."), "processed.iclip_maps.coverage."),
                    str_extract(basename(bleedthrough_paths), pattern = "bg|cr"),
                    sep = "."
                    )
# 
# # 2023-07-04_papa_cryptic_spliced.le_start.flank_500.coverage.bg.txt.gz -> le_start.bg
# spliced_nm <- str_split(basename(spliced_paths), "\\.",simplify = T) %>% 
#   apply(MARGIN = 1,
#         function(x) {paste(c(x[length(x)-5], x[length(x)-2]),
#                            collapse = ".")
#         })
# 
# # processed/iclip_maps/bleedthrough/2023-07-04_papa_cryptic_bleedthrough.le_start.flank_500.coverage.bg.txt.gz -> le_start.bg
# bleedthrough_nm <- str_split(basename(bleedthrough_paths), "\\.",simplify = T) %>% 
#   apply(MARGIN = 1,
#         function(x) {paste(c(x[length(x)-5], x[length(x)-2]),
#                            collapse = ".")
#         })



# read in coverages to combined dfs
# d3utr_average_coverage <- get_combined_coverages(set_names(d3utr_paths, d3utr_nm), 500)
spliced_average_coverage <- get_combined_coverages(set_names(spliced_paths, spliced_nm), 500)
bleedthrough_average_coverage <- get_combined_coverages(set_names(bleedthrough_paths, bleedthrough_nm), 500)


# d3utr - group specific tidying for plotting
# d3utr_average_coverage <- d3utr_average_coverage %>%
#   mutate(plot_type = if_else(type == "prox", "Proximal",
#                              "Distal"),
#          plot_type = factor(plot_type, levels = c("Proximal", "Distal")),
#          plot_cryptic = if_else(cryptic == "cr", "Cryptic", "Background"),
#          plot_cryptic = factor(plot_cryptic, levels = c("Background", "Cryptic"))
#          )

spliced_average_coverage <- spliced_average_coverage %>%
  mutate(plot_type = if_else(region_type == "le_start", "Exon Start", "PAS"),
         plot_type = factor(plot_type, levels = c("Exon Start", "PAS")),
         plot_cryptic = if_else(cryptic == "cr", "Cryptic", "Background"),
         plot_cryptic = factor(plot_cryptic, levels = c("Background", "Cryptic")))

bleedthrough_average_coverage <- bleedthrough_average_coverage %>%
  mutate(plot_type = if_else(region_type == "le_start", "Exon Start", "PAS"),
         plot_type = factor(plot_type, levels = c("Exon Start", "PAS")),
         plot_cryptic = if_else(cryptic == "cr", "Cryptic", "Background"),
         plot_cryptic = factor(plot_cryptic, levels = c("Background", "Cryptic")))


# "3'UTR-ALE" = d3utr_average_coverage,
event_lists <- list(
                    "AS-ALE" = spliced_average_coverage,
                    "Bleedthrough-ALE" = bleedthrough_average_coverage)


# More liberal confidence intervals (1*se) - faceted/side-by-side
iclip_maps_1se <- map2(.x = event_lists, .y = names(event_lists),
     ~ plot_coverage(.x, ci_se_mult = 1, title_lab = .y)
)

# get underlying dfs
iclip_dfs_maps_1se <- map2(.x = event_lists, .y = names(event_lists),
                           ~ plot_coverage_df(.x, ci_se_mult = 1)
)


# generate single event type maps (for spliced & bleedthrough)
# [2:length(event_lists)]
iclip_maps_single_le_start_1se <- map2(.x = event_lists,
                                       .y = names(event_lists),
                                       ~ plot_coverage(filter(.x, region_type == "le_start"),
                                                       ci_se_mult = 1,
                                                       title_lab = .y)
                                       )

iclip_maps_single_pas_1se <- map2(.x = event_lists,
                                  .y = names(event_lists),
                                  ~ plot_coverage(filter(.x, region_type == "pas"), ci_se_mult = 1, title_lab = .y)
                                  )

# # now make single panel plots for 3'UTR-APA
# iclip_maps_single_d3utr_1se <- c("prox", "dist") %>%
#   set_names() %>%
#   map(~ plot_coverage(filter(d3utr_average_coverage, type == .x),
#                       ci_se_mult = 1,
#                       title_lab = "3'UTR-APA")
#       )


if (!dir.exists("processed/iclip_maps/plots")) { dir.create("processed/iclip_maps/plots", recursive = T)}

# write to file (PNG and SVG)
walk2(.x = iclip_maps_1se,
      .y = names(iclip_maps_1se),
      ~ ggsave(filename = paste("2023-11-28_background_shsy5y_papa_cryptic_iclip_map.horiz_stack.fixed_ylim.1_se_ci.",
                                str_replace_all(.y, "'|-", "_"),
                                ".png",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = "png",
               height = 6,
               width = 18,
               units = "in",
               dpi = "retina")
)

walk2(.x = iclip_maps_1se,
      .y = names(iclip_maps_1se),
      ~ ggsave(filename = paste("2023-11-28_background_shsy5y_papa_cryptic_iclip_map.horiz_stack.fixed_ylim.1_se_ci.",
                                str_replace_all(.y, "'|-", "_"),
                                ".svg",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = svg,
               height = 6,
               width = 18,
               units = "in",
               dpi = "retina")
)

# write plot dfs to file
walk2(.x = iclip_dfs_maps_1se,
      .y = names(iclip_dfs_maps_1se),
      ~ write_tsv(.x, 
                  paste("processed/iclip_maps/plots/2023-11-28_background_shsy5y_papa_cryptic_iclip_df.horiz_stack.fixed_ylim.1_se_ci.",
                            str_replace_all(.y, "'|-", "_"),
                            ".tsv",
                            sep = ""),
                  col_names = T
                  )
      )


# write single panel d3'utrs to file
# walk2(.x = iclip_maps_single_d3utr_1se,
#       .y = names(iclip_maps_single_d3utr_1se),
#       ~ ggsave(filename = paste("2023-10-03_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.",
#                                 .y,
#                                 ".svg",
#                                 sep = ""),
#                plot = .x,
#                path = "processed/iclip_maps/plots/",
#                device = svg,
#                height = 7.5,
#                width = 22.5,
#                units = "in",
#                dpi = "retina")
#       )
# 
# walk2(.x = iclip_maps_single_d3utr_1se,
#       .y = names(iclip_maps_single_d3utr_1se),
#       ~ ggsave(filename = paste("2023-10-03_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.",
#                                 .y,
#                                 ".png",
#                                 sep = ""),
#                plot = .x,
#                path = "processed/iclip_maps/plots/",
#                device = "png",
#                height = 7.5,
#                width = 22.5,
#                units = "in",
#                dpi = "retina")
# )

# write single panel le starts to file
walk2(.x = iclip_maps_single_le_start_1se,
      .y = names(iclip_maps_single_le_start_1se),
      ~ ggsave(filename = paste("2023-11-28_background_shsy5y_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.le_start.",
                                str_replace_all(.y, "'|-", "_"),
                                ".svg",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = svg,
               height = 7.5,
               width = 22.5,
               units = "in",
               dpi = "retina")
        )

walk2(.x = iclip_maps_single_pas_1se,
      .y = names(iclip_maps_single_pas_1se),
      ~ ggsave(filename = paste("2023-11-28_background_shsy5y_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.pas.",
                                str_replace_all(.y, "'|-", "_"),
                                ".svg",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = svg,
               height = 7.5,
               width = 22.5,
               units = "in",
               dpi = "retina")
)


walk2(.x = iclip_maps_single_le_start_1se,
      .y = names(iclip_maps_single_le_start_1se),
      ~ ggsave(filename = paste("2023-11-28_background_shsy5y_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.le_start.",
                                str_replace_all(.y, "'|-", "_"),
                                ".png",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = "png",
               height = 7.5,
               width = 22.5,
               units = "in",
               dpi = "retina")
)

walk2(.x = iclip_maps_single_pas_1se,
      .y = names(iclip_maps_single_pas_1se),
      ~ ggsave(filename = paste("2023-11-28_background_shsy5y_papa_cryptic_iclip_map.single_panel.fixed_ylim.1_se_ci.pas.",
                                str_replace_all(.y, "'|-", "_"),
                                ".png",
                                sep = ""),
               plot = .x,
               path = "processed/iclip_maps/plots/",
               device = "png",
               height = 7.5,
               width = 22.5,
               units = "in",
               dpi = "retina")
)
