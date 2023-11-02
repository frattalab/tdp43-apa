library(tidyverse)
library(glue)
library(zoo)

#' convert peka _distribution table to long format with one row per position and kmer
peka_wide_to_long <- function(df, kmers, first_posn_idx = 14) {
  
  # all remaining columns after first are the position cols
  coord_cols <- colnames(df)[first_posn_idx:length(colnames(df))]
  
  df %>%
    filter(kmer %in% kmers) %>%
    pivot_longer(cols = all_of(coord_cols),
                 names_to = "rel_posn",
                 values_to = "rel_occur",
                 names_transform = list(rel_posn = as.integer)
    )    
  
}

#' convert n_pas_groups_by_window.tsv df to labels of counts in each direction for proximal and distal
npas_plot_label <- function(df) {
  
  df %>%
    # counts rows to columns for easier parsing
    pivot_wider(names_from = dirn, values_from = n_pas, names_prefix = "prox_") %>%
    mutate(proximal = glue("Proximal - up in contrast n = {prox_up_in_kd}, down in contrast n = {prox_down_in_kd}"),
           distal = glue("Distal - up in contrast n = {prox_down_in_kd}, down in contrast n = {prox_up_in_kd}")
    ) %>%
    # single rows per label for site type
    pivot_longer(cols = all_of(c("proximal", "distal")),
                 names_to = "pas",
                 values_to = "label") %>%
    select(experiment_name, contrast, group, pas, label)
  
}

df_plot_kmer_dbrn <- function(peka_df, count_df, join_cols = c("experiment_name", "contrast", "pas")) {
  
  joined <- left_join(peka_df, count_df, by = join_cols )
  
  joined %>%
    # sort now so can create a factor for plotting with proximal sites first in legend
    arrange(experiment_name, contrast, desc(pas)) %>%
    # want to plot up and down on the same graph, so need to artificially make down negative so can plot on same axes
    mutate(rel_occur = if_else(direction == "down",
                               -1 * rel_occur,
                               rel_occur),
           pas = factor(pas, levels = c("proximal", "distal")),
           label = fct_inorder(label),
           direction = factor(direction, levels = c("up", "down"))
    ) 
}

df_plot_pair_kmer_dbrn <- function(peka_df, count_df, join_cols = c("experiment_name", "contrast", "pas")) {
  
  joined <- left_join(peka_df, count_df, by = join_cols )
  
  joined %>%
    group_by(experiment_name) %>%
    mutate(pair = if_else((pas == "proximal" & direction == "up") | (pas == "distal" & direction == "down"),
                          "proximal_up_distal_down",
                          "proximal_down_distal_up"),
           pair = factor(pair, levels = c("proximal_up_distal_down", "proximal_down_distal_up"))
    ) %>%
    # sort now so can create a factor for plotting with proximal sites first in legend
    arrange(experiment_name, contrast, desc(pas)) %>%
    # want to plot up and down on the same graph, so need to artificially make down negative so can plot on same axes
    mutate(rel_occur = if_else(direction == "down",
                               -1 * rel_occur,
                               rel_occur),
           pas = factor(pas, levels = c("proximal", "distal")),
           label = fct_inorder(label),
           direction = factor(direction, levels = c("up", "down"))
    ) 
}

#' Construct a general kmer occurrence plot relative to PAS (with optional smoothing by rolling mean)
plot_kmer_dbrn <- function(plot_df, rolling_mean = F, rolling_k = 5, facet_w = "~ label", title = "", subtitle = "", n_row = 2, n_col = 1, base_size = 14, breaks = seq(-150, 150, 20)) {
  
  
  plot_base <- plot_df %>%
    ggplot(aes(x = rel_posn,
               y = rel_occur,
               colour = direction,
    )) +
    facet_wrap(facet_w,
               scales = "free_y",
               nrow = n_row,
               ncol = n_col) 
  
  if (rolling_mean) {
    plot_base <- plot_base + geom_line(aes(y = zoo::rollmean(rel_occur, k = rolling_k, na.pad = T)))
    
  } else {
    
    plot_base <- plot_base + geom_line()
    
  }
  
  plot_base +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(breaks = breaks) +
    scale_colour_manual(values = c("#1b9e77", "#7570b3")) +
    theme_classic(base_size = base_size) +
    labs(title = title,
         subtitle = subtitle,
         x = "Distance from polyA site (nt)",
         y = "Relative kmer occurrence")
  
}