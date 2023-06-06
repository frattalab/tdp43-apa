library(tidyverse)
library(glue)
library(zoo)

#' convert peka _distribution table to long format with one row per position and kmer
peka_wide_to_long <- function(df, kmers, first_posn_idx = 14) {
  
  # all remaining columns after first are the position cols
  coord_cols <- colnames(df)[14:length(colnames(df))]

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

#' Construct a general kmer occurrence plot relative to PAS (with optional smoothing by rolling mean)
plot_kmer_dbrn <- function(plot_df, rolling_mean = F, rolling_k = 5, facet_w = "~ label", title = "", subtitle = "", n_row = 2, n_col = 1, base_size = 14) {
  

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
  scale_x_continuous(breaks = seq(-150, 150, 20)) +
  scale_colour_manual(values = c("#1b9e77", "#7570b3")) +
  theme_classic(base_size = base_size) +
  labs(title = title,
       subtitle = subtitle,
       x = "Distance from polyA site (nt)",
       y = "Relative kmer occurrence")

}

dbrn_paths <- list.files(path = "data/peka_qapa",
           pattern = "_6mer_distribution_genome.tsv$",
           recursive = T,
           full.names = T) %>%
  set_names(str_remove(basename(.), "_6mer_distribution_genome.tsv$"))

l_dbrn_tbls <- map(dbrn_paths,
                   ~ read_tsv(.x, show_col_types = F) %>% rename(kmer = `...1`)
                   )

dbrn_tbl <- bind_rows(l_dbrn_tbls, .id = "experiment_comparison_name")
# split experiment_comparison_name into sub strings for easier categorisation
dbrn_tbl <- separate(dbrn_tbl,
         experiment_comparison_name,
         into = c("experiment_name", "contrast", "pas", NA, "direction"),
         sep = "\\.",
         remove = F
         )


shsy5y_dbrn_tbl <- filter(dbrn_tbl, str_detect(experiment_name, "shsy5y"))

# read in table with counts for number of PAS pairs in each category & experiment
# Will use to provide more informative legend names
npas_tbl <- read_tsv("data/peka_qapa/n_pas_groups_by_window.tsv")
# here only used PAS separated by at least 302 nt (extend 150 either side of PAS and no overlap)
npas_150_tbl <- filter(npas_tbl, window_size == "150nt", group %in% c("regulated", "background"))
# easier joining later
npas_150_tbl <- separate(npas_150_tbl, 
                         experiment_name_contrast,
                         into = c("experiment_name", "contrast"),
                         sep = "\\.")
shsy5y_npas_150_tbl <- filter(npas_150_tbl, str_detect(experiment_name, "shsy5y"))

# try a plot of just canonical TDP-43 motifs
shsy5y_eg_long <- peka_wide_to_long(shsy5y_dbrn_tbl, c("GUGUGU", "UGUGUG"))
shsy5y_plot_labs <- npas_plot_label(shsy5y_npas_150_tbl)
shsy5y_eg_long_plot_df <- df_plot_kmer_dbrn(shsy5y_eg_long, filter(shsy5y_plot_labs, group == "regulated"))

plot_kmer_dbrn(shsy5y_eg_long_plot_df, facet_w = "experiment_name ~ label", n_row = 4, n_col = 1) +
  labs(title = "kmers = GUGUGU & UGUGUG")

plot_kmer_dbrn(shsy5y_eg_long_plot_df,rolling_mean = T, rolling_k = 5, facet_w = "experiment_name ~ label",
               title = "kmers = GUGUGU & UGUGUG",
               subtitle = "smoothed with rolling mean (k = 5)",
               n_row = 4, n_col = 1)


plot_kmer_dbrn(shsy5y_eg_long_plot_df, facet_w = "experiment_name ~ label",
               title = "kmers = GUGUGU & UGUGUG", n_row = 2, n_col = 2)

plot_kmer_dbrn(shsy5y_eg_long_plot_df, rolling_mean = T, rolling_k = 5, facet_w = "experiment_name ~ label",               
               title = "kmers = GUGUGU & UGUGUG",
               subtitle = "smoothed with rolling mean (k = 5)",n_row = 2, n_col = 2)

all_eg_long <- peka_wide_to_long(dbrn_tbl, c("GUGUGU", "UGUGUG"))
all_plot_labs <- npas_plot_label(npas_150_tbl)
all_eg_long_plot_df <- df_plot_kmer_dbrn(all_eg_long, filter(all_plot_labs, group == "regulated"))

# no smoothing all datasets
plot_kmer_dbrn(all_eg_long_plot_df, facet_w = "experiment_name ~ label", n_row = 5, n_col = 2) +
  labs(title = "kmers = GUGUGU & UGUGUG")

# with smoothing all datasets
all_canon_plot <- plot_kmer_dbrn(all_eg_long_plot_df,
               rolling_mean = T,
               rolling_k = 5,
               facet_w = "experiment_name ~ label",
               title = "kmers = GUGUGU & UGUGUG",
               subtitle = "smoothed with rolling mean (k = 5)",
               n_row = 5,
               n_col = 2) +
  scale_y_continuous(labels = abs(seq(-4, 4, 1)),
                     breaks = seq(-4, 4, 1)
  )

all_canon_plot

if (!dir.exists("processed")) {dir.create("processed")}

ggsave(plot = all_canon_plot,
       filename = "2023-06-06_qapa_sig_peka_6mer_gu_ug_dibrn_line_plot.png",
       path = "processed",
       device = "png",
       height = 8,
       width = 12,
       units = "in",
       dpi = "retina")

