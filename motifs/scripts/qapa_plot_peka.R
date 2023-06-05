library(tidyverse)
library(glue)

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


# standard results cols are first 13 - want to pivot position columns (remaining cols)
coord_cols <- (colnames(shsy5y_dbrn_tbl)[14:length(shsy5y_dbrn_tbl)])

shsy5y_eg_long <- shsy5y_dbrn_tbl %>%
  filter(kmer %in% c("GUGUGU", "UGUGUG")) %>%
  pivot_longer(cols = all_of(coord_cols),
             names_to = "rel_posn",
             values_to = "rel_occur",
             names_transform = list(rel_posn = as.integer)
             )

# read in table with counts for number of PAS pairs in each category & experiment
# Will use to provide more informative legend names
npas_tbl <- read_tsv("data/peka_qapa/n_pas_groups_by_window.tsv")
# here only used PAS separated by at least 302 nt (extend 150 either side of PAS and no overlap)
npas_150_tbl <- filter(npas_tbl, window_size == "150nt", group %in% c("regulated", "background"))
shsy5y_npas_150_tbl <- filter(npas_150_tbl, str_detect(experiment_name_contrast, "shsy5y"))

npas_150_tbl

# create title labels for different plots 
# Experiment name - contrast - <proximal/distal> - up in KD n = x, down in KD n = x
shsy5y_plot_labs <- shsy5y_npas_150_tbl %>%
  # separate for joining later
  separate(experiment_name_contrast,into = c("experiment_name", "contrast"), sep = "\\.") %>%
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


shsy5y_eg_long_plot_df <- left_join(shsy5y_eg_long, filter(shsy5y_plot_labs, group == "regulated"), by = c("experiment_name", "contrast", "pas")) 


shsy5y_eg_long_plot_df %>%
  # sort now so can create a factor for plotting with proximal sites first in legend
  arrange(experiment_name, contrast, desc(pas)) %>%
  # want to plot up and down on the same graph, so need to artificially make down negative so can plot on same axes
  mutate(rel_occur = if_else(direction == "down",
                             -1 * rel_occur,
                             rel_occur),
         pas = factor(pas, levels = c("proximal", "distal")),
         label = fct_inorder(label),
         direction = factor(direction, levels = c("up", "down"))
         ) %>%
  # filter(experiment_name == "zanovello_shsy5y_chx") %>%
  ggplot(aes(x = rel_posn,
             y = rel_occur,
             colour = direction,
             linetype = kmer)) +
  facet_wrap("experiment_name ~ label",
             scales = "free_y",
             ncol = 1) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-150, 150, 20)) +
  scale_y_continuous(labels = abs(seq(-8, 8, 1)),
                     breaks = seq(-8, 8, 1)
                     ) +
  scale_colour_manual(values = c("#1b9e77", "#7570b3")) +
  theme_classic(base_size = 14) +
  labs(title = "kmers = GUGUGU & UGUGUG",
       x = "Distance from polyA site (nt)",
       y = "Relative kmer occurrence")
  
shsy5y_eg_long <- shsy5y_dbrn_tbl %>%
  filter(kmer %in% c("CGCGGA", "CGCGUG")) %>%
  pivot_longer(cols = all_of(coord_cols),
               names_to = "rel_posn",
               values_to = "rel_occur",
               names_transform = list(rel_posn = as.integer)
  )


# create title labels for different plots 
# Experiment name - contrast - <proximal/distal> - up in KD n = x, down in KD n = x
shsy5y_plot_labs <- shsy5y_npas_150_tbl %>%
  # separate for joining later
  separate(experiment_name_contrast,into = c("experiment_name", "contrast"), sep = "\\.") %>%
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


shsy5y_eg_long_plot_df <- left_join(shsy5y_eg_long, filter(shsy5y_plot_labs, group == "regulated"), by = c("experiment_name", "contrast", "pas")) 


shsy5y_eg_long_plot_df %>%
  # sort now so can create a factor for plotting with proximal sites first in legend
  arrange(experiment_name, contrast, desc(pas)) %>%
  # want to plot up and down on the same graph, so need to artificially make down negative so can plot on same axes
  mutate(rel_occur = if_else(direction == "down",
                             -1 * rel_occur,
                             rel_occur),
         pas = factor(pas, levels = c("proximal", "distal")),
         label = fct_inorder(label),
         direction = factor(direction, levels = c("up", "down"))
  ) %>%
  # filter(experiment_name == "zanovello_shsy5y_chx") %>%
  ggplot(aes(x = rel_posn,
             y = rel_occur,
             colour = direction,
             linetype = kmer)) +
  facet_wrap("experiment_name ~ label",
             scales = "free_y",
             ncol = 1) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-150, 150, 20)) +
  scale_y_continuous(labels = abs(seq(-8, 8, 1)),
                     breaks = seq(-8, 8, 1)
  ) +
  scale_colour_manual(values = c("#1b9e77", "#7570b3")) +
  theme_classic(base_size = 14) +
  labs(title = "kmers = CGCGGA & CGCGUG",
       x = "Distance from polyA site (nt)",
       y = "Relative kmer occurrence")

