library(tidyverse)


dbrn_paths <- list.files(path = "data/peka_qapa",
           pattern = "_6mer_distribution_genome.tsv$",
           recursive = T,
           full.names = T) %>%
  set_names(str_remove(basename(.), "_6mer_distribution_genome.tsv$"))

l_dbrn_tbls <- map(dbrn_paths,
                   ~ read_tsv(.x, n_max = 100) %>% rename(kmer = `...1`)
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

shsy5y_dbrn_tbl %>%
  filter(kmer == "AACGAG") %>%
  pivot_longer(cols = all_of(coord_cols),
             names_to = "rel_posn",
             values_to = "rel_occur",
             names_transform = list(rel_posn = as.integer)
             ) %>% 
  view()

# read in table with counts for number of PAS pairs in each category