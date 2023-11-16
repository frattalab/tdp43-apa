library(tidyverse)

# find paths to all kmer distrivution tables for different comparisons
dbrn_paths <- list.files(path = "data/peka_papa/2023-11-16_papa_cryptics_kmer6_window_250_distal_window_500_relpos_0",
                         pattern = "_UGUGUG_GUGUGU.tsv$",
                         recursive = T,
                         full.names = T)
dbrn_paths

# get the 'comparison name' from the name of bed file
comparison_names <- str_remove(basename(dbrn_paths), ".foreground_genome_20_UGUGUG_GUGUGU.tsv$")
# want to go two folders back to whether foreground or background seqs (cv_coverage_<foreground/background>)
# then extract whether foreground/background (final field after split by _)
region_names <- str_split_i(basename(dirname(dirname(dbrn_paths))), "_", 3)
# combine to construct final ID
comparison_ids <- paste(comparison_names, region_names, sep = ";")

# read in tables, combine into single df
dbrn_tbl <- dbrn_paths %>%
  set_names(comparison_ids) %>%
  map(~ read_tsv(.x,
                 col_names = c("rel_posn", "rel_occur"),
                 skip = 1, # skip defined header
                 show_col_types = F
                 )
      ) %>%
  bind_rows(.id = "comparison_id")





# construct comparison name - first element before ;
# region type - exonstart / pas / pas_proximal / pas_distal
# group = foreground/background
dbrn_tbl <- dbrn_tbl %>%
  separate(comparison_id, into = c("comparison_name", "group"), sep = ";", remove = F) %>%
  mutate(region_type = str_extract(comparison_name,  paste(c("pas$","exonstart", "pas_proximal", "pas_distal"), collapse = "|")),
         .after = comparison_name)

if (!dir.exists("processed/peka/papa")) {dir.create("processed/peka/papa", recursive = T)}

write_tsv(dbrn_tbl,
          "processed/peka/papa/2023-11-16_papa_cryptics_cvcoverage_window_500.gugugu_ugugug.distribution_genome.tsv",
          col_names = T)