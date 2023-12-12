library(tidyverse)

# Because had accidentally kept events cryptic in only intermediate curve KDs, need to remove from BED to get accurate junction count


jncs <- read_tsv("processed/2023-09-12_papa_as_ale_cryptic.le.end_shift.junctions.bed",
                 col_names = c("chrom", "start", "end", "name", "score", "strand")
                 )

# pull out le_id from name
# e.g. ENSG00000162390.18_4|ACOT11|spliced|cryptic -> ENSG00000162390.18_4
jncs <- mutate(jncs, le_id = str_split_i(name, "\\|", 1))


# get file of summarise cryptic events and categories AFTER removing intermediate 
cryptics_df <- read_tsv("../preprocessing/processed/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv")

# pull out ALE le_ids
spliced_le_ids <- cryptics_df %>%
  filter(simple_event_type == "spliced") %>%
  pull(le_id)

# remove any events only expressed in intermediate curve KDs
jncs_f <- jncs %>%
  filter(le_id %in% spliced_le_ids)

num_jnc <- distinct(jncs_f, .keep_all = T) %>%
  nrow()
# [1] 118

# 
num_ids <- distinct(jncs_f, le_id, .keep_all = T) %>%
  nrow()
# [1] 92

enframe(c("num_jncs" = num_jnc,
          "num_ids" = num_ids), "type", "n") %>%
  write_tsv("processed/2023-12-12_as_ale_cryptic_jncs_num_unique_events.tsv")
