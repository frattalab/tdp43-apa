library(tidyverse)

outdir <- "processed"

# find paths to per-event binding summaries
summary_paths <- list.files(path = "data",
                            pattern = "^2024-11-13_event\\.",full.names = T)

# extract event type from file name
summary_paths <- set_names(summary_paths, str_split_i(basename(summary_paths), "\\.", 4))

# read in binding info for cryptics
summary_dfs <- map(summary_paths,
    ~ read_tsv(.x, show_col_types = F) %>%
      filter(cryptic_status == "cr") %>%
      # extract le_id & gene name from the 'name' field
      separate(name, into = c("le_id", "gene_name", NA, NA), sep = "\\|", remove = F) 
    )

# pivot to wider format with column per region (and 1/0 for binding)
summary_dfs_wide <- imap(summary_dfs,
                         ~ if (.y %in% c("d3utr", "d3utr_proximal")) {
                           # have prox and dist le_ids represented per gene, so only use gene_name as id_col to get the 
                           .x %>%
                             pivot_wider(id_cols = gene_name, names_from = region_type, values_from = binding,names_prefix = "binding_")
                         } else {
                           .x %>%
                             pivot_wider(id_cols =  all_of(c("le_id", "gene_name")), names_from = region_type, values_from = binding,names_prefix = "binding_")
                         }
)

comb_summary <- summary_dfs_wide %>%
  bind_rows(.id = "simple_event_type")

# now need to unify le_ids for 3'UTR changes

# Extract a df of le_id | gene_name for 3'Ext and 3'Shortening events
id2name_d3utrs <- keep_at(summary_dfs, c("d3utr", "d3utr_proximal")) %>%
  imap(~ if (.y == "d3utr") {
    # select the distal event as the cryptic id
    .x %>%
      filter(str_detect(name, "distal")) %>%
      distinct(le_id, gene_name)
  } else {
    # extract proximal event as the cryptic id
    .x %>%
      filter(str_detect(name, "proximal")) %>%
      distinct(le_id, gene_name)
    }
  ) %>%
  bind_rows()

# add le_ids to summary df and unify le_id based on event type
comb_summary <- comb_summary %>%
  left_join(id2name_d3utrs, by = "gene_name", suffix = c("", ".d3utr")) %>%
  mutate(le_id = if_else(simple_event_type %in% c("d3utr", "d3utr_proximal"), le_id.d3utr, le_id)
         ) %>%
  select(-le_id.d3utr)

# Re-annotate event types with clean names
comb_summary <- comb_summary %>%
  mutate(event_type = case_when(simple_event_type == "spliced" ~ "ALE",
                                simple_event_type == "bleedthrough_uniq" ~ "IPA",
                                simple_event_type == "d3utr" ~ "3'Ext",
                                simple_event_type == "d3utr_proximal" ~ "3'Shortening",
                                TRUE ~ "NA")) %>%
  select(-simple_event_type) %>%
  relocate(event_type, .after = gene_name) %>%
  arrange(event_type, desc(binding_start), desc(binding_end))

comb_summary %>%
  a

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

write_tsv(comb_summary, file.path(outdir, "2024-11-26_cryptic_iclip_summary_cleaned_combined.tsv"))
