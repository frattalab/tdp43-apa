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
      filter(cryptic_status == "cr")
    )

# extract le_id from the 'name' field, pivot to wider format with column per region (and 1/0 for binding)
summary_dfs <- map(summary_dfs,
                   ~ .x %>%
                     separate(name, into = c("le_id", "gene_name", NA, NA), sep = "\\|", remove = F) %>%
                     pivot_wider(id_cols = all_of(c("le_id", "gene_name")), names_from = region_type, values_from = binding,names_prefix = "binding_")
                     )
# combine into a single df and re-annotate event types with simple values
comb_summary <- summary_dfs %>%
  bind_rows(.id = "simple_event_type") %>%
  mutate(event_type = case_when(simple_event_type == "spliced" ~ "ALE",
                                simple_event_type == "bleedthrough_uniq" ~ "IPA",
                                simple_event_type == "d3utr" ~ "3'Ext",
                                simple_event_type == "d3utr_proximal" ~ "3'Shortening",
                                TRUE ~ "NA")) %>%
  select(-simple_event_type) %>%
  relocate(event_type, .after = gene_name) 

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

write_tsv(comb_summary, file.path(outdir, "2024-11-26_cryptic_iclip_summary_cleaned_combined.tsv"))
