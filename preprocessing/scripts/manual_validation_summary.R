library(tidyverse)


# Create a final table of manually validated cryptic bleedthrough events

i3_mv <- read_tsv("data/PAPA/riboseq_manual_verification_of_i3_cortical_cryptic_bleedthroughs.tsv")
other_mv <- read_tsv("data/PAPA/cryptics_summary_bleedthrough_manual_validation.tsv")
complex_mv <- read_tsv("data/PAPA/cryptics_summary_complex_manual_curation.tsv")


# combine two dfs
joined_mv <- other_mv %>%
  select(-all_of(c("chromosome", "start", "end", "strand"))) %>%
  left_join(select(i3_mv,
                   -all_of(c("chromosome", "start", "end", "strand", "exper_cryp", "...12", "...13", "gene_name"))),
            by = "le_id", suffix = c(".all",".i3")) 

# geerate common manual validation and note columns for each source
joined_mv_clean <- joined_mv %>%
  mutate(event_manual_validation_f = if_else(is.na(event_manual_validation.all),
                                             event_manual_validation.i3,
                                             event_manual_validation.all),
         notes_f = if_else(is.na(notes.all),
                          notes.i3,
                           notes.all)
         ) %>% 
  select(-matches("\\.all$|\\.i3$")) %>%
  rename(event_manual_validation = event_manual_validation_f,
         notes = notes_f) %>%
  relocate(event_manual_validation, .before = cryptic_riboseq_reads)


mv_summary <- joined_mv_clean %>%
  count(event_manual_validation,sort = T) %>%
  mutate(frac = n / sum(n))

# A tibble: 5 × 3
# event_manual_validation     n   frac
# <chr>                   <int>  <dbl>
#   1 no                         29 0.527 
# 2 yes                        21 0.382 
# 3 ?                           2 0.0364
# 4 yes?                        2 0.0364
# 5 no?                         1 0.0182

write_tsv(joined_mv_clean, "processed/bleedthrough_manual_validation.tsv", col_names = T)
write_tsv(mv_summary, "processed/bleedthrough_summary_counts_manual_validation.tsv", col_names = T)


## Filter cryptics summary table for events not failing manual validation

cryp_summ <- read_tsv("processed/2023-12-10_cryptics_summary_all_events.tsv")

fail_ids <- joined_mv_clean %>% filter(event_manual_validation != "yes") %>% pull(le_id)

cryp_summ_f <- filter(cryp_summ, !(le_id %in% fail_ids))

write_tsv(cryp_summ_f, "processed/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv")
write_lines(fail_ids, "processed/2023-12-10_cryptics_manual_validation_fail_le_ids.txt")

# annotate multiple event types as complex
cryp_summ_f_c <- cryp_summ_f %>%
  mutate(simple_event_type = if_else(str_detect(simple_event_type, ","),
                                     "complex",
                                     simple_event_type)) 

write_tsv(cryp_summ_f, "processed/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation_complex.tsv")

# number of cryptic events
n_distinct(cryp_summ_f$le_id)
# [1] 227

# number for each category
cryp_summ_f %>%
  count(simple_event_type, sort = T) %>%
  mutate(frac = n / sum(n)) %>%
  write_tsv("processed/2023-12-10_cryptics_summary_counts_bleedthrough_manual_validation.tsv")

cryp_summ_f_c %>%
  count(simple_event_type, sort = T) %>%
  mutate(frac = n / sum(n)) %>%
  write_tsv("processed/2023-12-10_cryptics_summary_counts_complex_bleedthrough_manual_validation.tsv")

