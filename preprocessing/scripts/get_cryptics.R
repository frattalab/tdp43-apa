library(tidyverse)


df <- read_tsv("processed/PAPA/2023-12-10_i3_cortical_zanovello.all_datasets.dexseq_apa.results.processed.cleaned.tsv")


#
cryptics <- df %>%
  filter(padj < 0.05 & mean_PPAU_base < 0.1 & delta_PPAU_treatment_control > 0.1)


# summary info fo le_ids (aggregated across experiments)
cryptics_summ <- cryptics %>%
  group_by(le_id) %>%
  summarise(gene_name = unique(gene_name),
            simple_event_type = paste(sort(unique(simple_event_type)), collapse = ","),
            n_exper_cryptic = n_distinct(experiment_name),
            experiment_name = paste(sort(unique(experiment_name)), collapse = ","),
            annot_status = paste(sort(unique(annot_status)), collapse = ","),
            chromosome = unique(chromosome),
            start = paste(sort(unique(start)), collapse = ","),
            end = paste(sort(unique(end)), collapse = ","),
            strand = unique(strand)
            ) %>%
  arrange(desc(n_exper_cryptic), simple_event_type)

cryptics_summ

# want to check whether some 'spliced' events are in fact the proximal 3'UTR of a locus with a novel 3'UTR extension (and re-annotate if necessary)

# vector of le_ids that are novel 3'UTR extensions
d3utr_le_ids <- df %>%
  filter(simple_event_type == "distal_3utr_extension") %>%
  distinct(le_id) %>%
  pull()

# construct a 'putative le_id' for cryptic events - i.e. if it had a distal 3'UTR, this would be it's le_id
# gene_id prefix is same, but number suffix = le_id + 1
cryptics_summ <- cryptics_summ %>%
  mutate(distal_le_id = paste(str_remove_all(le_id, "_\\d+$"),
                              as.character(as.integer(str_split_i(le_id, "_", 2)) + 1),
                              sep = "_"))

# now check if putative distal_le_id is in fact a novel 3'UTR extension
# if so, update the simple event type annotation to proximal_3utr_extension
cryptics_summ <- cryptics_summ %>%
  mutate(has_distal_3utr = simple_event_type == "spliced" & distal_le_id %in% d3utr_le_ids) %>%
  # updare
  mutate(simple_event_type = if_else(has_distal_3utr,
                                     "proximal_3utr_extension", simple_event_type)) %>%
  # remove intermediate columns
  select(-distal_le_id, -has_distal_3utr)
  


cryptic_et_exper <- count(cryptics, simple_event_type, experiment_name) %>%
  arrange(simple_event_type, desc(n))

write_tsv(cryptics_summ, "processed/2023-12-10_cryptics_summary_all_events.tsv",col_names = T)

# annotate multiple event types as complex
cryptics_summ_c <- cryptics_summ %>%
  mutate(simple_event_type = if_else(str_detect(simple_event_type, ","),
                                     "complex",
                                     simple_event_type)) 

cryptics_summ_c %>%
  count(simple_event_type)
# simple_event_type           n
# <chr>                   <int>
#   1 bleedthrough               51
# 2 complex                     9
# 3 distal_3utr_extension      86
# 4 proximal_3utr_extension    20
# 5 spliced                    92

cryptics_summ_c %>%
  count(simple_event_type, annot_status) %>%
  arrange(simple_event_type)
# A tibble: 11 × 3
# simple_event_type       annot_status        n
# <chr>                   <chr>           <int>
#   1 bleedthrough            annotated          16
# 2 bleedthrough            novel              35
# 3 complex                 annotated           3
# 4 complex                 annotated,novel     5
# 5 complex                 novel               1
# 6 distal_3utr_extension   novel              86
# 7 proximal_3utr_extension annotated          17
# 8 proximal_3utr_extension annotated,novel     3
# 9 spliced                 annotated          25
# 10 spliced                 annotated,novel    10
# 11 spliced                 novel              57

# write full tsv with complex event type annotation
write_tsv(cryptics_summ_c, "processed/2023-12-10_cryptics_summary_all_events_complex.tsv",col_names = T)

# for each event type, write to separate tsv
for (event_type in unique(cryptics_summ_c$simple_event_type)) {
  
  write_tsv(filter(cryptics_summ_c, simple_event_type == event_type),
            paste("processed/2023-12-10_cryptics_summary_", event_type, ".tsv", sep = ""),
            col_names = T
            )
  
}

# write cryptic filtered PAPA TSV
write_tsv(cryptics,
          "processed/PAPA/2023-12-10_cryptics.i3_cortical_zanovello.all_datasets.dexseq_apa.results.processed.cleaned.tsv")




