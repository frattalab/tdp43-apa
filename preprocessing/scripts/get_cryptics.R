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

cryptic_et_exper <- count(cryptics, simple_event_type, experiment_name) %>%
  arrange(simple_event_type, desc(n))

write_tsv(cryptics_summ, "processed/2023-12-10_cryptics_summary_all_events.tsv",col_names = T)

# annotate multiple event types as complex
cryptics_summ_c <- cryptics_summ %>%
  mutate(simple_event_type = if_else(str_detect(simple_event_type, ","),
                                     "complex",
                                     simple_event_type)) 

# count number of events in each category
cryptics_summ_c_counts <- cryptics_summ_c %>%
  count(simple_event_type)

# A tibble: 4 × 2
# simple_event_type         n
# <chr>                 <int>
#   1 bleedthrough             51
# 2 complex                   9
# 3 distal_3utr_extension    86
# 4 spliced                 112

cryptics_summ_c %>%
  count(simple_event_type, annot_status) %>%
  arrange(simple_event_type)
# A tibble: 9 × 3
# simple_event_type     annot_status        n
# <chr>                 <chr>           <int>
#   1 bleedthrough          annotated          18
# 2 bleedthrough          novel              37
# 3 complex               annotated           4
# 4 complex               annotated,novel     6
# 5 complex               novel               5
# 6 distal_3utr_extension novel             104
# 7 spliced               annotated          48
# 8 spliced               annotated,novel    13
# 9 spliced               novel              58

# write full tsv with complex event type annotation
write_tsv(cryptics_summ_c, "processed/2023-12-10_cryptics_summary_all_events_complex.tsv",col_names = T)
# write tsv with counts for each event type
write_tsv(cryptics_summ_c_counts, "processed/2023-12-10_cryptics_summary_all_events_complex_counts.tsv", col_names = T)

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




