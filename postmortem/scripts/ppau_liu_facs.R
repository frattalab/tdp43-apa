library(tidyverse)

#' helper function to collapse duplicated values in a delimited string to non-redundant values
#' I wasn't very smart with gene_name column in PAPA output - often get duplicated gene name values...
collapse_names <- function(col, split=",") {
  
  apply(str_split(col,  split, simplify = T), # generate a matrix of strings separated by split 
        MARGIN = 1, # loop over rows
        FUN = function(x) {paste(unique(x)[unique(x) != ""],
                                 collapse = split)
        }, # collapse split gene names to unique values, then recombine if there are multiple IDs remaining
        simplify = T)
}



# path of matrix storing per-sample % poly A usage for each last exon
p_ppau <- "processed/2023-06-21_summarised_pas_quantification_cryptics_fix_tx2le.ppau.tsv"

# path to df mapping input transcripts (i.e. cryptic transcripts) to their assigned last exon isoform
p_tx2le <- "processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.tx2le.tsv"

# path to df mapping txs to last exon isoforms, focused only on input cryptic last exons
# just in case new annotation changes le_id assignment, will use this to map back to cryptic exons called with standard approach
p_input_tx2le <- "processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.input_tx2le.tsv"

# use final definition of cryptic events
p_cryp <- "../preprocessing/processed/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv"

p_orig_tx2le <- "data/PAPA/novel_ref_combined.tx2le.tsv"

p_le2name <- "processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.le2name.tsv"

p_sample_tbl <- "data/liu_facs_papa_sample_sheet.csv"

# due to some ID mapping errors and overlapping event types, performed manual curation of event types 
# (some due to e.g. overlapping isoforms being grouped together - focused on event type of cryptic by expression region)
p_event_type_na_mc <- "data/liu_facs_event_type_manual_curation.tsv"
p_event_type_complex_mc <- "data/cryptics_summary_complex_manual_curation.tsv"


ppau <- read_tsv(p_ppau)
tx2le <- read_tsv(p_tx2le)
le2name <- read_tsv(p_le2name) %>%
  # wasn't smart with gene_name collapsing, often get duplicates...
  mutate(gene_name = collapse_names(gene_name))
sample_tbl <- read_csv(p_sample_tbl)

event_type_na_mc <- read_tsv(p_event_type_na_mc)
event_type_complex_mc <- read_tsv(p_event_type_complex_mc)


## 
cryp <- read_tsv(p_cryp)
orig_tx2le <- read_tsv(p_orig_tx2le)


# get tx_ids of cryptic le_ids
orig_tx_cryp <- orig_tx2le %>% 
  filter(le_id %in% cryp$le_id) %>%
  pull(transcript_id)

# get a vector of cryptic last exon IDs (for current annotation)
cryp_le_ids <- tx2le %>%
  filter(transcript_id %in% orig_tx_cryp) %>%
  distinct(le_id) %>%
  pull()

# extract event type annotation
le2event <- cryp %>%
  select(gene_name, le_id, simple_event_type)
  

# add plot_le_id
le2event <- mutate(le2event, plot_le_id = paste(gene_name, str_split_i(le_id, "_", 2), sep = "_"))


# have manual curation of event types
# combine the two tables, filtering for passing events
event_type_complex_mc <- event_type_complex_mc %>%
  mutate(plot_le_id = paste(gene_name, str_split_i(le_id, "_", 2), sep = "_"))

event_type_complex_mc_p <- event_type_complex_mc %>%
  select(le_id, plot_le_id, event_manual_validation, simple_event_type = manual_simple_event_type) %>%
  # some events did not pass manual curation - filter out
  filter(event_manual_validation == "yes")

event_type_na_mc_p <- event_type_na_mc %>%
  filter(event_manual_validation == "yes") %>%
  select(plot_le_id, simple_event_type = manual_simple_event_type, event_manual_validation)

event_type_mc <- bind_rows(event_type_na_mc_p, event_type_complex_mc_p)

le2event <- bind_rows("original" = le2event, "manual" = event_type_mc, .id = "source") %>%
  group_by(plot_le_id) %>%
  # filter(any(source == "manual")) %>%
  # if le_id contains two annotations, manual will come first (i.e. manual curation), so prefer over complex annotation
  # if only in original, will keep whatever annotation provided there
  arrange(source, .by_group = T) %>%
  slice_head(n = 1) %>% 
  ungroup()

# above aslo has purpose of adding missign le_ids to tabel (due to mismapping)
# need to double check that manually curated events are remvoed from le2event

event_type_mc_f_ids <- c(event_type_complex_mc %>% filter(event_manual_validation != "yes") %>% pull(plot_le_id),
                         event_type_na_mc %>% filter(event_manual_validation != "yes") %>% pull(plot_le_id)
)

le2event <- filter(le2event, !plot_le_id %in% event_type_mc_f_ids)

count(le2event, simple_event_type)
# A tibble: 5 × 2
# simple_event_type                      n
# <chr>                              <int>
#   1 bleedthrough                        2495
# 2 bleedthrough,distal_3utr_extension   118
# 3 bleedthrough,spliced                 986
# 4 distal_3utr_extension               3251
# 5 spliced                             6992


# join in gene name for each le_id
# dfine a more readable isoform id - gene_name suffixed with isoform number
# annotate whether cryptic/nto cryptic
ppau <- left_join(ppau, distinct(le2name, .keep_all = T), by = "le_id") %>% 
  mutate(plot_le_id = paste(gene_name, str_split_i(le_id, "_", 2), sep = "_")) %>%
  relocate(gene_name, plot_le_id, .after = gene_id) %>%
  # don't need these columns
  select(-starts_with("mean"), -starts_with("delta")) %>%
  mutate(cryptic_status = if_else(le_id %in% cryp_le_ids, T, F)) %>%
  left_join(select(le2event, plot_le_id, ev_le_id = le_id, simple_event_type, source),
            by = "plot_le_id") %>%
  relocate(simple_event_type, .after = plot_le_id) %>%
  # make sure remove le_ids if manually curated and present in ppau df
  filter(!(plot_le_id %in% event_type_mc_f_ids))




# join in metadata from sample table (converting to long format with row per sample and isoform)
ppau <- ppau %>%
  pivot_longer(cols = -all_of(c("le_id", "gene_id", "gene_name", "cryptic_status", "plot_le_id", "simple_event_type", "ev_le_id", "source")),
               names_to = "sample_name",
               values_to = "ppau") %>%
  left_join(select(sample_tbl, -path, -fastq1, -fastq2),
            by = "sample_name")


# now want to calculate a sample-wise difference in PPAU between TDP-43 negative & positive samples (condition column)
ppau_delta_paired <- ppau %>%
  select(-sample_name) %>%
  # for each sample & le_id get ppau values in positive and negative as columns
  pivot_wider(id_cols = c("le_id", "gene_id", "gene_name", "patient_id", "disease", "gender", "cryptic_status", "plot_le_id", "simple_event_type"),
              names_from = "condition",
              values_from = "ppau"
  ) %>%
  # calculate PPAU difference between -ve & +ve for each sample
  # (positive = higher usage in TDP-43 negative vs positive)
  mutate(paired_delta_ppau_neg_pos = TDPnegative - TDPpositive)


# calculate mean & median delta PAUs across samples for each isoform
ppau_delta_paired_median_all <- ppau_delta_paired %>%
  group_by(plot_le_id) %>%
  summarise(mean_paired_delta_ppau = mean(paired_delta_ppau_neg_pos),
            median_paired_delta_ppau = median(paired_delta_ppau_neg_pos)
  ) %>%
  ungroup() %>%
  left_join(distinct(ppau_delta_paired, le_id, gene_id, gene_name, cryptic_status, plot_le_id), by = "plot_le_id") %>%
  relocate(all_of(c("mean_paired_delta_ppau", "median_paired_delta_ppau"
  )),
  .after = last_col()
  ) %>%
  # put cryptics first
  arrange(desc(cryptic_status))


# calculate mean & median delta PAUs across samples for each isoform and disease subtype (FTD or ALS-FTD)
ppau_delta_paired_median_disease <- ppau_delta_paired %>%
  group_by(plot_le_id, disease) %>%
  summarise(mean_paired_delta_ppau = mean(paired_delta_ppau_neg_pos),
            median_paired_delta_ppau = median(paired_delta_ppau_neg_pos)
            ) %>%
  ungroup() %>%
  # return annotation info
  left_join(distinct(ppau_delta_paired, le_id, gene_id, gene_name, cryptic_status, plot_le_id), by = "plot_le_id") %>%
  relocate(all_of(c("mean_paired_delta_ppau", "median_paired_delta_ppau"
  )),
  .after = last_col()
  ) %>%
  # put cryptics first
  arrange(desc(cryptic_status))

if (!dir.exists("processed")) {dir.create("processed", recursive = T)}

# output median & mean pairwise delta_PPAU tables
ppau_delta_paired_median_all %>%
  filter(cryptic_status) %>%
  write_tsv("processed/2023-12-12_liu_facs_decoys_delta_ppau.all_samples.cryptics.tsv", col_names = T)

ppau_delta_paired_median_all %>%
  write_tsv("processed/2023-12-12_liu_facs_decoys_delta_ppau.all_samples.all_ales.tsv.gz", col_names = T)

ppau_delta_paired_median_disease %>%
  filter(cryptic_status) %>%
  write_tsv("processed/2023-12-12_liu_facs_decoys_delta_ppau.subtype_split.cryptics.tsv", col_names = T)

ppau_delta_paired_median_disease %>%
  write_tsv("processed/2023-12-12_liu_facs_decoys_delta_ppau.subtype_split.all_ales.tsv.gz", col_names = T)

# output per sample
ppau_delta_paired %>%
  filter(cryptic_status) %>%
  write_tsv("processed/2023-12-12_liu_facs_decoys_per_sample_delta_ppau.cryptics.tsv", col_names = T)

ppau_delta_paired %>%
  write_tsv("processed/2023-12-12_liu_facs_decoys_per_sample_delta_ppau.all_ales.tsv.gz", col_names = T)


  