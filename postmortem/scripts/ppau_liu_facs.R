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
p_ppau <- "processed/liu_facs/2024-11-20_summarised_pas_quantification_cryptics.ppau.tsv"

# path to df mapping input transcripts (i.e. cryptic transcripts) to their assigned last exon isoform
p_tx2le <- "processed/decoys/2024-11-20_decoys_novel_ref_combined.quant.tx2le.tsv"

p_le2name <- "processed/decoys/2024-11-20_decoys_novel_ref_combined.quant.le2name.tsv"

p_sample_tbl <- "data/liu_facs/liu_facs_papa_sample_sheet.tdppos_first.csv"

p_cryptics_summary <- "data/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv"

p_dexseq <- "data/PAPA/2023-05-24_i3_cortical_zanovello.all_datasets.dexseq_apa.results.processed.cleaned.tsv"

ppau <- read_tsv(p_ppau)
tx2le <- read_tsv(p_tx2le)
le2name <- read_tsv(p_le2name) %>%
  # wasn't smart with gene_name collapsing, often get duplicates...
  mutate(gene_name = collapse_names(gene_name))


sample_tbl <- read_csv(p_sample_tbl)
cryptics_summary <- read_tsv(p_cryptics_summary)
dexseq <- read_tsv(p_dexseq)

# extract event type annotation for cryptic events
# annotate a 'plot_le_id' - gene name if the sole cryptic event, otherwise appending the suffix
cryptics_summary_simple <- cryptics_summary %>%
  select(le_id, gene_name, simple_event_type) %>%
   separate(le_id, into = c("gene_id", "le_number"), sep = "_", remove = F) %>%
   group_by(gene_name) %>%
   mutate(gene_ncryp = n_distinct(le_id),
          plot_le_id = if_else(gene_ncryp == 1, gene_name, paste(gene_name, le_number, sep = "_"))) %>%
   ungroup() %>%
   select(-gene_ncryp, -le_number, -gene_id)

# pull out event ID annotations for all events
# Where events are cryptic, rely on the summary table (which has cleaned up/merged event types where mutliple)
le2event <- dexseq %>%
  distinct(le_id, gene_name, simple_event_type) %>%
  left_join(select(cryptics_summary_simple, le_id, plot_le_id, simple_event_type), by = "le_id", suffix = c("", "_cryp")) %>%
  mutate(simple_event_type_upd = if_else(is.na(simple_event_type_cryp), simple_event_type, simple_event_type_cryp)) %>%
  distinct(le_id, plot_le_id, gene_name, simple_event_type = simple_event_type_upd)


# add annotations to PPAU matrix - event types with all events, cryptic status (and whether evaluated by DEXSeq)
ppau_meta <- ppau %>%
    left_join(le2event, by = "le_id") %>%
    mutate(evaluated_status = if_else(is.na(gene_name) & is.na(simple_event_type),
                                      F,
                                      T), # not in DEXSeq table (i.e. never evaluated for diff usage)
           cryptic_status = if_else(is.na(plot_le_id), F, T) # column joined from cryptic annotation
           ) 

# pull out group-level mean & delta PPAUs for later querying
ppau_meta_means <- ppau_meta %>%
  select(-starts_with("FTD"))

# df with only per-sample PPAUs (and annotation cols)
ppau_meta_sample <-  ppau_meta %>%
  select(-starts_with("mean_"), -starts_with("delta"))

# Convert to long-format - rows per sample and isoform. Add in per-sample metadata
ppau_meta_sample_long <- ppau_meta_sample %>%
  pivot_longer(cols = -all_of(c("le_id", "gene_id", "gene_name", "plot_le_id", "simple_event_type", "evaluated_status", "cryptic_status")),
              names_to = "sample_name",
              values_to = "ppau") %>%
  left_join(select(sample_tbl,
                   -all_of(c("path", "fastq1", "fastq2"))
                   ), by = "sample_name")


# Calculate a sample-wise difference in PPAU between TDP-43 negative & positive samples (condition column)
ppau_delta_paired <- ppau_meta_sample_long %>%
  select(-sample_name) %>%
  # for each sample & le_id get ppau values in positive and negative as columns
  pivot_wider(id_cols = c("le_id", "gene_id", "gene_name", "patient_id", "disease", "gender", "cryptic_status", "evaluated_status", "plot_le_id", "simple_event_type"),
              names_from = "condition",
              values_from = "ppau"
  ) %>%
  # calculate PPAU difference between -ve & +ve for each sample
  # (positive = higher usage in TDP-43 negative vs positive)
  mutate(paired_delta_ppau_neg_pos = TDPnegative - TDPpositive)

# Calculate mean and median delta PPAU across samples of each isoform
ppau_delta_paired_median_all <- ppau_delta_paired %>%
  group_by(le_id) %>%
  summarise(mean_paired_delta_ppau = mean(paired_delta_ppau_neg_pos),
            median_paired_delta_ppau = median(paired_delta_ppau_neg_pos)
  ) %>%
  ungroup() %>%
  # return annotation cols
  left_join(distinct(ppau_delta_paired,
                     across(c("le_id", "gene_id", "gene_name", "evaluated_status", "cryptic_status", "plot_le_id", "simple_event_type"))
                     ),
            by = "le_id"
            ) %>%
  relocate(all_of(c("mean_paired_delta_ppau", "median_paired_delta_ppau"
  )),
  .after = last_col()
  ) %>%
  # put cryptics first
  arrange(desc(cryptic_status))

# calculate mean & median delta PAUs across samples for each isoform and disease subtype (FTD or ALS-FTD)
ppau_delta_paired_median_disease <- ppau_delta_paired %>%
  group_by(le_id, disease) %>%
  summarise(mean_paired_delta_ppau = mean(paired_delta_ppau_neg_pos),
            median_paired_delta_ppau = median(paired_delta_ppau_neg_pos)
  ) %>%
  ungroup() %>%
  # return annotation cols
  left_join(distinct(ppau_delta_paired,
                     across(c("le_id", "gene_id", "gene_name", "cryptic_status", "evaluated_status", "disease",
                              "plot_le_id", "simple_event_type"))
  ),
  by = c("le_id", "disease")
  ) %>%
  relocate(all_of(c("mean_paired_delta_ppau", "median_paired_delta_ppau")),
  .after = last_col()
  ) %>%
  # put cryptics first
  arrange(desc(cryptic_status))

if (!dir.exists("processed/liu_facs/")) {dir.create("processed/liu_facs/", recursive = T)}
# output median & mean pairwise delta_PPAU tables

ppau_delta_paired_median_all %>%
  filter(cryptic_status) %>%
  write_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_delta_ppau.all_samples.cryptics.tsv", col_names = T)

ppau_delta_paired_median_all %>%
  write_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_delta_ppau.all_samples.all_ales.tsv.gz", col_names = T)

ppau_delta_paired_median_disease %>%
  filter(cryptic_status) %>%
  write_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_delta_ppau.subtype_split.cryptics.tsv", col_names = T)

ppau_delta_paired_median_disease %>%
  write_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_delta_ppau.subtype_split.all_ales.tsv.gz", col_names = T)

# output per sample
ppau_delta_paired %>%
  filter(cryptic_status) %>%
  write_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_per_sample_delta_ppau.cryptics.tsv", col_names = T)

ppau_delta_paired %>%
  write_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_per_sample_delta_ppau.all_ales.tsv.gz", col_names = T)

# means calculated using means per population (negative vs positive) ignoring per-sample differences
ppau_meta_means %>%
  filter(cryptic_status) %>%
  write_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_popn_delta_ppau.cryptics.tsv", col_names = T)

ppau_meta_means %>%
  write_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_popn_delta_ppau.all_ales.tsv.gz", col_names = T)
