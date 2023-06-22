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

p_dexseq <- "data/PAPA/2023-05-24_i3_cortical_zanovello.all_datasets.dexseq_apa.results.processed.cleaned.tsv"
p_orig_tx2le <- "data/PAPA/novel_ref_combined.tx2le.tsv"

p_le2name <- "processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.le2name.tsv"

p_sample_tbl <- "data/liu_facs_papa_sample_sheet.csv"

ppau <- read_tsv(p_ppau)
tx2le <- read_tsv(p_tx2le)
le2name <- read_tsv(p_le2name) %>%
  # wasn't smart with gene_name collapsing, often get duplicates...
  mutate(gene_name = collapse_names(gene_name))
sample_tbl <- read_csv(p_sample_tbl)

## 
dexseq <- read_tsv(p_dexseq)
orig_tx2le <- read_tsv(p_orig_tx2le)
dexseq_cryp <- dexseq %>% filter(padj < 0.05 & mean_PPAU_base < 0.1 & delta_PPAU_treatment_control > 0.1)

# get tx_ids of cryptic le_ids
orig_tx_cryp <- orig_tx2le %>% 
  filter(le_id %in% dexseq_cryp$le_id) %>%
  pull(transcript_id)

# get a vector of cryptic last exon IDs (for current annotation)
cryp_le_ids <- tx2le %>%
  filter(transcript_id %in% orig_tx_cryp) %>%
  distinct(le_id) %>%
  pull()

# extract event type annotation
le2event <- dexseq %>%
  distinct(le_id, simple_event_type) %>%
  group_by(le_id) %>%
  summarise(simple_event_type = paste(unique(sort(simple_event_type)), collapse = ",")) %>%
  ungroup()

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
# annotate whether cryptic/nto cryptic
ppau <- left_join(ppau, distinct(le2name, .keep_all = T), by = "le_id") %>% 
  relocate(gene_name, .after = gene_id) %>%
  # don't need these columns
  select(-starts_with("mean"), -starts_with("delta")) %>%
  mutate(cryptic_status = if_else(le_id %in% cryp_le_ids, T, F)) %>%
  left_join(le2event, by = "le_id")

# dfine a more readable isoform id - gene_name suffixed with isoform number
ppau <- ppau %>%
  group_by(gene_name) %>%
  arrange(le_id, .by_group = T) %>%
  mutate(le_n = row_number(),
         plot_le_id = paste(gene_name, le_n, sep = "_")) %>%
  ungroup() %>%
  select(-le_n)
  

# join in metadata from sample table (converting to long format with row per sample and isoform)
ppau <- ppau %>%
  pivot_longer(cols = -all_of(c("le_id", "gene_id", "gene_name", "cryptic_status", "plot_le_id", "simple_event_type")),
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
  group_by(le_id) %>%
  summarise(mean_paired_delta_ppau = mean(paired_delta_ppau_neg_pos),
            median_paired_delta_ppau = median(paired_delta_ppau_neg_pos),
            ) %>%
  ungroup() %>%
  left_join(distinct(ppau_delta_paired, le_id, gene_id, gene_name, cryptic_status, plot_le_id), by = "le_id") %>%
  relocate(all_of(c("mean_paired_delta_ppau", "median_paired_delta_ppau"
                    )),
           .after = last_col()
  )


# calculate mean & median delta PAUs across samples for each isoform
ppau_delta_paired_median_disease <- ppau_delta_paired %>%
  group_by(le_id, disease) %>%
  summarise(mean_paired_delta_ppau = mean(paired_delta_ppau_neg_pos),
            median_paired_delta_ppau = median(paired_delta_ppau_neg_pos),
            rank_mean_up = min_rank(desc(mean_paired_delta_ppau)),
            rank_median_up = min_rank(desc(median_paired_delta_ppau))
  ) %>%
  ungroup() %>%
  arrange(rank_median_up, rank_mean_up)


# get an order of le_ids according to ranks of median delta between negative and positive
# (this puts events with most consistent enrichment at the top of the heatmap)
ppau_order_cryp <- ppau_delta_paired_median_all %>%
  filter(cryptic_status) %>%
  # rank deltas from largest to smallest, giving same rank if identical delta
mutate(rank_mean_up = min_rank(desc(mean_paired_delta_ppau)),
       rank_median_up = min_rank(desc(median_paired_delta_ppau))) %>%
  arrange(rank_median_up, rank_mean_up)

ppau_order_cryp_median <- ppau_order_cryp %>%
  arrange(rank_median_up) %>%
  pull(le_id)

ppau_order_cryp_median_gn <- ppau_order_cryp %>%
  arrange(rank_median_up) %>%
  pull(plot_le_id)

ppau_delta_paired_median_all %>%
  filter(cryptic_status) %>%
  summary(paired_delta_ppau_neg_pos)

# rejoin gene information 
# ppau_delta_paired_median_all <- ppau_delta_paired_median_all %>%
#   left_join(distinct(ppau_delta_paired, le_id, gene_id, gene_name), by = "le_id") %>%
#   relocate(all_of(c("mean_paired_delta_ppau", "median_paired_delta_ppau",
#                     "rank_mean_up", "rank_median_up")),
#            .after = last_col()
#            )


ppau_delta_paired_cryp <- ppau_delta_paired %>%
  filter(cryptic_status)

ppau_delta_paired_cryp %>%
  mutate(le_id = factor(le_id, levels = rev(ppau_order_cryp_median)),
         plot_le_id = factor(plot_le_id, levels = rev(ppau_order_cryp_median_gn))
         ) %>%
  ggplot(aes(x = patient_id,
             y = le_id,
             fill = paired_delta_ppau_neg_pos)) +
  geom_tile() +
  scale_fill_gradientn(name = "Sample-wise dPPAU (TDPneg - TDPpos)",
                       colours = c("#998ec3", "#f7f7f7", "#f1a340"),
                       limits = c(-1, 1)) +
  theme_bw() +
  labs(title = "Cell line cryptic last exons - relative usage change in FACS nuclei",
       x = "Sample ID",
       y = "Gene name") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y = element_text(size=rel(0.82)))

ppau_delta_paired_cryp %>%
  mutate(le_id = factor(le_id, levels = rev(ppau_order_cryp_median)),
         plot_le_id = factor(plot_le_id, levels = rev(ppau_order_cryp_median_gn))
  ) %>%
  ggplot(aes(x = patient_id,
             y = plot_le_id,
             fill = paired_delta_ppau_neg_pos)) +
  geom_tile() +
  scale_fill_gradientn(name = "Sample-wise dPPAU (TDPneg - TDPpos)",
                       colours = c("#998ec3", "#f7f7f7", "#f1a340"),
                       limits = c(-1, 1)) +
  theme_bw() +
  labs(title = "Cell line cryptic last exons - relative usage change in FACS nuclei",
       x = "Sample ID",
       y = "Gene name") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y = element_text(size=rel(0.82)))



ppau_order_cryp %>%
  filter(median_paired_delta_ppau > 0.05) %>%
  pull(plot_le_id)


ppau_order_cryp %>%
  filter(mean_paired_delta_ppau > 0.05)


# use ggside to add gender (& possibly gene expression) as side tiles to heatmap, event type
# https://cran.r-project.org/web/packages/ggside/vignettes/ggside_basic_usage.html

