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

p_le2name <- "processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.le2name.tsv"

p_sample_tbl <- "data/liu_facs_papa_sample_sheet.csv"

ppau <- read_tsv(p_ppau)
tx2le <- read_tsv(p_tx2le)
le2name <- read_tsv(p_le2name) %>%
  # wasn't smart with gene_name collapsing, often get duplicates...
  mutate(gene_name = collapse_names(gene_name))
sample_tbl <- read_csv(p_sample_tbl)

# get a vector of cryptic last exon IDs
cryp_le_ids <- tx2le %>%
  distinct(le_id) %>% 
  pull()

# join in gene name for each le_id
ppau <- left_join(ppau, distinct(le2name, .keep_all = T), by = "le_id") %>% 
  relocate(gene_name, .after = gene_id) %>%
  # don't need these columns
  select(-starts_with("mean"), -starts_with("delta"))

# join in metadata from sample table (converting to long format with row per sample and isoform)
ppau <- ppau %>%
  pivot_longer(cols = -all_of(c("le_id", "gene_id", "gene_name")),
               names_to = "sample_name",
               values_to = "ppau") %>%
  left_join(select(sample_tbl, -path, -fastq1, -fastq2), by = "sample_name")


# now want to calculate a sample-wise difference in PPAU between TDP-43 negative & positive samples (condition column)
ppau_delta_paired <- ppau %>%
  select(-sample_name) %>%
  # for each sample & le_id get ppau values in positive and negative as columns
  pivot_wider(id_cols = c("le_id", "gene_id", "gene_name", "patient_id", "disease", "gender"),
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
median_ppau_order <- ""

# rank deltas from largest to smallest, giving same rank if identical delta
mutate(rank_mean_up = min_rank(desc(mean_paired_delta_ppau)),
       rank_median_up = min_rank(desc(median_paired_delta_ppau))) %>%
  arrange(rank_median_up, rank_mean_up)


# rejoin gene information 
ppau_delta_paired_median_all <- ppau_delta_paired_median_all %>%
  left_join(distinct(ppau_delta_paired, le_id, gene_id, gene_name), by = "le_id") %>%
  relocate(all_of(c("mean_paired_delta_ppau", "median_paired_delta_ppau",
                    "rank_mean_up", "rank_median_up")),
           .after = last_col()
           )





# use ggside to add gender (& possibly gene expression) as side tiles to heatmap
# https://cran.r-project.org/web/packages/ggside/vignettes/ggside_basic_usage.html

