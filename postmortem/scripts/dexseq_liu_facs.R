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

process_dexseq_df <- function(df, le2name, med_ppau) {
  # add gene name to dexseq df
  df <- left_join(df, le2name, by = "le_id")
  # add median ppaus to dexseq df
  df <- left_join(df, select(med_ppau, -gene_name), by = "le_id")
  
  # sort by significance for easier exploration
  df <- arrange(df, gene.qvalue, padj)
  
  # remove redudnatn symbol labels for some events
  mutate(df,
         gene_name = collapse_names(gene_name),
         sig_005 = padj < 0.05,
         sig_01 = padj < 0.1)

}


res_df <- read_tsv("processed/liu_facs/2023-11-21_liu_facs_covariate_patient_only_TDPneg_TDPpos.results.tsv")
le2name <- read_tsv("processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.le2name.tsv")
med_ppau <- read_tsv("processed/2023-09-11_liu_facs_decoys_delta_ppau.all_samples.all_ales.tsv.gz")

# event type annots
res_df <- process_dexseq_df(res_df, le2name, med_ppau)
# pvalue histogram
pval_hist <- res_df %>%
  ggplot(aes(x = pvalue)) +
  geom_histogram(bins = 1000) +
  theme_bw(base_size = 20) +
  labs(subtitle = "Patient as covariate",
       x = "P-value (DEXSeq)",
       y = "Count")

# get number of cryptic events that are sig at each threshold
sig_cryp_counts <- res_df %>%
  filter(padj < 0.1) %>%
  pivot_longer(starts_with("sig_"),names_to = "sig_threshold", values_to = "sig", names_prefix = "sig_") %>%
  # for each threshold, select rep row for gene based on whether cryptic or not
  group_by(sig_threshold, gene_name) %>%
  slice_max(cryptic_status) %>%
  ungroup() %>%
  count(sig_threshold, cryptic_status, sig)

# get cryptic genes at each threshold
cryp_gn_005 <- res_df %>% filter(cryptic_status, sig_005) %>% pull(gene_name)
cryp_gn_01 <- res_df %>% filter(cryptic_status, sig_01) %>% pull(gene_name)

write_lines(cryp_gn_005, "processed/liu_facs/2023-11-21_liu_facs_covariate_patient_only_TDPneg_TDPpos.cryptics_sig_005.txt")
write_lines(cryp_gn_01, "processed/liu_facs/2023-11-21_liu_facs_covariate_patient_only_TDPneg_TDPpos.cryptics_sig_01.txt")
write_tsv(sig_cryp_counts, "processed/liu_facs/2023-11-21_liu_facs_covariate_patient_only_TDPneg_TDPpos.sig_gene_counts_cryptic_status.tsv")

ggsave("processed/liu_facs/2023-11-21_liu_facs_covariate_patient_only_TDPneg_TDPpos.pvalue_histogram.png",
       pval_hist,
       device = "png",
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina")