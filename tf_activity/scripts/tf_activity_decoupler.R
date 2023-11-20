library(tidyverse)
library(fgsea)
library(dorothea)
library(decoupleR)
library(janitor)
source("scripts/helpers.R")

# output of dl_collectri.R
collectri_hs <- read_tsv("data/2023-11-15_collectri_homosapiens.tsv")
collectri_hs <- rename(collectri_hs, tf = source)

dorothea_hs_abc <- filter(dorothea_hs, confidence %in% c("A", "B", "C"))

# 1. Assess different scoring methods and target gene sets in an ELK1 KD
# Previously used Dorothea DB, now a new collectRI database with updated targets (+ for SIX3, TLX1)
# For ELK1, want to know whether these curated target lists can predict ELK1 TF activity as perturbed in this dataset

# a) compare dorothea & collectRI - do they identify TF as putatively repressed?
# TFome wide, is ELK1 considered rhe most deregulated

elk1_kd <- read_csv("data/knockTF/KnockTF-Search-ELK1-GSE46333.csv")
elk1_kd <- elk1_kd %>% 
  janitor::clean_names(case = "snake")

elk1_kd %>%
  filter(target_gene == "ELK1")
# - 1.75 FC

# add -log10 signed pvalue, -log10 pvalue * fc as columns
elk1_kd <- elk1_kd %>%
  add_signed_pval(fc_col = "log2fc", pval_col = "p_value") %>%
  add_log2fold_pval(fc_col = "log2fc", pval_col = "p_value")


# get matrix of diff exprn results for running decoupler methdos
elk1_kd_mtx <- elk1_kd %>%
  select(gene_name = target_gene, log2fc, signed_pvalue,pvalScaledLog2FoldChange) %>%
  drop_na(gene_name) %>%
  column_to_rownames(var = "gene_name") %>%
  as.matrix()

# try dorothea & collectri targets on all 3 stats
elk1_kd_decoupler_ulm_collectri <- c("log2fc", "signed_pvalue", "pvalScaledLog2FoldChange") %>%
  set_names() %>%
  map(~ run_ulm(elk1_kd_mtx[, .x, drop=FALSE],
               network = collectri_hs,
               .source = "tf",
               .target = "target",
               .mor = "mor")
        ) %>%
  bind_rows(.id = "score_metric")

elk1_kd_decoupler_ulm_dorothea <- c("log2fc", "signed_pvalue", "pvalScaledLog2FoldChange") %>%
  set_names() %>%
  map(~ run_ulm(elk1_kd_mtx[, .x, drop=FALSE],
                network = dorothea_hs_abc,
                .source = "tf",
                .target = "target",
                .mor = "mor")
  ) %>%
  bind_rows(.id = "score_metric")


# ELk1 only just reaches nominal significance threshold with dorothea targets
# -log10(pvalue) * log2fold is consistently most extreme however.


# b) compare scoring methods
# what about consensus score?
# by default, ulm, mlm & norm_wsum
elk1_kd_decoupler_consensus_collectri <- c("log2fc", "signed_pvalue", "pvalScaledLog2FoldChange") %>%
  set_names() %>%
  map(~ decouple(elk1_kd_mtx[, .x, drop=FALSE],
                 network = collectri_hs,
                 .source = "tf",
                 .target = "target")
  ) %>%
  bind_rows(.id = "score_metric")

elk1_kd_decoupler_consensus_dorothea <- c("log2fc", "signed_pvalue", "pvalScaledLog2FoldChange") %>%
  set_names() %>%
  map(~ decouple(elk1_kd_mtx[, .x, drop=FALSE],
                 network = dorothea_hs_abc,
                 .source = "tf",
                 .target = "target")
  ) %>%
  bind_rows(.id = "score_metric")


# Again, ELK1 rarely reaches nominal significance.
# Maybe this isn't a useful activity, because we know ELK1 LoF can be compensated

# 2. Liu FACS
# Since no obvious leader, going to try all gene sets + consensus score (as gives all 3 in one go)

# deseq size-factor normalised counts
liu_facs_counts <- read_tsv("data/liu_facs/2023-09-07_liu_facs_salmon_summarised_counts_deseq_normalised.tsv")

# mtx ready for decoupler
liu_facs_counts_mtx <- liu_facs_counts %>%
  select(-all_of(c("gene_id", "FTD_S3_unsorted","FTD_S4_unsorted"))) %>%
  drop_na(gene_name) %>%
  # for now remove duplicated genes (only 11 - seem to have diff gene IDs)
  group_by(gene_name) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  # log transform the observed counts, adding a pseudocount of 1
  mutate(across(-gene_name, ~ log2(.x + 1))) %>%
  column_to_rownames(var = "gene_name") %>%
  as.matrix()

# run consensus for each DB on dorothea data
liu_facs_decoupler_consensus <- list(collectri = collectri_hs, dorothea_abc = dorothea_hs_abc) %>%
  map(~ decouple(liu_facs_counts_mtx,
                 network = .x,
                 .source = "tf",
                 .target = "target",args = list(wsum = list(times = 1000)))) %>%
  bind_rows(.id = "database")


#1. Calculate ranks by score/p-value
### Split by DB + metric

# pvalue histograms for each db + score
liu_facs_decoupler_consensus %>%
  ggplot(aes(x = p_value)) +
  facet_wrap("database ~ statistic", scales = "free_y") +
  geom_histogram(bins = 200) + 
  theme_bw(base_size = 20)

# pvalues only appear globally well-behaved in case of collectri mlm, ulm + dorothea mlm
# wsum, corr_wsum & norm_wsum are quite sparse, probably due to number of permutations used

#dorothea pvalues are incredibly 0 inflated (i.e. probably invalid)
# NB: not sure whether to expect a inflation near zero, as that suggests dataset is enriched for many TFs being changed

#2. plot the activity scores by patient + subtype - expect to be differential with TDP status
liu_facs_decoupler_consensus %>%
  filter(source == "ELK1") %>%
  mutate(tdp_status = if_else(str_ends(condition, "_positive"),
                              "TDPpos", "TDPneg"),
         sample_id = str_remove(condition, "_TDP_43_(positive|negative)$")
         ) %>%
  ggplot(aes(x = sample_id, y = score, fill = tdp_status)) +
  facet_wrap("database ~ statistic", scales = "free_y") +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
         
