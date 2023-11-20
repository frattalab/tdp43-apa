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

