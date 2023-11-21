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
liu_facs_deseq <- read_tsv("../postmortem/processed/gene_exprn/deseq2/2023-20-11_deseq2_liu_facs_results.tsv")

# mtx ready for decoupler
liu_facs_deseq_mtx <- liu_facs_deseq %>%
  group_by(gene_name) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(gene_name, stat, log2FoldChange, log2FoldChangeShrink, pvalue, padj) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
  
liu_facs_deseq_decoupler_consensus <- list(collectri = collectri_hs, dorothea_abc = dorothea_hs_abc) %>%
  map(~ decouple(liu_facs_deseq_mtx[, "stat", drop=FALSE],
                 network = .x,
                 .source = "tf",
                 .target = "target",args = list(wsum = list(times = 1000)))) %>%
  bind_rows(.id = "database")
print("deseq stat done")

#1. Calculate ranks by score/p-value
### Split by DB + metric

# pvalue histograms for each db + score
liu_facs_deseq_decoupler_consensus %>%
  ggplot(aes(x = p_value)) +
  facet_wrap("database ~ statistic", scales = "free_y") +
  geom_histogram(bins = 100) + 
  theme_bw(base_size = 20)


#2. plot the activity scores by patient + subtype - expect to be differential with TDP status
liu_facs_deseq_decoupler_consensus %>%
  ggplot(aes(x = score)) +
  facet_wrap("database ~ statistic", scales = "free_y") +
  geom_density() + 
  theme_bw(base_size = 20)

# add rank within db + statistic
liu_facs_deseq_decoupler_consensus <- liu_facs_deseq_decoupler_consensus %>%
  group_by(database, statistic) %>%
  arrange(desc(abs(score)), .by_group = T) %>%
  mutate(rank = row_number(),
         frac_rank = rank / n()) %>%
  ungroup()


elk1_hela_scores_liu_plot <- liu_facs_deseq_decoupler_consensus %>%
  filter(source == "ELK1") %>%
  ggplot(aes(x = statistic, y = frac_rank * 100, label = round(p_value, 3))) +
  facet_wrap("~ database") +
  geom_point() + 
  geom_text(nudge_y = 5) +
  scale_y_continuous(limits = c(0,100),
                     breaks = seq(0,100,10)) +
  theme_bw(base_size = 20) +
  labs(title = "ELK1 inferred TF activity in Liu FACS",
       subtitle = "text label = raw p-value",
       x = "",
       y = "Score percentile"
       ) +
  theme(axis.text.x = element_text(angle = 90))

ggsave("processed/2023-11-21_liu_facs_elk1_decoupler_tf_activity_percentile.png",
       height = 10,
       width = 10,
       units = "in",
       dpi = "retina")
         
