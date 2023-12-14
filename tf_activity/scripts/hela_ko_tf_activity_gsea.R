library(tidyverse)
library(fgsea)
library(dorothea)
library(decoupleR)
source("scripts/helpers.R")
set.seed(123)

# target genes with an ELK1 chIP seq peak within 1kb of TSS
chipatlas_elk1 <- read_tsv("data/chip_atlas/2023-11-15_chipatlas_tss_1kb_ELK1.tsv")
ferguson_deseq <- read_csv("data/tdp43_kd_collection/deseq2_outputs/ferguson.hela.DESEQ2_results.csv")
ferguson_counts <- read_csv("data/tdp43_kd_collection/deseq2_outputs/ferguson.hela.DESEQ2_normedcounts.csv")
ferguson_splicing <- read_csv("data/tdp43_kd_collection/splicing/ferguson_collected_splicing.csv")

# any genes with 'significant' splicing change
ferguson_diff_spliced <- ferguson_splicing %>%
  filter(probability_changing > 0.95) %>%
  distinct(gene_name) %>%
  pull()

#
ferguson_cryptic <- ferguson_splicing %>%
  filter(baseline_PSI < 0.05 & mean_dpsi_per_lsv_junction > 0.20 & probability_changing > 0.95) %>%
  distinct(gene_name) %>%
  pull()


# output of dl_collectri.R
collectri_hs <- read_tsv("data/2023-11-15_collectri_homosapiens.tsv")

# shoudl only be two datasets 
cols_hela <- colnames(chipatlas_elk1)[str_detect(colnames(chipatlas_elk1), "HeLa")]
length(cols_hela) == 2

# subset for hela targets
hela_targets_elk1 <- chipatlas_elk1 %>%
  select(Target_genes, contains("HeLa")) %>%
  filter(!!sym(cols_hela[1]) != 0 | !!sym(cols_hela[2]) != 0)

# subset for target genes called in both HeLa datasets
hela_targets_elk1_both <- hela_targets_elk1 %>%
  filter(!!sym(cols_hela[1]) != 0 & !!sym(cols_hela[2]) != 0) 




# ELK4
chipatlas_elk4 <- read_tsv("data/chip_atlas/2023-11-29_chipatlas_tss_1kb_ELK4.tsv")
cols_hela_elk4 <- colnames(chipatlas_elk4)[str_detect(colnames(chipatlas_elk4), "HeLa")]

# targets where number of datasets
elk4_num_datasets <- chipatlas_elk4 %>%
  select(target = Target_genes, ends_with("HeLa")) %>%
  pivot_longer(cols = ends_with("HeLa"), names_to = "dataset") %>%
  group_by(target) %>%
  summarise(n_datasets = sum(value != 0))


hela_targets_elk4_all <- filter(chipatlas_elk4, Target_genes %in% pull(filter(elk4_num_datasets, n_datasets == 3), target))


# get overlapping and unique targets
hela_targets_elk1_elk4_shared <- intersect(hela_targets_elk1_both$Target_genes, hela_targets_elk4_all$Target_genes)
hela_targets_elk1_unique <- setdiff(hela_targets_elk1_both$Target_genes, hela_targets_elk1_elk4_shared)
hela_targets_elk4_unique <- setdiff(hela_targets_elk4_all$Target_genes, hela_targets_elk1_elk4_shared)

length(hela_targets_elk1_elk4_shared) / length(hela_targets_elk1_both$Target_genes)
length(hela_targets_elk1_elk4_shared) / length(hela_targets_elk4_all$Target_genes)


# define a list of target gene IDs for use with fgsea
chipseq_hela_target_lists <- list(chipseq_hela_both_elk1 = pull(hela_targets_elk1_both, Target_genes),
                                  chipseq_hela_elk1_elk4 = hela_targets_elk1_elk4_shared,
                                  chipseq_hela_unique_elk1 = hela_targets_elk1_unique,
                                  chipseq_hela_unique_elk4 = hela_targets_elk4_unique,
                                  chipseq_hela_all_elk4 = pull(hela_targets_elk4_all, Target_genes) )

# remove genes with NA padj value, keep one row per gene name
ferguson_deseq <- clean_deseq_df(ferguson_deseq)

# calculate signed -log10 pvalue (multiplied by sign of log2Fold change)
# also handles exact matches of pvalues (genes with larger FCs are prioritised)
ferguson_deseq <- add_signed_pval(ferguson_deseq,)

# also calculate absolute stat and -log10 pvalue
ferguson_deseq <- ferguson_deseq %>%
  mutate(abs_stat = abs(stat),
         abs_signed_pvalue = abs(signed_pvalue))

# create a list of genes ranked by various metrics
# stat, signed_pvalue, pval 
ferguson_deseq_ranks <- c("stat", "signed_pvalue") %>%
  set_names() %>%
  map(~ get_ranked_gene_list(ferguson_deseq, score_col = .x))

# repeat for absolute cols
ferguson_deseq_ranks_abs <- c("abs_stat", "abs_signed_pvalue") %>%
  set_names() %>%
  map(~ get_ranked_gene_list(ferguson_deseq, score_col = .x))

# list of genes with diff spliced removed (+ TARDBP) (i.e. effect on expression likely accounted for by splicing)
ferguson_deseq_ranks_nospl <- map(ferguson_deseq_ranks,
                                  ~ .x[!names(.x) %in% c(ferguson_diff_spliced, "TARDBP")])

ferguson_deseq_ranks_abs_nospl <- map(ferguson_deseq_ranks_abs,
                                  ~ .x[!names(.x) %in% c(ferguson_diff_spliced, "TARDBP")])


# now can run fgsea using all 3 different stats (to see how stable results are)

gsea_ferguson_chipseq <- map(ferguson_deseq_ranks,
                             ~ fgsea(pathways = chipseq_hela_target_lists,
                                     stats = .x,
                                     eps = 0,
                                     nproc = 2),
                             .progress = T
                             ) %>%
  bind_rows(.id = "score_type")

gsea_ferguson_nospl_chipseq <- map(ferguson_deseq_ranks_nospl,
                             ~ fgsea(pathways = chipseq_hela_target_lists,
                                     stats = .x,
                                     eps = 0,
                                     nproc = 2),
                             .progress = T
) %>%
  bind_rows(.id = "score_type")

# repeat for abs
gsea_ferguson_chipseq_abs <- map(ferguson_deseq_ranks_abs,
                             ~ fgsea(pathways = chipseq_hela_target_lists,
                                     stats = .x,
                                     eps = 0,
                                     scoreType = "pos",
                                     nproc = 2),
                             .progress = T
) %>%
  bind_rows(.id = "score_type")

gsea_ferguson_nospl_chipseq_abs <- map(ferguson_deseq_ranks_abs_nospl,
                                   ~ fgsea(pathways = chipseq_hela_target_lists,
                                           stats = .x,
                                           eps = 0,
                                           scoreType = "pos",
                                           nproc = 2),
                                   .progress = T
) %>%
  bind_rows(.id = "score_type")


# log2fold * pvalue seems the most unstable,
# stat and signed-pvalue are pretty consistent in relative terms. CHIP-seq both is by far the weakest by significance, but still consistent directionality

# get list of ELK1 dorothea targets split by direction for fgsea input
dorothea_elk1_target_list <- dorothea_hs %>%
  filter(tf == "ELK1") %>%
  dorothea_to_gsea()

#  do not split by mode of regulation
dorothea_elk1_target_list_abs <- dorothea_hs %>%
  filter(tf == "ELK1") %>%
  dorothea_to_gsea(split_by_mor = F)


# Subset Dorothea to ELK1 HeLA CHIP-seq tagrets (so have expected mode of regulation for these genes)
dorothea_chip_elk1_target_list <- hela_targets_elk1 %>%
  select(target = Target_genes) %>%
  mutate(tf = "ELK1") %>%
  inner_join(dorothea::dorothea_hs, by = c("tf", "target")) %>%
  dorothea_to_gsea()


dorothea_chip_elk4_target_list <- hela_targets_elk4_all %>%
  select(target = Target_genes) %>%
  mutate(tf = "ELK4") %>%
  inner_join(dorothea::dorothea_hs, by = c("tf", "target")) %>%
  dorothea_to_gsea()


names(dorothea_elk1_target_list) <- paste("dorothea", names(dorothea_elk1_target_list), sep = "_")
names(dorothea_elk1_target_list_abs) <- paste("dorothea", names(dorothea_elk1_target_list_abs), sep = "_")
names(dorothea_chip_elk1_target_list) <- paste("dorothea_chipseq", names(dorothea_chip_elk1_target_list), sep = "_")
names(dorothea_chip_elk4_target_list) <- paste("dorothea_chipseq", names(dorothea_chip_elk4_target_list), sep = "_")

# repeat for collectri
collectri_elk1_target_list <- collectri_hs %>%
  filter(source == "ELK1") %>%
  rename(tf = source) %>%
  dorothea_to_gsea()

collectri_elk1_target_list_abs <- collectri_hs %>%
  filter(source == "ELK1") %>%
  rename(tf = source) %>%
  dorothea_to_gsea(split_by_mor = F)


collectri_chip_elk1_target_list <- hela_targets_elk1 %>%
  select(target = Target_genes) %>%
  mutate(source = "ELK1") %>%
  inner_join(collectri_hs, by = c("source", "target")) %>%
  rename(tf = source) %>%
  dorothea_to_gsea()
  
  
collectri_chip_elk4_target_list <- hela_targets_elk4_all %>%
  select(target = Target_genes) %>%
  mutate(source = "ELK4") %>%
  inner_join(collectri_hs, by = c("source", "target")) %>%
  rename(tf = source) %>%
    dorothea_to_gsea()

names(collectri_elk1_target_list) <- paste("collectri", names(collectri_elk1_target_list), sep = "_")
names(collectri_elk1_target_list_abs) <- paste("collectri", names(collectri_elk1_target_list_abs), sep = "_")
names(collectri_chip_elk1_target_list) <- paste("collectri_chipseq", names(collectri_chip_elk1_target_list), sep = "_")
names(collectri_chip_elk4_target_list) <- paste("collectri_chipseq", names(collectri_chip_elk4_target_list), sep = "_")

# final combined list of all targets
elk1_target_list <- c(dorothea_elk1_target_list,
                      collectri_elk1_target_list,
                      chipseq_hela_target_lists,
                      collectri_chip_elk1_target_list,
                      collectri_chip_elk4_target_list,
                      dorothea_chip_elk1_target_list,
                      dorothea_chip_elk4_target_list
                      )

elk1_target_list_abs <- c(dorothea_elk1_target_list_abs,
                          collectri_elk1_target_list_abs,
                          chipseq_hela_target_lists
                          )


# run for final gene-sets
gsea_ferguson_all <- map(ferguson_deseq_ranks,
                         ~ fgsea(pathways = elk1_target_list,
                                 stats = .x,
                                 eps = 0,
                                 minSize = 3),
                         .progress = T
) %>%
  bind_rows(.id = "score_type")

gsea_ferguson_nospl_all <- map(ferguson_deseq_ranks_nospl,
                         ~ fgsea(pathways = elk1_target_list,
                                 stats = .x,
                                 eps = 0,
                                 minSize = 3),
                         .progress = T
) %>%
  bind_rows(.id = "score_type")

gsea_ferguson_all_abs <- map(ferguson_deseq_ranks,
                         ~ fgsea(pathways = elk1_target_list_abs,
                                 stats = .x,
                                 eps = 0,
                                 minSize = 3,
                                 scoreType = "pos"),
                         .progress = T
) %>%
  bind_rows(.id = "score_type")

gsea_ferguson_nospl_all_abs <- map(ferguson_deseq_ranks_nospl,
                               ~ fgsea(pathways = elk1_target_list_abs,
                                       stats = .x,
                                       eps = 0,
                                       minSize = 3,
                                       scoreType = "pos"),
                               .progress = T
) %>%
  bind_rows(.id = "score_type")

# As before:
# ranking by signed pvalue gives modest significant enrichment for downreg genes with ChIP-seq, but stat a strong enrichment
# Dorothea & collectri targets, which include sign/expected mode of regulation, show no significant association with up/downregulation
# Dorothea expected pos are biased towards pos NES, suggesting that ELK1's activity has increased. Does that suggest ChIP-seq genes are broadly downregulated/expected to be repressed?

# TODO: remove diff spliced genes (expression change likely due to TDP-43's splicing function) + TDP-43 (expression down due to deletion) & see effect on observed enrichment
# TODO: try univariate linear model with collect RI/dorothea


# Notes:
#

# padjs and extract out ELK1

# repeat, but remove differentially spliced
if (!dir.exists("processed/ferguson_hela")) {dir.create("processed/ferguson_hela", recursive = T)}

save.image("processed/ferguson_hela/2023-11-29_hela_ko_tf_activity_gsea.Rdata")
