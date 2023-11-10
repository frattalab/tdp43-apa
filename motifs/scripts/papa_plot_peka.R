library(tidyverse)
source("scripts/fncs_plot_peka.R")


dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-03_papa_cryptics_kmer6_window_250_distal_window_500.cleaned_6mer_distribution_genome.tsv")

dbrn_tbl_long_gu <- peka_wide_to_long(dbrn_tbl, c("GUGUGU"), first_posn_idx = 10)
dbrn_tbl_long_ug <- peka_wide_to_long(dbrn_tbl, c("UGUGUG"), first_posn_idx = 10)

dbrn_tbl_long_both <- peka_wide_to_long(dbrn_tbl, c("UGUGUG", "GUGUGU"), first_posn_idx = 10, sum_occur = T, sum_group_cols = c( "comparison_name", "rel_posn"))

plot_kmer_dbrn(mutate(dbrn_tbl_long_gu, direction = "cryptic"),
               rolling_mean = T,
               facet_w = "~ comparison_name",
               n_col = 2,
               n_row = 3,
               breaks = seq(-500,500,100),
               rolling_k = 10
               )

plot_kmer_dbrn(mutate(dbrn_tbl_long_ug, direction = "cryptic"),
               rolling_mean = T,
               facet_w = "~ comparison_name",
               n_col = 2,
               n_row = 3,
               breaks = seq(-500,500,100),
               rolling_k = 10
)

plot_kmer_dbrn(mutate(dbrn_tbl_long_both, direction = "cryptic"),
               rolling_mean = T,
               facet_w = "~ comparison_name",
               n_col = 2,
               n_row = 3,
               breaks = seq(-500,500,100),
               rolling_k = 10
)
