library(tidyverse)
source("scripts/fncs_plot_peka.R")


dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-02_papa_cryptics_kmer6_window_200_distal_window_250.cleaned_6mer_distribution_genome.tsv")

dbrn_tbl_long <- peka_wide_to_long(dbrn_tbl, c("GUGUGU"), first_posn_idx = 10)


plot_kmer_dbrn(mutate(dbrn_tbl_long, direction = "cryptic"),
               rolling_mean = T,
               facet_w = "~ comparison_name",
               n_col = 2,
               n_row = 3,
               breaks = seq(-250,250,25)
               )
