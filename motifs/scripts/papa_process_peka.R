library(tidyverse)

# find paths to all kmer distrivution tables for different comparisons
dbrn_paths <- list.files(path = "data/peka_papa/2023-11-27_papa_cryptics_fixed_kmer6_window_250_distal_window_500_relpos_0",
                         pattern = "_6mer_distribution_genome.tsv$",
                         recursive = T,
                         full.names = T) %>%
  set_names(str_remove(basename(.), ".all_6mer_distribution_genome.tsv$"))

# read in kmer distribution tables for all experiments
l_dbrn_tbls <- map(dbrn_paths,
                   ~ read_tsv(.x, show_col_types = F) %>% rename(kmer = `...1`,
                                                                 PEKA_score = `PEKA-score`,
                                                                 pvalue = `p-value`
                                                                 )
)

#combine into single table, 
dbrn_tbl <- bind_rows(l_dbrn_tbls, .id = "comparison_name")

# remove position cols for easier exploration
simp_dbrn_tbl <- dbrn_tbl %>%
  select(kmer,
         comparison_name,
         kmer,
         mtxn,
         artxn,
         aroxn,
         etxn,
         PEKA_score,
         pvalue)

if (!dir.exists("processed/peka/papa")) {dir.create("processed/peka/papa", recursive = T)}

write_tsv(dbrn_tbl,
          "processed/peka/papa/2023-11_27_papa_cryptics_kmer6_window_250_distal_window_500_relpos_0.cleaned_6mer_distribution_genome.tsv",
          col_names = T)

write_tsv(simp_dbrn_tbl,
          "processed/peka/papa/2023-11-27_papa_cryptics_kmer6_window_250_distal_window_500_relpos_0.cleaned_6mer_distribution_genome_simple.tsv",
          col_names = T)