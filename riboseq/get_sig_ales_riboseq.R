library(tidyverse)

df <- read_tsv("processed/2023-10-18_i3_cryptic_genes_riboseq_deseq_summary.tsv")

# get all diff translated ALEs
df_sig <- df %>% filter(diff_translated)

df_sig_ale <- df_sig %>%
  filter(simple_event_type == "spliced") %>%
  arrange(desc(abs(log2FoldChangeShrink)))


# remove TRAPPC12 as know overlaps with a cassette exon
df_sig_ale <- df_sig_ale %>%
  filter(gene_name != "TRAPPC12")


if (!dir.exists("processed")) {dir.create("processed", recursive = T)}

write_tsv(df_sig_ale, "processed/2024-04-04_i3_cryptic_ale_genes_sig_riboseq_summary.tsv")