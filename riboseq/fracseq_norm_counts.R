library(tidyverse)
library(DESeq2)

# isoform-level counts for fracseq data
fracseq_counts <- read_tsv("processed/fracseq/2024-04-30_summarised_pas.counts.tsv")

# At least for part of the analysis, want to try to compare the proportion of total isoform expression in each fraction
# As comparing across samples/fractions, absolutely essential to normalise for differences in library depth and composition to prevent confounding
# Will use the standard DESeq/DEXSeq size factor method to do this

#
le2gene <- fracseq_counts %>%
  distinct(le_id, gene_id)

# convert to matrix with isoform ids as rownames, cols as samples
fracseq_counts_mtx <- fracseq_counts %>%
  select(-gene_id) %>%
  column_to_rownames("le_id")

# do separately for esc and npc cells
fracseq_counts_mtx_npc <- fracseq_counts_mtx %>%
  select(contains("NPC_short_read"))

fracseq_counts_mtx_esc <- fracseq_counts_mtx %>%
  select(contains("ESC_short_read"))


# Named vector of size factors for each sample in matrix
sfs_npc <- DESeq2::estimateSizeFactorsForMatrix(fracseq_counts_mtx_npc)
sfs_esc <-  DESeq2::estimateSizeFactorsForMatrix(fracseq_counts_mtx_esc)

# Divide each count column (sample) by its corresponding size factor
fracseq_counts_norm_npc <- sweep(x = fracseq_counts_mtx_npc,
                                MARGIN = 2, # operate on columns
                                STATS = sfs_npc,
                                FUN = '/')


fracseq_counts_norm_esc <- sweep(x = fracseq_counts_mtx_esc,
                                 MARGIN = 2, # operate on columns
                                 STATS = sfs_esc,
                                 FUN = '/')


# Tidy up dfs for output
fracseq_counts_norm_npc <- fracseq_counts_norm_npc %>%
  rownames_to_column("le_id") %>%
  left_join(le2gene, by = "le_id") %>%
  relocate(gene_id, .after = le_id)

fracseq_counts_norm_esc <- fracseq_counts_norm_esc %>%
  rownames_to_column("le_id") %>%
  left_join(le2gene, by = "le_id") %>%
  relocate(gene_id, .after = le_id)

# convert size factors into dfs
sfs_npc_df <- sfs_npc %>%
  enframe("sample_name", "size_factor")

sfs_esc_df <- sfs_esc %>%
  enframe("sample_name", "size_factor")


# write to disk
write_tsv(fracseq_counts_norm_npc, "processed/fracseq/2024-04-30_summarised_pas.counts.normalised.npc.tsv")
write_tsv(fracseq_counts_norm_esc, "processed/fracseq/2024-04-30_summarised_pas.counts.normalised.esc.tsv")
write_tsv(sfs_npc_df, "processed/fracseq/2024-04-30_size_factors.npc.tsv")
write_tsv(sfs_esc_df, "processed/fracseq/2024-04-30_size_factors.esc.tsv")
