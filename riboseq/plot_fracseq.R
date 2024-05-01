library(tidyverse)

normed_counts_npc <- read_tsv("processed/fracseq/2024-04-30_summarised_pas.counts.normalised.npc.tsv")
le2name <- read_tsv("../postmortem/processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.le2name.tsv")
sample_tbl <- read_csv("data/fracseq/ritter_short_read_fracseq_sample_sheet_minus_esc_hp_rep3.csv")

# quant done with rnaseq-single-steps, which uses basenames of files as sample names
# need to prepend SRR accession to sample name to ensure matches
sample_tbl <- mutate(sample_tbl,
                     sample_name_counts = paste(unit, sample_name, sep = "_")
                     )

# miniimal metadata for each sample
sample_tbl_meta <- select(sample_tbl, sample_name_counts, group, cell_type)

normed_counts_npc <- normed_counts_npc %>%
  left_join(le2name, by = "le_id") %>%
  relocate(gene_name, .after = gene_id)

# longer format (single row per sample)
normed_counts_npc_long <- pivot_longer(normed_counts_npc,
             cols = contains("NPC_short_read"), names_to = "sample_name_counts", values_to = "count") %>%
  left_join(sample_tbl_meta, by = "sample_name_counts") 

# Extract replicate label (so can assess proportion within replciates)
normed_counts_npc_long <- normed_counts_npc_long %>%
  mutate(replicate = str_extract_all(sample_name_counts, "rep[0-9]$", simplify = T)) %>%
  relocate(count, .after = everything())

# add a cleaned group label (for plotting)
normed_counts_npc_long <- normed_counts_npc_long %>%
mutate(plot_group = str_remove_all(group, "^NPC_"),
       plot_group = str_remove_all(plot_group, "^(light|heavy)_polyribosome_"),
       plot_group = factor(plot_group, levels = c("cytosol",
                                                  "monosome",
                                                  "2-4_ribosomes",
                                                  "5_ribosomes")
       )
)


# subset to ELK1, SIX3 and TLX1
normed_counts_npc_long_3exts <- normed_counts_npc_long %>%
  filter(gene_name %in% c("ELK1", "SIX3", "TLX1")) %>%
  mutate(event_label = if_else(str_ends(le_id, "_2"), "Cryptic", "Proximal"))


# calculate propn of total isoform counts for each -some fraction (replicate-wise)
normed_counts_npc_long_3exts <- normed_counts_npc_long_3exts %>%
  group_by(replicate, event_label) %>%
  mutate(rep_fraction = count / sum(count)) %>%
  ungroup()


normed_counts_npc_long_3exts %>%
  group_by(gene_name, plot_group, replicate) %>%
  summarise(gene_counts = sum(count)) %>%
  ungroup() %>%
  ggplot(aes(x = plot_group, y = gene_counts, shape = replicate)) +
  facet_wrap("~ gene_name", ncol = 1, scales = "free_y")+
  geom_point(position = position_dodge(width = 0.5), size = 3)

# For each isoform, what proportion of total expression across fractions originates from each fraction?
# Although have tried to normalise with DESeq, this could still partly be driven by gene expression diffs between fractions
# key here is to assess relative difference between the two isoforms

normed_counts_npc_long_3exts %>%
  filter(gene_name %in% c("ELK1")) %>%
  ggplot(aes(x = plot_group, y = rep_fraction, colour = event_label, shape = replicate)) +
  facet_wrap("~ gene_name", ncol = 1) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  theme_bw(base_size = 14) +
  labs(x = "Fraction",
       y = "Fraction of total expression") +
  theme(axis.text.x = element_text(angle = 90))


# plot relative difference in proportion between the two isoforms in each replicate?
# i.e. each fraction + replicate, normalise to propn of total expression for proximal ('annotated') UTR
# retains/not disrupted by profile of relative GE between fractions
# But allows to compare across fractions the relative enrichment/depletion of each isoform
# Goal here is to look for preferential recruitment of the cryptic isoform relative to the control isoform


# What if just plot PPAU (i.e. percent of total iso expression) within sample for each fraction?
# Tells us of total ELK1 RNA in that fraction, this proportion of it is the cryptic
# Compare the proportion across fractions, relatively higher proportion of total RNA within 1 fraction suggests preferential recruitment
# BUT, PPAU normalised to total abundance in that fraction. not comparing between the two isoforms?


