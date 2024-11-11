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
  theme()


# plot relative difference in proportion between the two isoforms in each replicate?
# i.e. each fraction + replicate, normalise to propn of total expression for proximal ('annotated') UTR
# retains/not disrupted by profile of relative GE between fractions
# But allows to compare across fractions the relative enrichment/depletion of each isoform
# Goal here is to look for preferential recruitment of the cryptic isoform relative to the control isoform


# What if just plot PPAU (i.e. percent of total iso expression) within sample for each fraction?
# Tells us of total ELK1 RNA in that fraction, this proportion of it is the cryptic
# Compare the proportion across fractions, relatively higher proportion of total RNA within 1 fraction suggests preferential recruitment
# BUT, PPAU normalised to total abundance in that fraction. not comparing between the two isoforms?


# TPMs & PPAUs

ppau <- read_tsv("processed/fracseq/2024-04-30_summarised_pas.ppau.tsv")
# add gene_name
ppau <- left_join(ppau, le2name, by = "le_id") %>%
  relocate(gene_name, .after = gene_id)

ppau_long <- pivot_longer(ppau,
             cols = -all_of(c("le_id", "gene_id", "gene_name")),
             names_to = "sample_name",
             values_to = "ppau")

ppau_3exts <- filter(ppau_long, gene_name %in% c("ELK1", "SIX3", "TLX1"))

# extract fraction, replicate number and cell type from sample_names
ppau_3exts <- ppau_3exts %>%
  mutate(sample_name_clean = str_remove_all(sample_name, "^SRR[0-9]*_"),
         replicate = str_extract_all(sample_name_clean, "rep[0-9]$", simplify = T),
         cell_type = str_split_i(sample_name_clean, "_", 1),
         group = case_when(str_detect(sample_name_clean, "cytosol") ~ "cytosol",
                           str_detect(sample_name_clean, "monosome") ~ "monosome",
                           str_detect(sample_name_clean, "2-4_ribosomes") ~ "2-4_ribosomes",
                           str_detect(sample_name_clean, "5_ribosomes") ~ "5_ribosomes"
                           ),
         group = factor(group, levels = c("cytosol",
                                          "monosome",
                                          "2-4_ribosomes",
                                          "5_ribosomes")
                        )
         ) %>%
  relocate(ppau, .after = everything())

ppau_3exts %>%
  filter(gene_name == "ELK1",
         str_ends(le_id, "_2")) %>%
  mutate(ppau = ppau * 100) %>%
  ggplot(aes(x = group, y = ppau, group = replicate, shape = replicate)) +
  facet_wrap("~ cell_type") +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(0,20, 2)
                     ) +
  labs(x = "Fraction",
       y = "Cryptic PAS usage (%)",
       shape = "Replicate") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# Normalise to cytosol usage
ppau_3exts_npc_cyto_normed <- ppau_3exts %>%
  filter(cell_type == "NPC") %>%
  pivot_wider(values_from = ppau, names_from = group, names_prefix = "ppau_") %>%
  filter(gene_name == "ELK1",
         str_ends(le_id, "_2")) %>%
  group_by(replicate) %>%
  # for each replicate, add default column containing cytosol usage
  mutate(all_ppau_cytosol = sort(unique(ppau_cytosol), na.last = NA)) %>%
  ungroup() %>%
  # each column, divide against cytosol usage
  mutate(across(starts_with("ppau_"), ~ .x / all_ppau_cytosol)) %>%
  arrange(replicate, sample_name_clean) %>% 
  select(-all_ppau_cytosol) %>%
  pivot_longer(cols = starts_with("ppau_"), names_prefix = "ppau_", names_to = "group", values_to = "cytosol_normed_ppau") %>%
  # introduces NAs for fraction + other fraction combos
  drop_na(cytosol_normed_ppau)


ppau_3exts_npc_cyto_normed %>%
  mutate(cytosol_normed_ppau = cytosol_normed_ppau * 100) %>%
  mutate(group = case_when(str_detect(sample_name_clean, "cytosol") ~ "cytosol",
                    str_detect(sample_name_clean, "monosome") ~ "monosome",
                    str_detect(sample_name_clean, "2-4_ribosomes") ~ "2-4_ribosomes",
                    str_detect(sample_name_clean, "5_ribosomes") ~ "5_ribosomes"),
         group = factor(group, levels = c("cytosol",
                                          "monosome",
                                          "2-4_ribosomes",
                                          "5_ribosomes")
                        )
         ) %>%
  ggplot(aes(x = group,
             y = cytosol_normed_ppau,
             group = replicate,
             shape = replicate)) +
  facet_wrap("~ cell_type") +
  geom_line() +
  geom_point() +
  labs(x = "Fraction",
       y = "Cytosol normalised PAS usage (%)",
       shape = "Replicate") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# bar plot of all fractions
ppau_3exts_npc_cyto_normed %>%
  mutate(cytosol_normed_ppau = cytosol_normed_ppau * 100) %>%
  mutate(group = case_when(str_detect(sample_name_clean, "cytosol") ~ "cytosol",
                           str_detect(sample_name_clean, "monosome") ~ "monosome",
                           str_detect(sample_name_clean, "2-4_ribosomes") ~ "2-4_ribosomes",
                           str_detect(sample_name_clean, "5_ribosomes") ~ "5_ribosomes"),
         group = factor(group, levels = c("cytosol",
                                          "monosome",
                                          "2-4_ribosomes",
                                          "5_ribosomes"),
         ),
         replicate = str_remove_all(replicate, "^rep")
         ) %>%
  ggplot(aes(x = replicate,
             y = cytosol_normed_ppau,
             fill = group,
             )) +
  geom_col(position = "dodge") +
  labs(title = "ELK1 cryptic PAS",
       subtitle = "Individual fractions",
       x = "Replicate",
       y = "Cytosol normalised PAS usage (%)",
       fill = "Fraction") +
  scale_y_continuous(limits = c(0,200),
                     breaks = seq(0,200,25)) +
  scale_fill_brewer(type = "qual",palette = "Paired") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# plot individual fractions on x axis, summarising replicates with a bar and jittering individual replicates

# clean up names for plotting
ppau_3exts_npc_cyto_normed_plot_df <- ppau_3exts_npc_cyto_normed %>%
  # mutate(cytosol_normed_ppau = cytosol_normed_ppau * 100) %>%
  mutate(group = case_when(str_detect(sample_name_clean, "cytosol") ~ "Cytosol",
                           str_detect(sample_name_clean, "monosome") ~ "Monosome",
                           str_detect(sample_name_clean, "2-4_ribosomes") ~ "Light polysome",
                           str_detect(sample_name_clean, "5_ribosomes") ~ "Heavy polysome"),
         group = factor(group, levels = c("Cytosol",
                                          "Monosome",
                                          "Light polysome",
                                          "Heavy polysome"),
         ),
         replicate = str_remove_all(replicate, "^rep")
  )

# plot
npc_cyto_normed_elk1_bar <- ggpubr::ggbarplot(ppau_3exts_npc_cyto_normed_plot_df,
                  x = "group",
                  y = "cytosol_normed_ppau",
                  add = c("median"),
                  ggtheme = theme_bw(base_size = 14)
                  ) +
  geom_point(data = ppau_3exts_npc_cyto_normed_plot_df,
             aes(x = group,
                 y = cytosol_normed_ppau,
                 shape = replicate),
             position = position_dodge(width = 0.4),
             size = 3
             ) +
  theme(legend.position = "top") +
  labs(x = "Fraction",
       y = "Cytosol-normalised ELK1 cryptic PAS usage",
       shape = "Replicate")
  
npc_cyto_normed_elk1_bar
# ppau_3exts_npc_cyto_normed_plot_df %>%
#   ggplot(aes(x = group,
#              y = cytosol_normed_ppau,
#              fill = group,
#   )) +
#   geom_col(position = "dodge") +
#   labs(title = "ELK1 cryptic PAS",
#        subtitle = "Individual fractions",
#        x = "Replicate",
#        y = "Cytosol normalised PAS usage (%)",
#        fill = "Fraction") +
#   scale_y_continuous(limits = c(0,2),
#                      breaks = seq(0,2,0.25)) +
#   scale_fill_brewer(type = "qual",palette = "Paired") +
#   theme_bw(base_size = 14) +
#   theme(legend.position = "bottom")



#
# group all 'some' fractions together vs cytosol
#

tpm <- read_tsv("processed/fracseq/2024-04-30_summarised_pas.tpm.tsv")
# add gene_name
tpm <- left_join(tpm, le2name, by = "le_id") %>%
  relocate(gene_name, .after = gene_id)

tpm_long <- pivot_longer(tpm,
                          cols = -all_of(c("le_id", "gene_id", "gene_name")),
                          names_to = "sample_name",
                          values_to = "tpm")

tpm_long_3exts <- filter(tpm_long, gene_name %in% c("ELK1", "SIX3", "TLX1")) %>%
  mutate(sample_name_clean = str_remove_all(sample_name, "^SRR[0-9]*_"),
         replicate = str_extract_all(sample_name_clean, "rep[0-9]$", simplify = T),
         cell_type = str_split_i(sample_name_clean, "_", 1),
         fraction = case_when(str_detect(sample_name_clean, "cytosol") ~ "cytosol",
                              str_detect(sample_name_clean, "monosome") ~ "monosome",
                              str_detect(sample_name_clean, "2-4_ribosomes") ~ "2-4_ribosomes",
                              str_detect(sample_name_clean, "5_ribosomes") ~ "5_ribosomes"),
         fraction = factor(fraction, levels = c("cytosol",
                                          "monosome",
                                          "2-4_ribosomes",
                                          "5_ribosomes")
                           )
  )

# pool the expression of isoforms across -some fractions and cytosol (per replicate)
sum_tpm_some_cyto_3exts <- tpm_long_3exts %>%
  mutate(group = if_else(fraction == "cytosol", "cytosol", "ribosome")) %>%
  group_by(cell_type, replicate, gene_name, le_id, group) %>%
  summarise(sum_tpm = sum(tpm)) %>%
  ungroup()

# calculate PPAu for each fraction (per replicate)
ppau_some_cyto_3exts <- sum_tpm_some_cyto_3exts %>%
  group_by(cell_type, replicate, group, gene_name) %>%
  mutate(ppau = (sum_tpm / sum(sum_tpm)) * 100) %>%
  ungroup()
  

# Normalise to cytosol usage for ELK1
ppau_npc_elk1_some_cyto_cyto_normed <- ppau_some_cyto_3exts %>%
  filter(cell_type == "NPC") %>%
  select(-sum_tpm) %>%
  pivot_wider(values_from = ppau, names_from = group, names_prefix = "ppau_") %>%
  filter(gene_name == "ELK1",
         str_ends(le_id, "_2")) %>%
  # each ppau column, divide against cytosol usage
  mutate(across(starts_with("ppau_"), ~ (.x / ppau_cytosol) * 100)) %>%
  arrange(replicate) %>% 
  pivot_longer(cols = starts_with("ppau_"), names_prefix = "ppau_", names_to = "group", values_to = "cytosol_normed_ppau")
  

# bar plot of cytosol normalised usage
ppau_npc_elk1_some_cyto_cyto_normed %>%
  ggplot(aes(x = replicate,
             y = cytosol_normed_ppau,
             fill = group,
             replicate)) +
  geom_col(position = "dodge") +
  labs(title = "ELK1 cryptic PAS",
       subtitle = "ribosome = monosome + polysome",
       x = "Replicate",
       y = "Cytosol normalised PAS usage (%)",
       fill = "Fraction") +
  scale_y_continuous(limits = c(0,180),
                     breaks = seq(0,200,20)) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# same, but plotting fractioons on x
ppau_npc_elk1_some_cyto_cyto_normed %>%
  ggplot(aes(fill = replicate,
             y = cytosol_normed_ppau,
             x = group,
             replicate)) +
  geom_col(position = "dodge") +
  labs(title = "ELK1 cryptic PAS",
       subtitle = "ribosome = monosome + polysome",
       x = "Replicate",
       y = "Cytosol normalised PAS usage (%)",
       fill = "Fraction") +
  scale_y_continuous(limits = c(0,160),
                     breaks = seq(0,200,20)) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")


# Repeat but for cytosol, monosome + 2-4/5+ polysomes
sum_tpm_polysome_cyto_3exts <- tpm_long_3exts %>%
  mutate(group = if_else(fraction %in% c("2-4_ribosomes", "5_ribosomes"),
                         "polysome",
                         fraction)) %>%
  group_by(cell_type, replicate, gene_name, le_id, group) %>%
  summarise(sum_tpm = sum(tpm)) %>%
  ungroup()

# calculate PPAu for each fraction (per replicate)
ppau_polysome_cyto_3exts <- sum_tpm_polysome_cyto_3exts %>%
  group_by(cell_type, replicate, group, gene_name) %>%
  mutate(ppau = (sum_tpm / sum(sum_tpm)) * 100) %>%
  ungroup()

# Normalise to cytosol usage for ELK1
ppau_npc_elk1_polysome_cyto_cyto_normed <- ppau_polysome_cyto_3exts %>%
  filter(cell_type == "NPC") %>%
  select(-sum_tpm) %>%
  pivot_wider(values_from = ppau, names_from = group, names_prefix = "ppau_") %>%
  filter(gene_name == "ELK1",
         str_ends(le_id, "_2")) %>%
  # each ppau column, divide against cytosol usage
  mutate(across(starts_with("ppau_"), ~ (.x / ppau_cytosol) * 100)) %>%
  arrange(replicate) %>% 
  pivot_longer(cols = starts_with("ppau_"), names_prefix = "ppau_", names_to = "group", values_to = "cytosol_normed_ppau")

ppau_npc_elk1_polysome_cyto_cyto_normed

# bar plot of polysome pooled, 
ppau_npc_elk1_polysome_cyto_cyto_normed %>%
  ggplot(aes(x = replicate,
             y = cytosol_normed_ppau,
             fill = group,
             replicate)) +
  geom_col(position = "dodge") +
  labs(title = "ELK1 cryptic PAS",
       subtitle = "polysome = 2-4 & 5+ ribsomes pooled",
       x = "Replicate",
       y = "Cytosol normalised PAS usage (%)",
       fill = "Fraction") +
  scale_y_continuous(limits = c(0,200),
                     breaks = seq(0,200,25)) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")


# same plot, but plotting fractions on x axis
ppau_npc_elk1_polysome_cyto_cyto_normed %>%
  ggplot(aes(x = group,
             y = cytosol_normed_ppau,
             fill = replicate,
             )) +
  geom_col(position = "dodge") +
  labs(title = "ELK1 cryptic PAS",
       subtitle = "polysome = 2-4 & 5+ ribsomes pooled",
       x = "Replicate",
       y = "Cytosol normalised PAS usage (%)",
       fill = "Fraction") +
  scale_y_continuous(limits = c(0,200),
                     breaks = seq(0,200,25)) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# median bar plot with replicates labelled as different points
npc_cyto_normed_elk1_pooled_bar <- ggpubr::ggbarplot(ppau_npc_elk1_polysome_cyto_cyto_normed,
                  x = "group",
                  y = "cytosol_normed_ppau",
                  add = c("median"),
                  ggtheme = theme_bw(base_size = 14)
) +
  geom_point(data = ppau_npc_elk1_polysome_cyto_cyto_normed,
             aes(x = group,
                 y = cytosol_normed_ppau,
                 shape = replicate),
             position = position_dodge(width = 0.4),
             size = 3
  ) +
  theme(legend.position = "top") +
  labs(x = "Fraction",
       y = "Cytosol-normalised ELK1 cryptic PAS usage",
       shape = "Replicate")

npc_cyto_normed_elk1_pooled_bar

# Save plots to disk
outdir <- "processed/fracseq"
if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

ggsave(file.path(outdir,
                 "2024-11-11_fracseq.npc.elk1.cytosol_normed.median_bar.all_fractions.png"),
       plot = npc_cyto_normed_elk1_bar, width = 150, height = 150, units = "mm", dpi = "retina")


ggsave(file.path(outdir,
                 "2024-11-11_fracseq.npc.elk1.cytosol_normed.median_bar.polysome_pooled.png"),
       plot = npc_cyto_normed_elk1_pooled_bar, width = 150, height = 150, units = "mm", dpi = "retina")


