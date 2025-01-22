library(tidyverse)

# table containing per-sample junction counts for all cryptic SJs
spliced_counts_ale <- read_tsv("processed/nygc/spliced_counts_ale_all.tsv")

# table containing 'detected' (nreads >=2) counts across TDP path status - useful for subsetting to relevant junctions
expression_by_pathology_ale <- read_tsv("processed/nygc/expression_by_pathology_ale_all.tsv")

# calculate path:not_path detection ratio (with view to selective most enriched/detected)
expression_by_pathology_ale <- expression_by_pathology_ale %>%
  mutate(path_fold_enrichment = (fraction_path + 0.001) / (fraction_not_path + 0.001)) %>%
  arrange(desc(selective), desc(path_fold_enrichment))

# pick selective SJs and 2 most enriched events (previously identified, also have high path detection between STMN2 and SYNJ2)
selected_ales <- expression_by_pathology_ale %>%
  filter(selective | gene_name %in% c("HS6ST3", "ARHGAP32"))

selected_spliced_counts_ale <- filter(spliced_counts_ale, paste_into_igv_junction %in% pull(selected_ales, paste_into_igv_junction))




###
stbl_cols_to_drop <- c("strandedness", "bam")
nygc_patr_stats <- read_tsv("data/nygc/patrs/2025-01-18_polya_softclips/pas_clusters/per_sample/all_samples.patrs.valid_stats.tsv")
nygc_patr_stbl <- read_csv("data/nygc/patrs/2025-01-18_polya_softclips/patrs.nygc_sample_table.csv")
nygc_patr_stats <- nygc_patr_stats %>%
  left_join(select(nygc_patr_stbl, -any_of(stbl_cols_to_drop)),by = "sample_name")

# 
nygc_metadata <- read_tsv("data/nygc/NYGC_all_RNA_samples_support.tsv")
nygc_patr_stats <- left_join(nygc_patr_stats, select(nygc_metadata, sample_name = sample, platform, prep), by = "sample_name")
count(nygc_patr_stats, platform, prep)

kd_patr_stats <- read_tsv("data/nygc/patrs/tdp_ko_collection/pas_clusters/per_sample/all_samples.patrs.valid_stats.tsv")
kd_patr_stbl <- read_csv("data/nygc/patrs/tdp_ko_collection/combined_sample_table.csv")
kd_patr_stats <- kd_patr_stats %>%
  left_join(select(kd_patr_stbl, -any_of(stbl_cols_to_drop)), by = "sample_name")

kd_metadata <- read_tsv("../misc/data/2023-11-22_paper_tdp43_collection_library_statistics.tsv")

# before joining/adding in library stats, need a unified dataset column
meta_names <- kd_metadata %>%
  filter(experiment_name != "Klim i3 Motor") %>%
  pull(experiment_name) %>%
  sort()

stbl_names <- sort(unique(kd_patr_stbl$dataset))
df_names <- tibble(dataset = stbl_names, experiment_name = meta_names)

# prep kd_metadata to add to stats df
kd_metadata <- kd_metadata %>%
  # remove Klim...
  inner_join(df_names, by = "experiment_name") %>%
  select(dataset, experiment_name, mean_aligned_read_length, mean_library_size)

kd_patr_stats <- left_join(kd_patr_stats, kd_metadata, by = "dataset")

comb_patr_stats <- bind_rows(nygc = nygc_patr_stats, kds = kd_patr_stats, .id = "origin")

# initial exploratory plot, labelling outlier values. Start with PATR per million
#
comb_patr_stats %>%
  mutate(plot_lab = if_else(patr_percentage > 0.2, sample_name, "")) %>%
  ggplot(aes(x = origin, y = log2(patr_per_million))) +
  geom_boxplot() +
  geom_text(aes(label = plot_lab), nudge_x = 0.1)

kd_patr_stats %>%
  filter(sample_name == "TDP_5")

# plot percentage values in hope can keep to non-log scale (easier to interpret) 
# for plotting purposes shrink the combined boxplot with jitter
comb_patr_stats %>%
  mutate(plot_lab = if_else(patr_percentage > 0.2, paste(sample_name, "- true value = ", round(patr_percentage, 3), "%", sep = " "), ""),
         plot_patr_percentage = if_else(patr_percentage > 0.2, 0.2, patr_percentage)) %>%
  ggplot(aes(x = origin, y = plot_patr_percentage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.25, seed = 123), alpha = 0.5)+
    geom_text(aes(label = plot_lab),nudge_x = 0.25) +
  labs(x = "Dataset",
       y = "% library that are valid PATRs") +
  theme_bw(base_size = 14)

comb_patr_stats %>%
  filter(patr_percentage < 0.2) %>% # remove outlier
  # mutate(plot_lab = if_else(patr_percentage > 0.2, paste(sample_name, "- true value = ", round(patr_percentage, 3), "%", sep = " "), ""),
  #        plot_patr_percentage = if_else(patr_percentage > 0.2, 0.2, patr_percentage)) %>%
  ggplot(aes(x = origin, y = patr_percentage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.25, seed = 123), alpha = 0.5)+
  scale_y_continuous(limits = c(0, 0.1),
                     breaks = seq(0, 0.2, 0.01)#seq(0, 0.2, 0.025)
  ) +
  labs(x = "Dataset",
       y = "% library that are valid PATRs") +
  theme_bw(base_size = 14)

# KD collection - not frequently higher than the NYGC data - should check if associated with dataset
# Certainly appears to be two groups for the NYGC data - what variables are associated?

kd_patr_stats %>%
  # handle pesky humphrey sample
  mutate(plot_lab = if_else(patr_percentage > 0.2, paste(sample_name, "- true value = ", round(patr_percentage, 3), "%", sep = " "), ""),
         plot_patr_percentage = if_else(patr_percentage > 0.2, 0.2, patr_percentage)) %>%
  ggplot(aes(x = dataset, y = plot_patr_percentage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.25, seed = 123), alpha = 0.5) +
  scale_y_continuous(limits = c(0, 0.2),
                     breaks = seq(0, 0.2, 0.01)#seq(0, 0.2, 0.025)
  ) +
  theme_bw(base_size = 14) +
  labs(y = "% library that are valid PATRs") +
  theme(axis.text.x = element_text(angle = 90))

kd_patr_stats %>%
  # handle pesky humphrey sample
  mutate(plot_lab = if_else(patr_percentage > 0.2, paste(sample_name, "- true value = ", round(patr_percentage, 3), "%", sep = " "), ""),
         plot_patr_percentage = if_else(patr_percentage > 0.2, 0.2, patr_percentage)) %>%
  ggplot(aes(x = dataset, y = plot_patr_percentage, colour = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 123), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 123), alpha = 0.5) +
  scale_y_continuous(limits = c(0, 0.2),
                     breaks = seq(0, 0.2, 0.01)#seq(0, 0.2, 0.025)
                     ) +
  theme_bw(base_size = 14) +
  labs(y = "% library that are valid PATRs") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "top")

# Driven by the high depth, high read length i3 cortical samples - as much as 5x the minimum by increasing the read length (I think doubling read length in seddighis vs brown)
# Pull in library parameters from supplemental table, and see if association with read-length, read depth etc.?

kd_patr_stats %>%
  # remove humphrey outlier
  filter(patr_percentage < 0.2) %>%
  mutate(mean_aligned_read_length = fct_inseq(as.character(mean_aligned_read_length)),
         # plot_patr_percentage = if_else(patr_percentage > 0.12, 0.12, patr_percentage)
         ) %>%
  ggplot(aes(x = mean_aligned_read_length, y = patr_percentage, colour = experiment_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 123), alpha = 0.75) +
  scale_y_continuous(limits = c(0, 0.1),
                     breaks = seq(0, 0.2, 0.01)#seq(0, 0.2, 0.025)
  ) +
  theme_bw(base_size = 14) +
  labs(y = "% library that are valid PATRs",
       x = "Average read length (nt)") +
  theme(legend.position = "top")



# NYGC inspection - 

# Inspect the association of PATR occurrence with available variables
nygc_patr_stats %>%
  filter(disease_tissue) %>%
  # minimal variables of interest
  select(sample_name, patr_percentage, region, disease, tissue_clean, tdp_path, platform, prep) %>%
  # Pivot to long format
  pivot_longer(
    cols = -all_of(c("sample_name", "patr_percentage")),
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(aes(x = value, y = patr_percentage)) +
  facet_wrap(~variable, scales = "free_x") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2, seed = 123), alpha = 0.25, size = 0.8) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90))

# clear separation on platform/prep - note that platform and prep are redundant, so repeat removing prep
nygc_patr_stats %>%
  filter(disease_tissue) %>%
  # minimal variables of interest
  select(sample_name, patr_percentage, region, disease, tissue_clean, tdp_path, platform) %>%
  # Pivot to long format
  pivot_longer(
    cols = -all_of(c("sample_name", "patr_percentage")),
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(aes(x = value, y = patr_percentage)) +
  facet_wrap(~variable, scales = "free_x") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2, seed = 123), alpha = 0.25, size = 0.8) +
  theme_bw(base_size = 14) +
  labs(subtitle = "HiSeq 2500 = 125bp and Manual KAPA Total RNA\nNovaSeq = 100bp and Automated KAPA Total RNA") +
  theme(axis.text.x = element_text(angle = 90))

# repeat, but colour other variables by platform 
nygc_patr_stats %>%
  filter(disease_tissue) %>%
  # minimal variables of interest
  select(sample_name, patr_percentage, region, disease, tissue_clean, tdp_path, platform) %>%
  # Pivot to long format
  pivot_longer(
    cols = -all_of(c("sample_name", "platform", "patr_percentage")),
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(aes(x = value, y = patr_percentage, colour = platform)) +
  facet_wrap(~variable, scales = "free_x") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 123),#position_jitter(width = 0.2, seed = 123), 
              alpha = 0.5, size = 0.8) +
  theme_bw(base_size = 14) +
  labs(title = "Sequencing platform and read length systematically affects PATR detection",
       subtitle = "HiSeq 2500 = 125bp and Manual KAPA Total RNA (n = 427)\nNovaSeq = 100bp and Automated KAPA Total RNA (n = 713)",
       x = NULL,
       y = "% library that are valid PATRs") +
  theme(axis.text.x = element_text(angle = 90))

# simple plot of just platform against each-other
nygc_patr_stats %>%
  filter(disease_tissue) %>%
  ggplot(aes(x = platform, y = patr_percentage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2, seed = 123), alpha = 0.25, size = 0.8) +
  scale_y_continuous(limits = c(0, 0.06),
                     breaks = seq(0, 0.6, 0.01)) +
  theme_bw(base_size = 14) +
  labs(title = "Sequencing platform and read length systematically affects PATR detection",
       subtitle = "HiSeq 2500 = 125bp and Manual KAPA Total RNA (n = 427)\nNovaSeq = 100bp and Automated KAPA Total RNA (n = 713)",
       x = "Platform",
       y = "% library that are valid PATRs")


nygc_patr_stats %>%
  filter(disease_tissue) %>%
  count(platform, sort = T) %>%
  mutate(frac = n / sum(n))


nygc_patr_stats %>%
  filter(disease_tissue) %>%
  group_by(platform) %>%
  summarise(median_perc = median(patr_percentage))

kd_patr_stats %>%
  group_by(dataset) %>%
  summarise(median_perc = median(patr_percentage)) %>%
  arrange(desc(median_perc))


### PolyASite overlap
nygc_pasite_overlap <- read_tsv("data/nygc/patrs/2025-01-18.two_class_simple.patr_clusters.polyasite_overlap.tsv")
kd_pasite_overlap <- read_tsv("data/nygc/patrs/tdp_ko_collection.two_class_simple.patr_clusters.polyasite_overlap.tsv")

nygc_pasite_overlap <- nygc_pasite_overlap %>%
  mutate(file_name = str_remove_all(file_name, ".polya_clusters.bed$")) %>%
  rename(sample_name = file_name)

kd_pasite_overlap <- kd_pasite_overlap %>%
  mutate(file_name = str_remove_all(file_name, ".polya_clusters.bed$")) %>%
  rename(sample_name = file_name)

comb_pasite_overlap <- bind_rows(nygc = nygc_pasite_overlap, kds = kd_pasite_overlap, .id = "origin")

comb_pasite_overlap %>%
  ggplot(aes(x = origin, y = percent_overlapping)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.1) +
  scale_y_continuous(limits = c(0,100))


comb_pasite_overlap %>%
  ggplot(aes(x = origin, y = percent_overlapping)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.3, seed = 123), alpha = 0.25) +
  scale_y_continuous(limits = c(0,100)) +
  labs(y = "% PATR PAS overlapping with PolyASite") +
    theme_bw(base_size = 14)

# No read filtering, systematically higher for KDs. Outliers in KDs are the brown et al SH & SK samples, where polyA+/rrna depletion is not directly reported and cannot be easily

# Simple - count PATRs overlapping with annotated PAS - compare proportion of library vs KDs?
# Need to modify shell script to report total count of PAS