library(tidyverse)

#' read in and combine tables containing per-sample PATR statistics 
process_combine_tbls <- function(stats_path, stbl_path, overlap_path, 
                                 stbl_cols_to_drop = c("strandedness", "bam")) {
  
  # read in per sample 'stats' table containing PATR counts passing filters
  # combine with sample_name metadata from the pipeline run  
  nygc_patr_stats <- read_tsv(stats_path) %>%
    left_join(
      read_csv(stbl_path) %>% select(-any_of(stbl_cols_to_drop)),
      by = "sample_name"
    )
  
  # Add in statistics of PAS cluster overlap with PolyASite db
  # Output a single combined tibble (1 row per sample_name)
  read_tsv(overlap_path) %>%
    mutate(file_name = str_remove_all(file_name, ".polya_clusters.bed$")) %>%
    rename(sample_name = file_name) %>%
    left_join(nygc_patr_stats, by = "sample_name") %>%
    mutate(
      overlap_patr_percentage = (overlapping_score_sum / patr_count) * 100,
      overlap_total_percentage = (overlapping_score_sum / primary_count) * 100,
      overlap_mean_count = overlapping_score_sum / overlapping
    ) %>%
    relocate(starts_with("overlap_"), .after = patr_per_million)
}

####


nygc_patr_stats <- process_combine_tbls(
  "data/nygc/patrs/2025-01-18_polya_softclips/pas_clusters/per_sample/all_samples.patrs.valid_stats.tsv",
  "data/nygc/patrs/2025-01-18_polya_softclips/patrs.nygc_sample_table.csv",
  "data/nygc/patrs/2025-01-18.two_class_simple.patr_clusters.polyasite_overlap.wcounts.tsv"
)

kd_patr_stats <- process_combine_tbls("data/nygc/patrs/tdp_ko_collection/pas_clusters/per_sample/all_samples.patrs.valid_stats.tsv",
                                      "data/nygc/patrs/tdp_ko_collection/combined_sample_table.csv",
                                      "data/nygc/patrs/tdp_ko_collection.two_class_simple.patr_clusters.polyasite_overlap.wcounts.tsv")


slamseq_patr_stats <- process_combine_tbls("data/nygc/patrs/slamseq_second_i3_cortical/pas_clusters/per_sample/all_samples.patrs.valid_stats.tsv",
                                           "data/nygc/patrs/slamseq_second_i3_cortical/i3cortical_slamseq_sample_table.csv",
                                           "data/nygc/patrs/slamseq_second_i3_cortical.two_class_simple.patr_clusters.polyasite_overlap.wcounts.tsv"
                                           )
outdir <- "processed/nygc/patrs/"

## add metadat fro the nygc and kd datasets

nygc_metadata <- read_tsv("data/nygc/NYGC_all_RNA_samples_support.tsv")
nygc_patr_stats <- left_join(nygc_patr_stats, select(nygc_metadata, sample_name = sample, platform, prep), by = "sample_name")

# KD collection is a bit more involved...
kd_metadata <- read_tsv("../misc/data/2023-11-22_paper_tdp43_collection_library_statistics.tsv")

# before joining/adding in library stats, need a unified dataset column
meta_names <- kd_metadata %>%
  filter(experiment_name != "Klim i3 Motor") %>%
  pull(experiment_name) %>%
  sort()

stbl_names <- sort(unique(kd_patr_stats$dataset))
df_names <- tibble(dataset = stbl_names, experiment_name = meta_names)

# prep kd_metadata to add to stats df
kd_metadata <- kd_metadata %>%
  # remove Klim...
  inner_join(df_names, by = "experiment_name") %>%
  select(dataset, experiment_name, mean_aligned_read_length, mean_library_size)

kd_patr_stats <- left_join(kd_patr_stats, kd_metadata, by = "dataset")


# generate combined df
comb_patr_stats <- bind_rows(slamseq = slamseq_patr_stats, nygc = nygc_patr_stats, kds = kd_patr_stats, .id = "origin")
comb_patr_stats <- comb_patr_stats %>%
  mutate(plot_origin = case_when(origin == "nygc" ~ "NYGC",
                                 origin == "kds" ~ "Cells",
                                 T ~ "SLAM-seq"
                                 )
         )

### PLOTS

kd_nygc_db_overlap_boxplot <- comb_patr_stats %>%
  # filter out the outlier sample
  filter(overlap_total_percentage < 0.2,
         plot_origin != "SLAM-seq") %>%
  mutate(plot_alpha = if_else(origin == "kds", 0.5, 0.01)) %>%
  ggplot(aes(x = plot_origin, y = overlap_total_percentage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(alpha = plot_alpha), position = position_jitter(width = 0.35, seed = 123),show.legend = F, size = 1
              #alpha = 0.1
              ) +
  labs(y = "PATRs overlapping with PolyASite\n(% library)",
       x = "Dataset") +
  theme_bw(base_size = 14)

kd_nygc_db_overlap_boxplot

# colour by prep for NYGC data
comb_patr_stats %>%
  # filter out the outlier sample
  filter(overlap_total_percentage < 0.2) %>%
  ggplot(aes(x = origin, y = overlap_total_percentage, colour = platform)) +
  geom_boxplot(position = "dodge", outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 123), alpha = 0.25) +
  labs(title = "Proportion of library that are PATRs (PolyASite PAS)",
       y = "PATRs overlapping with PolyASite (% library)") +
  theme_bw(base_size = 14)


# unsure of status of Brown Sh & SKs - what happens if remove?
comb_patr_stats %>%
  filter(!dataset %in% c("brown_shsy5y","brown_sknbe2")) %>%
  filter(overlap_total_percentage < 0.2) %>%
  ggplot(aes(x = origin, y = overlap_total_percentage, colour = platform)) +
  geom_boxplot(position = "dodge", outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 123), alpha = 0.25) +
  labs(title = "Proportion of library that are PATRs (PolyASite PAS)",
       y = "PATRs overlapping with PolyASite (% library)") +
  theme_bw(base_size = 14)

#  average support for PATR PAS overlapping PolyASite?

comb_patr_stats %>%
  # filter(overlap_total_percentage < 0.2) %>%
  ggplot(aes(x = plot_origin, y = overlap_mean_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.3, seed = 123), alpha = 0.25) +
  scale_y_continuous(limits = c(1,11),
                     breaks = seq(1,11)
  ) +
  theme_bw(base_size = 14) +
  labs(y = "Mean PATR count for detected PolyASite PAS",
       x = "Dataset")

comb_patr_stats %>%
    filter(!dataset %in% c("brown_shsy5y","brown_sknbe2")) %>%
  # filter(overlap_total_percentage < 0.2) %>%
  ggplot(aes(x = origin, y = overlap_mean_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.3, seed = 123), alpha = 0.25) +
  scale_y_continuous(limits = c(1,11),
                     breaks = seq(1,11)
  ) +
  theme_bw(base_size = 14) +
  labs(title = "PATR read support for detected PolyASite PAS")


# rRNA depletion - similar proportion of library that is PATRs for cells and tissue, with depletion vs polyA+ (even stronger with ambiguous Brown SH & SKs)
# Average Read support for annotated PAS is similar for cells vs tissue - suspect this is a cell type proportion issue?

# Sanity check - what factors are associated with read detection?

# read length

comb_patr_stats %>%
  filter(origin == "kds") %>%
  mutate(plot_len = fct_inseq(as.character(mean_aligned_read_length))) %>%
  filter(overlap_total_percentage < 0.2) %>%
  ggplot(aes(x = plot_len, y = overlap_total_percentage, colour = experiment_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 123))

comb_patr_stats %>%
  filter(origin == "kds") %>%
  filter(!dataset %in% c("brown_shsy5y","brown_sknbe2")) %>%
  mutate(plot_len = fct_inseq(as.character(mean_aligned_read_length))) %>%
  filter(overlap_total_percentage < 0.2) %>%
  ggplot(aes(x = plot_len, y = overlap_total_percentage, colour = experiment_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 123))


comb_patr_stats %>%
  filter(origin == "nygc") %>%
  mutate(mean_aligned_read_length = if_else(platform == "HiSeq 2500", 125, 100),
         plot_len = fct_inseq(as.character(mean_aligned_read_length))) %>%
  ggplot(aes(x = plot_len, y = overlap_total_percentage, colour = platform)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 123))


# read depth
comb_patr_stats %>%
  filter(origin == "kds") %>%
  filter(overlap_total_percentage < 0.2) %>%
  mutate(plot_libsize = primary_count / 1e6) %>%
  ggplot(aes(x = plot_libsize, y = overlap_total_percentage, colour = experiment_name)) +
  geom_point() +
  geom_smooth(mapping = aes(x = plot_libsize, y = overlap_total_percentage),inherit.aes = F,
                            method = "lm")
  scale_x_continuous(limits = c(0,250)) +
  scale_y_continuous(limits = c(0, 0.07))

comb_patr_stats %>%
  filter(origin == "kds") %>%
  filter(overlap_total_percentage < 0.2) %>%
  mutate(plot_libsize = primary_count / 1e6) %>%
  ggplot(aes(x = plot_libsize, y = overlap_total_percentage, colour = experiment_name)) +
  facet_wrap(~ experiment_name, scales = "free") +
  geom_point() +
  # scale_x_continuous(limits = c(0,250)) +
  # scale_y_continuous(limits = c(0, 0.07)) +
  geom_smooth(method = "lm")

comb_patr_stats %>%
  filter(origin == "slamseq") %>%
  filter(overlap_total_percentage < 0.2) %>%
  mutate(plot_libsize = primary_count / 1e6) %>%
  ggplot(aes(x = plot_libsize, y = overlap_total_percentage)) +
  geom_point() +
  # scale_x_continuous(limits = c(0,250)) +
  # scale_y_continuous(limits = c(0, 0.07)) +
  geom_smooth(method = "lm")


comb_patr_stats %>%
  filter(origin == "nygc") %>%
  mutate(plot_libsize = primary_count / 1e6) %>%
  ggplot(aes(x = plot_libsize, y = overlap_total_percentage)) +
  facet_wrap(~ platform) +
  geom_point() +
  scale_x_continuous(limits = c(0,250)) +
  scale_y_continuous(limits = c(0, 0.015)) +
  geom_smooth(method = "lm")

comb_patr_stats %>%
  filter(overlap_total_percentage < 0.2) %>%
  mutate(plot_libsize = primary_count / 1e6) %>%
  ggplot(aes(x = plot_libsize, y = overlap_total_percentage, colour = origin)) +
  # facet_wrap(~ experiment_name, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm")
  

### SAVING TO DISK

if (!dir.exists(outdir)) {dir.create(outdir)}

ggsave(filename = "2025-01-24_kds_nygc.patr_polysite_overlap_percentage.boxplot.png",
       plot = kd_nygc_db_overlap_boxplot,
       device = "png",
       path = outdir,
       height = 100,
       width = 100,
       units = "mm",
       dpi = "retina"
       )

ggsave(filename = "2025-01-24_kds_nygc.patr_polysite_overlap_percentage.boxplot.png",
       plot = kd_nygc_db_overlap_boxplot,
       device = "pdf",
       path = outdir,
       height = 100,
       width = 100,
       units = "mm",
       dpi = "retina"
)
