library(tidyverse)

fish_counts <- read_tsv("data/fish_counts_cleaned.tsv")

# tidyup image & condition col (remove prefix)
fish_counts <- mutate(fish_counts, 
                      image = as.numeric(str_remove_all(image, "^CTRL|TDP")),
                      condition = str_remove_all(condition, "[0-9]")
                      )

# colnames(fish_counts)

# pivot longer the nuclear & extranuclear columns to cols of probe & count_type
fish_counts_longer <- pivot_longer(fish_counts,
             cols = ends_with("nuclear"),
             names_sep = "_",
             names_to = c("probe", "region"),
             values_to = "n_foci")

# Now calculate processed values
# 1. mean total foci counts per condition and replicate
# 2. ratios of extranuclear:nuclear foci counts per condition and replicate

# 10 images for each replicate
# sum the counts across images, normalise to sum of nuclei counts
fish_counts_mean <- fish_counts_longer %>%
  group_by(condition, replicate, probe) %>%
  summarise(sum_nuclei = sum(n_nuclei),
            sum_foci = sum(n_foci),
            mean_cell = sum_nuclei / sum_foci) %>%
  ungroup()

# for each probe and replicate, calculate the nuclear:extra-nuclear foci-count ratio
fish_counts_mean_region <- fish_counts_longer %>%
  group_by(condition, replicate, probe, region) %>%
  summarise(sum_foci = sum(n_foci)) %>%
  ungroup() %>%
  # convert foci counts to columns and calculate the ratios
  pivot_wider(names_from = region, values_from = sum_foci, names_prefix = "sum_foci_") %>%
  mutate(extranuc_nuc_ratio = sum_foci_extranuclear / sum_foci_nuclear,
         nuc_extranuc_ratio = sum_foci_nuclear / sum_foci_extranuclear)

# combine mean counts & ratios into a single df
fish_counts_mean <- left_join(fish_counts_mean, fish_counts_mean_region, by = c("condition", "replicate", "probe"))

# output tables to file
if (!dir.exists("processed")) {dir.create("processed")}

write_tsv(fish_counts_longer, "processed/2023-10-13_fish_counts_cleaned_longer.tsv", col_names = T)
write_tsv(fish_counts_mean, "processed/2023-10-13_fish_counts_processed.tsv", col_names = T)