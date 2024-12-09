library(tidyverse)
set.seed(123)

cryptics_bed <- read_tsv("processed/zeng_2024/supplementary_s5.cryptic_pas.all.bed", col_names = c("chr", "start", "end", "name", "score", "strand"))
cryptics_bed

cryptics_bed <- cryptics_bed %>%
  separate(name, into = c("APA_ID", "region", "gene_name", "pas_usage_control", "pas_usage_kd"), sep = "\\|", remove = F, convert = T)


cryptics_bed <- cryptics_bed %>%
  mutate(delta_pas_usage = pas_usage_kd - pas_usage_control,
         quantile_delta_pas_usage = ntile(delta_pas_usage, 4))

# counts for each event type and quantile of delta usage
quantile_counts <- cryptics_bed %>%
  count(region, quantile_delta_pas_usage)


# sample without replacement n rows for each event type and usage quantile
cryptics_bed_sampled <- cryptics_bed %>%
  group_by(region, quantile_delta_pas_usage) %>%
  slice_sample(n = 3, replace = FALSE) %>%
  ungroup()

# repeat but remove those with no assigned gene name
cryptics_bed_sampled_na <- cryptics_bed %>%
  filter(gene_name != "nan") %>%
  group_by(region, quantile_delta_pas_usage) %>%
  slice_sample(n = 3, replace = FALSE) %>%
  ungroup()


# write sampled dfs to tsv 
write_tsv(cryptics_bed_sampled, file.path("processed", "zeng_2024", "sampled_events_region_usage.all.cryptic_pas.tsv"))
write_tsv(cryptics_bed_sampled_na, file.path("processed", "zeng_2024", "sampled_events_region_usage.rm_na_genename.cryptic_pas.tsv"))

# Summarise global condition-wise usages and deltas
cryptics_exprn_summary <- cryptics_bed %>%
  select(pas_usage_control, pas_usage_kd, delta_pas_usage) %>%
  summary() %>%
  as.data.frame() %>%
  as_tibble() %>%
  select(-Var1, column_name = Var2, statistic = Freq)

# counts of events in each region
cryptics_region_counts <- cryptics_bed %>%
  count(region, sort = T)


write_tsv(cryptics_exprn_summary, file.path("processed", "zeng_2024", "expression_summary.all.cryptic_pas.tsv"))
write_tsv(cryptics_region_counts, file.path("processed", "zeng_2024", "event_counts.region.cryptic_pas.tsv"))


