library(tidyverse)

# read in formatted NYGC metdata used for analysis, plotting
nygc_metadata <- read_tsv("processed/nygc/NYGC_all_RNA_samples_support_formatted.tsv")

# raw metdata (for )
raw_nygc_metadata <- read_tsv("data/nygc/NYGC_all_RNA_samples_support.tsv")

# extract gender for each individual
indiv_gender <- raw_nygc_metadata %>%
  select(individual, sex) %>%
  distinct(.keep_all = T)

# add a simple categorisation of disease
nygc_metadata <- nygc_metadata %>%
  mutate(simple_disease = case_when(str_detect(disease,"^ALS") ~ "ALS",
                                    str_detect(disease, "^FTD") ~ "FTD",
                                    TRUE ~ "Control"))

# add sex to metadata
nygc_metadata <- left_join(nygc_metadata, indiv_gender, by = "individual")

# Get a table of the number of samples and individuals across dataset
n_individuals <- n_distinct(nygc_metadata$individual)
n_samples <- n_distinct(nygc_metadata$sample)
dataset_wide_counts <- tibble( type = c("individuals", "samples"), n = c(n_individuals, n_samples))

# number of indivuals split by sex
nygc_metadata %>%
  distinct(individual, .keep_all = T) %>%
  count(sex)

# number of indivuals for each disease subcategory
nygc_metadata %>%
  distinct(individual, .keep_all = T) %>%
  count(simple_disease)

# number of indivuals for each disease subcategory and sex
nygc_metadata %>%
  distinct(individual, .keep_all = T) %>%
  count(simple_disease, sex)


nygc_metadata %>%
  distinct(individual, .keep_all = T) %>%
  group_by(simple_disease) %>%
  summarise(median_age = median(age, na.rm=T),
            iqr_age = IQR(age, na.rm = T))



# get a table with sample counts by TDP path subtype
dataset_wide_counts_path <- count(nygc_metadata, tdp_path)

# Number of samples for each tissue, subdivided by disease subtype
tissue_subtype_counts <- count(nygc_metadata, tissue_clean, disease)
tissue_counts <- count(nygc_metadata, tissue_clean)


write_tsv(dataset_wide_counts, "processed/nygc/2023-12-20_nygc_dataset_wide_counts.tsv")
write_tsv(dataset_wide_counts_path, "processed/nygc/2023-12-20_nygc_dataset_wide_pathology_counts.tsv")
write_tsv(tissue_subtype_counts, "processed/nygc/2023-12-20_nygc_tissue_subtype_counts.tsv")
write_tsv(tissue_counts, "processed/nygc/2023-12-20_nygc_tissue_total_counts.tsv")

# repeat but filtering for disease tissues only
nygc_metadata_disease <- filter(nygc_metadata, disease_tissue)

# Get a table of the number of samples and individuals across dataset
n_individuals_disease <- n_distinct(nygc_metadata_disease$individual)
n_samples_disease <- n_distinct(nygc_metadata_disease$sample)
dataset_wide_counts_disease <- tibble(type = c("individuals", "samples"),
                                      n = c(n_individuals_disease, n_samples_disease))

# get a table with sample counts by TDP path subtype
dataset_wide_counts_path_disease <- count(nygc_metadata_disease, tdp_path)

# Number of samples for each tissue, subdivided by disease subtype
tissue_subtype_counts_disease <- count(nygc_metadata_disease, tissue_clean, disease)
tissue_counts_disease <- count(nygc_metadata_disease, tissue_clean)


write_tsv(dataset_wide_counts, "processed/nygc/2023-12-20_nygc_disease_counts.tsv")
write_tsv(dataset_wide_counts_path, "processed/nygc/2023-12-20_nygc_disease_pathology_counts.tsv")
write_tsv(tissue_subtype_counts, "processed/nygc/2023-12-20_nygc_disease_tissue_subtype_counts.tsv")
write_tsv(tissue_counts, "processed/nygc/2023-12-20_nygc_disease_tissue_total_counts.tsv")

### Supplementary table with NYGC sample characteristics


raw_nygc_metadata_filtered <- raw_nygc_metadata %>%
  left_join(distinct(nygc_metadata, sample, disease_tissue, tdp_path, disease, simple_disease), by = "sample", suffix = c(".raw", ".clean")) %>%
  # subset to analysed samples only
  filter(sample %in% nygc_metadata$sample) %>%
  # subset to disease tissues analysed
  filter(disease_tissue) 

# summarise distribution of numeric values
disease_metadata_summ_stats <- raw_nygc_metadata_filtered %>%
  group_by(disease.clean) %>%
  summarise(
    across(
      c(rin, pmi, age, onset), # Specify the columns you want to summarise
      list(
        median = ~median(.x, na.rm = TRUE),
        q1 = ~quantile(.x, 1/4, na.rm = T),
        q3 = ~quantile(.x, 3/4, na.rm = T),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE),
        missing = ~sum(is.na(.x))
      ),
      .names = "{.fn}_{.col}" # Specify the naming pattern for new columns
    )
  )

dplyr::last_dplyr_warnings()

disease_metadata_summ_stats <- disease_metadata_summ_stats %>%
  mutate(across(ends_with("rin"), ~ round(.x, 3)),
         across(ends_with("pmi"), ~ round(.x, 3))
         )

disease_metadata_summ_stats <- disease_metadata_summ_stats %>%
  mutate(
    rin_summary = glue::glue("{median_rin} ({q1_rin}, {q3_rin}, missing = {missing_rin})"),
    pmi_summary = glue::glue("{median_pmi} ({q1_pmi}, {q3_pmi}, missing = {missing_pmi})"),
    age_summary = glue::glue("{median_age} ({q1_age}, {q3_age}, missing = {missing_age})"),
    onset_summary = glue::glue("{median_onset} ({q1_onset}, {q3_onset}, missing = {missing_onset})")
  ) %>%
  select(disease.clean,ends_with("_summary")) %>%
  rename_with(
    .fn = ~ str_remove(.x, "_summary"), 
    .cols = ends_with("_summary")       
  )

disease_metadata_summ_stats

# calculate number of samples and individuals across disease types
disease_summ_counts <- raw_nygc_metadata_filtered %>%
  group_by(disease.clean) %>%
  summarise(n_individuals = n_distinct(individual),
            n_samples = n_distinct(sample)
            )
# individual counts split by sex
disease_summ_counts_sex <-  raw_nygc_metadata_filtered %>%
  distinct(disease.clean, individual, sex) %>%
  group_by(disease.clean) %>%
  summarise(n_male = sum(sex == "Male"),
            n_female = sum(sex == "Female")
            ) %>%
  ungroup()


# add counts to 
disease_metadata_summ_stats <- disease_summ_counts %>%
  left_join(disease_summ_counts_sex, by = "disease.clean", relationship = "one-to-one") %>%
  relocate(n_samples, .after = everything()) %>%
  left_join(disease_metadata_summ_stats, by = "disease.clean", relationship = "one-to-one")

# final clean up
disease_metadata_summ_stats <- disease_metadata_summ_stats %>%
  rename(subtype = disease.clean, age_of_onset = onset) %>%
  mutate(subtype = factor(subtype, levels = c("Control", "ALS-non-TDP", "ALS-TDP", "FTD-non-TDP", "FTD-TDP"))) %>%
  arrange(subtype)


write_tsv(disease_metadata_summ_stats, "processed/nygc/2024-11-27_nygc_metadata_summary.tdp_subtypes.tsv")
  




