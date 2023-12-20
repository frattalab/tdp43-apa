library(tidyverse)

# read in formatted NYGC metdata used for analysis, plotting
nygc_metadata <- read_tsv("processed/nygc/NYGC_all_RNA_samples_support_formatted.tsv")

# Get a table of the number of samples and individuals across dataset
n_individuals <- n_distinct(nygc_metadata$individual)
n_samples <- n_distinct(nygc_metadata$sample)
dataset_wide_counts <- tibble( type = c("individuals", "samples"), n = c(n_individuals, n_samples))

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
