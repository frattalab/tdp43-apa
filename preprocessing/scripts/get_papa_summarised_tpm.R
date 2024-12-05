library(tidyverse)


# output of get_num_pas.py - used to get le_ids of interest (as a vector)
le_ids <- read_tsv("processed/curation/cryptic_annot_comparison/2024-09-03_le_id_pas_counts.tsv", col_select = le_id) %>%
pull(le_id)

# top level directory storing PAPA outputs for all runs
# i.e. <papa_results_dir>/<dataset_name>/ --> PAPA outputs
papa_results_dir <- "data/PAPA/i3_cortical_upf1_zanovello_novel"

# directory storing PAPA sample tables (with same subdirectory name and structure as papa_results_dir)
# i.e. <papa_results_dir>/<dataset_name>/sample_table.csv
sample_tables_dir <- "data/PAPA/sample_tables"

# names of subdirectories corresponding to experiments where PATRs were extracted
# (produces matrix where splits are columns, access as a single vector for the one row
dataset_names <- str_split("brown_i3_cortical brown_shsy5y brown_skndz humphrey_i3_cortical seddighi_i3_cortical zanovello_i3_cortical_upf1_tdp_tdpkd_upf1ctl_vs_tdpctl_upf1ctl zanovello_shsy5y_chx_kd_ctl_vs_ctl_ctl zanovello_shsy5y_curve_0075 zanovello_skndz_curve_1", " ", simplify=T)[1, ]

# output directory for wide-format mean and median TPM matrices
outdir <- "processed/curation/cryptic_annot_comparison/"

# Generate vectors of paths to dataset-specific directories
papa_result_paths <- file.path(papa_results_dir, dataset_names)
sample_table_paths <- file.path(sample_tables_dir, dataset_names)

# Check if all generated directories exist
stopifnot(sapply(papa_result_paths, dir.exists))
stopifnot(sapply(sample_table_paths, dir.exists))

# Construct file paths to isoform-level TPM matrices
papa_result_files <- file.path(papa_result_paths, "differential_apa/summarised_pas_quantification.tpm.tsv") %>%
  set_names(dataset_names)

# Check if all constructed file paths exist
stopifnot(sapply(papa_result_files, file.exists))

# Locate sample_tables for each dataset
sample_table_files <- sapply(sample_table_paths, function(dir) {
  csv_files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  stopifnot(length(csv_files) == 1)  # Ensure there's exactly one .csv file per directory
  return(csv_files)
}, USE.NAMES = FALSE) %>%
  set_names(dataset_names)

# stopifnot(sapply(sample_table_files, file.exists))

# Extract KD sample names for each dataset

# Function to read in the CSV file, find the 2nd occurring value in the 'condition' column (i.e. TDP-43 KD),
# Return a vector of the corresponding sample_names.
extract_sample_names <- function(file_path) {
  # Read the CSV file
  df <- read_csv(file_path, show_col_types = F)
  
  # Find the 2nd occurring value in the 'condition' column
  second_condition <- unique(df$condition)[2]
  
  # Filter rows where 'condition' equals the 2nd occurring value
  filtered_df <- df %>% filter(condition == second_condition)
  
  # Return a vector of the 'sample_name' column
  sample_names <- filtered_df$sample_name
  
  return(sample_names)
}

# Iterate over the sample_table_files vector and apply the function using map from purrr
sample_name_vectors <- purrr::map(sample_table_files, extract_sample_names)
sample_name_vectors


# Function to extract specific columns from the TSV files based on sample names and le_ids
extract_columns_from_tsv <- function(sample_names, tsv_file) {
  # Generate a vector of column names to extract
  columns_to_extract <- c("le_id", sample_names)
  
  # Read in the TSV file, only extracting the specified columns
  df <- read_tsv(tsv_file, col_select = all_of(columns_to_extract), show_col_types = F)
  
  # Subset the data to include only rows where le_id is in the le_ids vector
  filtered_df <- df %>% filter(le_id %in% le_ids)
  
  return(filtered_df)
}

# Extract per-isoform TPM values across KD samples for each dataset
tpm_dfs <- purrr::map2(sample_name_vectors, papa_result_files, extract_columns_from_tsv)


tpm_dfs$brown_i3_cortical


# Calculate summarised TPM values for each le_id from a wide-format TPM matrix
calculate_mean_median <- function(df) {
  df_long <- df %>%
    pivot_longer(cols = -le_id, names_to = "sample_name", values_to = "tpm") %>%
    group_by(le_id) %>%
    summarize(
      mean_tpm = mean(tpm, na.rm = TRUE),
      median_tpm = median(tpm, na.rm = TRUE)
    )
  
  return(df_long)
}

# Calculate mean and median TPM values across KD samples for each dataset
tpm_stats <- purrr::map(tpm_dfs, calculate_mean_median)
tpm_stats$brown_i3_cortical

# Generate wide format matrices with summarised TPM per dataset

# Combine the data frames for mean TPM values
mean_tpm_combined <- tpm_stats %>%
  purrr::map(~ .x %>% select(le_id, mean_tpm)) %>%  # Extract le_id and mean_tpm
  bind_rows(.id = "dataset") %>%                   # Combine the data frames with dataset identifier
  pivot_wider(names_from = dataset, values_from = mean_tpm)  # Pivot to wide format

# Combine the data frames for median TPM values
median_tpm_combined <- tpm_stats %>%
  purrr::map(~ .x %>% select(le_id, median_tpm)) %>%  # Extract le_id and median_tpm
  bind_rows(.id = "dataset") %>%                      # Combine the data frames with dataset identifier
  pivot_wider(names_from = dataset, values_from = median_tpm)  # Pivot to wide format


if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

write_tsv(mean_tpm_combined,
          file.path(outdir, "2024-09-04_le_id.mean_tpm.kd_samples.tsv"),
          col_names = T)

write_tsv(median_tpm_combined,
          file.path(outdir, "2024-09-04_le_id.median_tpm.kd_samples.tsv"),
          col_names = T)
