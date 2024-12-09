library(tidyverse)
library(nullranges)
library(cobalt)
requireNamespace("ks")
set.seed(123)

# Function to add pseudocount and log2 transform using tidyverse
add_pseudocount_and_log2 <- function(data, columns, pseudocount = 1.01) {
  data %>%
    mutate(across(all_of(columns), ~log2(.x + pseudocount)))
}


# Function to perform repeated matchRanges using replicate
repeat_matchRanges <- function(focal_df, pool_df, covariate_form, n_iterations = 1, method = "stratified", replace = FALSE) {
  
  # Coerce the input to data.frames
  focal_df_converted <- as.data.frame(focal_df)
  pool_df_converted <- as.data.frame(pool_df)
  
  replicate(n_iterations, {
    matchRanges(
      focal = focal_df_converted,
      pool = pool_df_converted,
      covar = covariate_form,
      method = method,
      replace = replace
    )
  }, simplify = FALSE)
  
}


# Function to extract IDs from each matchRanges object and write to a text file
write_id_files <- function(obj_list, dir_path, prefix) {
  
  # Ensure the directory exists
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Extract the names of the objects in the list
  obj_names <- names(obj_list)
  
  # Iterate through the list and write files
  walk2(
    obj_list, obj_names,
    ~ {
      # Extract the vector of IDs
      ids <- matched(.x)[, "le_id"]
      
      # Construct the file name
      file_name <- file.path(dir_path, paste0(prefix, .y, ".txt"))
      
      # Write the vector of IDs to the file
      write_lines(ids, file_name)
    }
  )
}


outdir <- "processed/curation/cryptic_annot_comparison/"
mean_tpm_combined <- read_tsv(file.path(outdir, "2024-09-04_le_id.mean_tpm.kd_samples.tsv"))
median_tpm_combined <- read_tsv(file.path(outdir, "2024-09-04_le_id.median_tpm.kd_samples.tsv"))
pas_counts <- read_tsv(file.path(outdir, "2024-09-03_le_id_pas_counts.tsv"))

# Left join mean_tpm_combined with pas_counts
mean_tpm_joined <- left_join(pas_counts, mean_tpm_combined, by = "le_id")

# Left join median_tpm_combined with pas_counts
median_tpm_joined <- left_join(pas_counts, median_tpm_combined, by = "le_id")

# Filter for annotated, non-cryptic events
mean_tpm_annot <- mean_tpm_joined %>%
  filter(annot_status == "annotated" & cryptic_status == 0)
  
median_tpm_annot <- median_tpm_joined %>%
  filter(annot_status == "annotated" & cryptic_status == 0)

# filter for cryptics
mean_tpm_cryp <- mean_tpm_joined %>%
  filter(cryptic_status == 1)

median_tpm_cryp <- median_tpm_joined %>%
filter(cryptic_status == 1)

# First attempt at performing covariate matching and sampling
# Want to do without replacement, such that each annotated event is unique (not double counting overlaps)
# rejection sampling is reportedly fastest without replacement method, and ~ 35x more annotated than cryptic events, so unlikely to be unstable due to limited difference in size

# covariates = number of unique PAS (separated by at least 12nt in either dirn) + all datasets
paste(colnames(median_tpm_cryp), collapse = " ")
# "le_id annot_status cryptic_status simple_event_type count count_clpsd12 brown_i3_cortical brown_shsy5y brown_skndz humphrey_i3_cortical seddighi_i3_cortical zanovello_i3_cortical_upf1_tdp_tdpkd_upf1ctl_vs_tdpctl_upf1ctl zanovello_shsy5y_chx_kd_ctl_vs_ctl_ctl zanovello_shsy5y_curve_0075 zanovello_skndz_curve_1"

# Basically everything from count_clpsd12 onwards - create formulat
covariate_cols <- colnames(median_tpm_cryp)[6:ncol(median_tpm_cryp)]
dataset_cols <- covariate_cols[2:length(covariate_cols)]
covariate_form <- reformulate(covariate_cols)
covariate_form



# Log transform the TPM columns
median_tpm_cryp <- add_pseudocount_and_log2(median_tpm_cryp, dataset_cols)
median_tpm_annot <- add_pseudocount_and_log2(median_tpm_annot, dataset_cols)


# attempt both methods for without replacement sampling
# Error in rejectSample(nfps, npps) : 
# kernel density estimates by ks::kde are negative, cannot perform rejection sampling
# Suspect due to very small log2 TPM + 1 values when untransformed TPM = 0 or very close (i.e. a boundary effect)

# median_tpm_annot_matched_rej <- repeat_matchRanges(focal_df = median_tpm_cryp,
#                    pool_df = median_tpm_annot,
#                    covariate_form = covariate_form,
#                    n_iterations = 5,
#                    method = "rejection",
#                    replace = F)

# sample annotated IDs without replacement, attempting to balance for expression across datasets and number of unique PAS
median_tpm_annot_matched_strat <- repeat_matchRanges(focal_df = median_tpm_cryp,
                                                         pool_df = median_tpm_annot,
                                                         covariate_form = covariate_form,
                                                         n_iterations = 1000)

# output matched annotated IDs to file, one per iteration
median_tpm_annot_matched_strat <- set_names(median_tpm_annot_matched_strat, as.character(seq(1:length(median_tpm_annot_matched_strat))))
outdir_ids <- file.path(outdir, "ids")
write_id_files(median_tpm_annot_matched_strat, outdir_ids, prefix = "iteration_")

# save global env to disk for further analysis
save.image(file.path(outdir, "2024-09-13_get_matched_annot_pas.RData"))





