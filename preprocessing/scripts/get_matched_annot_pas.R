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

# Function to calculate Szymkiewicz–Simpson coefficient ('overlap coefficient')
calculate_szymkiewicz_simpson <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  min_length <- min(length(set1), length(set2))
  return(intersection / min_length)
}

# Evaluate uniqueness of IDs in a series of matchRanges objects in pairwise fashion (Szymkiewicz–Simpson coefficient or 'overlap' coefficient)
evaluate_uniqueness_pairwise <- function(matchRanges_output) {
  
  n_iterations <- length(matchRanges_output)
  
  # Extract identifiers for each iteration
  identifiers <- lapply(matchRanges_output, function(x) matched(x)[, "le_id"])
  
  # Initialize results matrix
  results <- matrix(0, nrow = n_iterations, ncol = n_iterations)
  
  # Calculate Szymkiewicz–Simpson coefficient for each pair of iterations
  for (i in 1:n_iterations) {
    for (j in 1:n_iterations) {
      if (i != j) {
        results[i, j] <- calculate_szymkiewicz_simpson(identifiers[[i]], identifiers[[j]])
      }
    }
  }
  
  # Convert results to a tibble for easier handling
  results_tbl <- as_tibble(results)
  colnames(results_tbl) <- paste0("iteration_", 1:n_iterations)
  results_tbl <- results_tbl %>% 
    mutate(iteration = paste0("iteration_", 1:n_iterations)) %>%
    relocate(iteration)
  
  return(results_tbl)
}

# Function to evaluate uniqueness using Szymkiewicz–Simpson coefficient against union of other iterations
evaluate_uniqueness_union <- function(matchRanges_output) {
  n_iterations <- length(matchRanges_output)
  
  # Extract identifiers for each iteration
  identifiers <- lapply(matchRanges_output, function(x) matched(x)[, "le_id"])
  
  # Initialize results vector
  results <- numeric(n_iterations)
  
  # Calculate Szymkiewicz–Simpson coefficient for each iteration against union of others
  for (i in 1:n_iterations) {
    other_iterations <- setdiff(1:n_iterations, i)
    union_others <- Reduce(union, identifiers[other_iterations])
    results[i] <- calculate_szymkiewicz_simpson(identifiers[[i]], union_others)
  }
  
  # Convert results to a tibble for easier handling
  results_tbl <- tibble(
    iteration = paste0("iteration_", 1:n_iterations),
    coefficient = results
  )
  
  return(results_tbl)
}

# calculate summary statistics (median, iqr) from the output matrix of evaluate_uniqueness_pairwise
summarise_pairwise_uniqueness <- function(uniqueness_results) {
  uniqueness_results %>%
    pivot_longer(cols = -iteration, names_to = "compared_to", values_to = "coefficient") %>%
    # remove the self comparison 
    filter(iteration != compared_to) %>%
    group_by(iteration) %>%
    summarise(
      mean_coefficient = mean(coefficient),
      sd_coefficient = sd(coefficient),
      median_coefficient = median(coefficient),
      iqr_coefficient = IQR(coefficient),
      q1_coefficient = quantile(coefficient, 0.25),
      q3_Coefficient = quantile(coefficient, 0.75)
    )
}

# evaluate balancing from matchRanges object using cobalt
evaluate_balancing_cobalt <- function(mgr) {
  bal.tab(
    f.build("set", covariates(mgr)),
    data = matchedData(mgr),
    distance = "ps",       # name of column containing propensity score
    focal = "focal",       # name of focal group in set column
    which.treat = "focal", # compare everything to focal
    s.d.denom = "all"      # how to adjust standard deviation
  )
}

# Function to process the input list of bal.tab objects (generated from matchRanges objects) and create a dataframe
process_balance_data <- function(bal) {
  # Outer map loop over elements in bal
  iteration_data <- imap(bal, function(bal_element, iteration) {
    comparison_names <- names(bal_element$Pair.Balance) %>% set_names()
    
    # Inner map loop to extract std mean diffs for each 'comparison'
    comparison_data <- map(comparison_names, ~ bal_element$Pair.Balance[[.x]]$Balance %>%
                             as_tibble(rownames = "variable") %>%
                             mutate(comparison = .x)
    )
    
    # Combine comparison data for this iteration into single df
    bind_rows(comparison_data, .id = "comparison")
  })
  
  # Combine all iterations into single df
  bind_rows(iteration_data, .id = "iteration")
}

# generate boxplot version of 'Love plot' from a list of bal.tab objects
create_enhanced_love_plot <- function(bal, variable_df, variable_col = "variable_clean",
                                      comparison_vals = c("matched vs. focal", "pool vs. focal")) {
  # Process the data
  plot_data <- process_balance_data(bal)
  # return(plot_data)
  
  # subset to set comparisons wish to visualise
  plot_data <- filter(plot_data, comparison %in% comparison_vals)
  
  #
  plot_data <- inner_join(variable_df, plot_data, by = "variable")
  
  # Create the plot
  ggplot(plot_data, aes(x = Diff.Un, y = !!sym(variable_col))) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    geom_boxplot(aes(color = comparison), alpha = 0.7) +
    # facet_wrap(~ Comparison, scales = "free_x", ncol = 1) +
    # scale_color_manual(values = c("Unadjusted" = "red", "Adjusted" = "blue"),
    #                    name = "Sample") +
    labs(title = "Covariate Balance",
         x = "Standardized Mean Differences",
         y = "") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
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
# median_tpm_annot_matched_rej <- repeat_matchRanges(focal_df = median_tpm_cryp,
#                    pool_df = median_tpm_annot,
#                    covariate_form = covariate_form,
#                    n_iterations = 5,
#                    method = "rejection",
#                    replace = F)


median_tpm_annot_matched_strat <- repeat_matchRanges(focal_df = median_tpm_cryp,
                                                         pool_df = median_tpm_annot,
                                                         covariate_form = covariate_form,
                                                         n_iterations = 100)


# evaluate per-iteration uniqueness of sampled ids in pairwise manner
# overlap_index/Szymkiewicz–Simpson coefficient ranges from 0-1, where 0 = completely unique, 1 = complete overlap
uniq_pw_median_strat <- evaluate_uniqueness_pairwise(median_tpm_annot_matched_strat)
uniq_pw_median_strat

# repeat, but evaluate uniqueness relative to all other iterations
uniq_un_median_strat <- evaluate_uniqueness_union(median_tpm_annot_matched_strat)
uniq_un_median_strat

# get median + iqr of pairwise comparisons for each sample
summ_uniq_pw_median_strat <- summarise_pairwise_uniqueness(uniq_pw_median_strat)
summ_uniq_pw_median_strat

# evaluate balancing of covariates from each iteration of matchRanges (standardised mean differences)
bal_median_strat <- map(median_tpm_annot_matched_strat, evaluate_balancing_cobalt)
bal_median_strat[[1]]

# example love.plot using only a single iteration
love.plot(bal_median_strat[[1]], drop.distance = T)

# Generate a love.plot using values from all iterations (visualising distribution as a box plot)
# Exploring efficacy of matching and whether consistently an improvement on the initial pool of events

# process_balance_data(head(bal_median_strat)) %>% View()

# add cleaned covariate column names for plotting
covariate_cols_clean <- c("Number of PAS", "Brown i3 Cortical", "Brown SH-SY5Y", "Brown SKNDZ", "Humphrey i3 Cortical", "Seddighi i3 Cortical", "Zanovello i3 Cortical", "Zanovello SH-SY-5Y (CHX)", "Zanovello SH-SY5Y (Curve)", "Zanovello SK-N-BE(2)")
covar_name_df <- tibble(variable = covariate_cols, variable_clean = covariate_cols_clean)
covar_name_df

# boxplot of standardised mean differences for all iterations (matched subset & all annotated)
create_enhanced_love_plot(bal_median_strat, covar_name_df)
