library(tidyverse)
library(nullranges)
library(cobalt)
set.seed(123)

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
  }, .progress = T)
  
  # Combine all iterations into single df
  bind_rows(iteration_data, .id = "iteration")
}

# generate boxplot version of 'Love plot' from a list of bal.tab objects
# TODO: separate generating plot df from plotting
create_enhanced_love_plot <- function(bal, variable_df, comparison_df,
                                      variable_col = "variable_clean",
                                      comparison_vals = c("matched vs. focal", "pool vs. focal"),
                                      comparison_col = "comparison_clean") {
  # Process the data
  plot_data <- process_balance_data(bal)
  # return(plot_data)
  
  # subset to set comparisons wish to visualise
  plot_data <- filter(plot_data, comparison %in% comparison_vals)
  
  #
  plot_data <- inner_join(variable_df, plot_data, by = "variable") %>%
    inner_join(comparison_df, by = "comparison")
  
  # Create the plot
  ggplot(plot_data, aes(x = Diff.Un, y = !!sym(variable_col))) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    geom_boxplot(aes(color = !!sym(comparison_col)), alpha = 0.7, outlier.size = 0.5, linewidth = 0.25) +
    # facet_wrap(~ Comparison, scales = "free_x", ncol = 1) +
    # scale_color_manual(values = c("Unadjusted" = "red", "Adjusted" = "blue"),
    #                    name = "Sample") +
    labs(title = "Covariate Balance",
         x = "Standardized Mean Differences",
         y = "") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
}

# load in environment from get_matched_annot_pas.R (all matchRanges objects etc.)
outdir <- "processed/curation/cryptic_annot_comparison/"
load(file.path(outdir, "2024-09-13_get_matched_annot_pas.RData"))

# median_tpm_annot_matched_strat - list of matchedData objects (output of repeated (x1000) runs of matchRanges, balancing covariates by sampling without replacement

## HOW UNIQUE ARE EXTRACTED IDS FROM EACH ITERATION?

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

# calculate overall summary statistics of all pairwise comparisons
summ_all_uniq_pw_median_strat <- uniq_pw_median_strat %>%
  pivot_longer(cols = -iteration, names_to = "iteration2", values_to = "coefficient") %>%
  # remove self comparisons (where not evaluated & set to 0)
  filter(iteration != iteration2) %>%
  # Remove other 'symmetric' comparisons (e.g. 1 vs 2 and 2 vs 1) 
  # Generate combined ID columns by sorting two IDs alphanumerically
  mutate(iter_min = pmin(iteration, iteration2),
         iter_max = pmax(iteration, iteration2)
    ) %>%
  distinct(iter_min, iter_max, .keep_all = T) %>%
  # drop intermediate columns
  select(-iter_min, -iter_max) %>%
  summarise(
    mean_coefficient = mean(coefficient),
    sd_coefficient = sd(coefficient),
    median_coefficient = median(coefficient),
    iqr_coefficient = IQR(coefficient),
    q1_coefficient = quantile(coefficient, 0.25),
    q3_Coefficient = quantile(coefficient, 0.75)
  )


## HOW EFFECTIVE IS THE COVARIATE BALANCING FOR EACH ITERATION?

# evaluate balancing of covariates from each iteration of matchRanges (standardised mean differences)
bal_median_strat <- map(median_tpm_annot_matched_strat, evaluate_balancing_cobalt, .progress = T)
# bal_median_strat[["1"]]

# example love.plot using only a single iteration
# love.plot(bal_median_strat[["1"]], drop.distance = T)

# Generate a love.plot using values from all iterations (visualising distribution as a box plot)
# Exploring efficacy of matching and whether consistently an improvement on the initial pool of events

# add cleaned covariate column names for plotting
covariate_cols_clean <- c("Number of PAS", "TPM Brown i3 Cortical", "TPM Brown SH-SY5Y", "TPM Brown SK-N-BE(2)", "TPM Humphrey i3 Cortical", "TPM Seddighi i3 Cortical", "TPM Zanovello i3 Cortical", "TPM Zanovello SH-SY-5Y (CHX)", "TPM Zanovello SH-SY5Y (Curve)", "TPM Zanovello SK-N-BE(2)")
covar_name_df <- tibble(variable = covariate_cols, variable_clean = covariate_cols_clean)
covar_name_df

# boxplot of standardised mean differences for all iterations (matched subset & all annotated)
compar_name_df <- tibble(comparison = c("matched vs. focal", "pool vs. focal"),
                         comparison_clean = c("Matched Annotated\nvs Cryptic", "All Annotated\nvs Cryptic"))


annot_cryptic_love_boxplot <- create_enhanced_love_plot(bal_median_strat, covar_name_df, compar_name_df, comparison_vals = compar_name_df$comparison)

# Update labelling & and add annotation
annot_cryptic_love_boxplot <- annot_cryptic_love_boxplot +
  labs(colour = "Comparison",
       title = NULL) +
  scale_colour_manual(values = c( "#7570b3", "#1b9e77")) +
  scale_x_continuous(limits = c(-1.2,1.2), breaks = seq(-1,1,0.5)) +
  annotate(geom = "text", x = -1, y = "TPM Zanovello SK-N-BE(2)", label = "Higher in\nCryptic set", size = 3) +
  
  annotate(geom = "text", x = 1, y = "TPM Zanovello SK-N-BE(2)", label = "Higher in\nAnnotated set", size = 3)

annot_cryptic_love_boxplot


if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

ggsave(filename = file.path(outdir, "2024-11-14_covariate_balance.love_plot.boxplot.png"),
       plot = annot_cryptic_love_boxplot,
       dpi = "retina",
       units = "mm", width = 150, height = 100)

ggsave(filename = file.path(outdir, "2024-11-14_covariate_balance.love_plot.boxplot.pdf"),
       plot = annot_cryptic_love_boxplot,
       dpi = "retina",
       units = "mm", width = 150, height = 100)


