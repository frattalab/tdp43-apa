library(tidyverse)
library(ggprism)

# Define the function to calculate two-sided p-value
calculate_p_value <- function(treatment_df, control_df, group_col, stat_col) {
  
  # Get the distinct values of 'group_col' from the treatment dataframe
  unique_groups <- unique(treatment_df[[group_col]])
  
  # Initialize a result list to store p-values for each distinct value of 'group_col'
  results <- list()
  
  # Loop over each distinct 'group_col' value
  for (group_val in unique_groups) {
    
    # Subset the treatment and control dataframes for the current group_col value
    treatment_subset <- treatment_df[treatment_df[[group_col]] == group_val, ]
    control_subset <- control_df[control_df[[group_col]] == group_val, ]
    
    # Extract the test statistic (stat_col) for the treatment group
    T_treatment <- treatment_subset[[stat_col]]
    
    # Extract the test statistics (stat_col) for the control groups
    T_control <- control_subset[[stat_col]]
    
    # Calculate the mean of the control group's test statistics
    T_control_mean <- mean(T_control)
    
    # Calculate the absolute deviation for the treatment group
    T_treatment_dev <- abs(T_treatment - T_control_mean)
    
    # Calculate the absolute deviations for the control groups
    T_control_dev <- abs(T_control - T_control_mean)
    
    # Calculate the two-sided p-value
    p_value <- mean(T_control_dev >= T_treatment_dev)
    
    # Store the p-value for the current group_col value
    # p_values[[as.character(group_val)]] <- p_value
    # Store the results for the current group_col value
    results[[as.character(group_val)]] <- data.frame(
      group_value = group_val,
      T_treatment = T_treatment,
      T_control_mean = T_control_mean,
      T_treatment_dev_abs = T_treatment_dev,
      p_value = p_value
    )
  }    
    
  
  # # Convert the list to a data frame for easier interpretation
  # result_df <- data.frame(
  #   group = names(p_values),  # Dynamically assign the group column name
  #   p_value = unlist(p_values)
  # )
  # 

  # Combine the results into a single data frame
  result_df <- do.call(rbind, results)

  # # Rename the 'group' column to the passed group_col name
  names(result_df)[1] <- group_col

  
  return(result_df)
  
}

# TSV of number of events with support across range of distance thresholds
cryptic_patr_counts <- read_tsv("processed/curation/cryptic_annot_comparison/cryptics.patr_all.nearest_threshold.counts.tsv")
annot_patr_counts <- read_tsv("processed/curation/cryptic_annot_comparison/annotated_iterations.patr_all.nearest_threshold.counts.tsv")

# Get counts for events analysed (cryptic/each iteration of annotated)
# Using matrices of le_ids and their support (1/0) at each distance threshold
cryptic_patr_n_ids <- read_tsv("processed/curation/cryptic_annot_comparison/cryptics.patr_all.nearest_threshold.matrix.tsv") %>%
  nrow()

annot_patr_n_ids <- read_tsv("processed/curation/cryptic_annot_comparison/annotated_iterations.patr_all.nearest_threshold.matrix.tsv.gz",
                             col_select = iteration) %>%
  count(iteration)


outdir <- "processed/curation/cryptic_annot_comparison/"


# pivot annotated to a longer format (with one row per distance threshold and iteration combo) 
annot_patr_counts_long <- annot_patr_counts %>%
  pivot_longer(
    cols = starts_with("distancethresh_"),  # Select columns starting with 'distancethresh_'
    names_to = "distancethresh",             # Name the key column that contains the column names
    names_prefix = "distancethresh_",       # Remove the 'distancethresh_' prefix from column names
    values_to = "n"                         # Rename the values column to 'n'
  ) %>%
  # calculate total events with support at each threshold for each iteration
  left_join(annot_patr_n_ids, by = "iteration", suffix = c("", "_tot")) %>%
  mutate(distancethresh = fct_inseq(distancethresh),
         frac = n / n_tot)

# calc frac support for cryptic events (& standardise ormat with annot_patr_counts_long)
cryptic_patr_counts_long <- cryptic_patr_counts %>%
  mutate(distancethresh = str_remove_all(distancethresh, "^distancethresh_"),
         distancethresh = fct_inseq(distancethresh),
         n_tot = cryptic_patr_n_ids,
         frac = n / n_tot)

annot_patr_counts_long
cryptic_patr_counts_long


# plot the distribution of annotated % support vs cryptic support at each threshold
annot_vs_cryp_plot <- ggplot(data = annot_patr_counts_long, aes(x = distancethresh, y = frac * 100)) +
  # Plot distribution as boxplot or violin (choose one based on your preference)
  geom_boxplot(aes(group = distancethresh), fill = "lightgrey", alpha = 0.5, outlier.shape = NA) +
  # Optionally, you can use geom_violin instead of boxplot if desired
  # geom_violin(aes(group = distancethresh), fill = "lightblue", alpha = 0.5) +
  # Overlay the cryptic_patr_counts_long data with a slight horizontal adjustment using geom_jitter
  geom_point(data = cryptic_patr_counts_long, 
             aes(x = distancethresh, y = frac * 100), 
             size = 2, 
             position = position_nudge(x = 0.1),  # Nudge horizontally by 0.1
             # shape = 21,   # Shape of the points
             colour = "#d95f02") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  theme_bw(base_size = 10) +
  labs(title = "PolyA-tail containing read (PATR) support" ,
       subtitle = "Orange = cryptic events, boxplot = annotated event samples",
       x = "Distance Threshold (nt)",
       y = "Events with PATR support (%)")

annot_vs_cryp_plot

# Calculate empirical, two-sided p-value comparing % support at each threshold for cryptics vs annotated
# (two-sided = abs difference from the mean in annotated, count number of times abs annotated diff from mean is greater)
cryptic_annot_p_value <- calculate_p_value(cryptic_patr_counts_long, annot_patr_counts_long, group_col = "distancethresh", stat_col = "frac") %>%
  mutate(distancethresh = fct_inseq(distancethresh))

# Make plot with p-value added
annot_vs_cryp_plot_p <- annot_vs_cryp_plot + add_pvalue(cryptic_annot_p_value, x = "distancethresh", xmin = "distancethresh", y.position = 80,
                                label = "p = {p_value}", remove.bracket = T,label.size = 2.5
                              )
annot_vs_cryp_plot_p



if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

ggsave(filename = file.path(outdir, "2024-09-16_patr_support.cryptic_vs_annotated.all_thresholds.box_plot.png"),
       plot = annot_vs_cryp_plot,
       dpi = "retina",
       units = "mm", width = 125, height = 100)

ggsave(filename = file.path(outdir, "2024-09-16_patr_support.cryptic_vs_annotated.all_thresholds.box_plot.pvals.png"),
       plot = annot_vs_cryp_plot_p,
       dpi = "retina",
       units = "mm", width = 125, height = 100)




