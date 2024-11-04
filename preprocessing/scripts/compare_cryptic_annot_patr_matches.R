library(tidyverse)
library(ggprism)

# Calculate a group-wise empirical p-value
# Null hypothesis = no difference in values between the control group and the treatment group i.e. they come from the same distribution
# test statistic = number (fraction) of times control samples are further from the control mean (absolute) than the treatment mean
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
    
    # Extract the test statistics (stat_col) for the control group samples
    T_control <- control_subset[[stat_col]]
    
    # Calculate the mean of the control group's test statistics
    T_control_mean <- mean(T_control)
    
    # Calculate the absolute deviation for the treatment group
    T_treatment_dev <- abs(T_treatment - T_control_mean)
    
    # Calculate the absolute deviations for the control groups
    T_control_dev <- abs(T_control - T_control_mean)
    
    # Calculate the two-sided p-value
    # (comparison produces True/False (1/0) -> mean = fraction of events that are True)
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

# distance thresholds to plot
plot_distancethresh <- as.character(c(10, 25, 50, 100, 200))

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
         frac = n / n_tot) %>%
  filter(distancethresh %in% plot_distancethresh)

# calc frac support for cryptic events (& standardise ormat with annot_patr_counts_long)
cryptic_patr_counts_long <- cryptic_patr_counts %>%
  mutate(distancethresh = str_remove_all(distancethresh, "^distancethresh_"),
         distancethresh = fct_inseq(distancethresh),
         n_tot = cryptic_patr_n_ids,
         frac = n / n_tot) %>%
  filter(distancethresh %in% plot_distancethresh)


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
  scale_y_continuous(limits = c(40,75), breaks = seq(0,100,10)) +
  theme_bw(base_size = 10) +
  labs(
    # title = "PolyA-tail containing read (PATR) support" ,
    #    subtitle = "Orange = cryptic events, boxplot = annotated event samples",
       x = "Distance Threshold (nt)",
       y = "Events with PATR support (%)")

annot_vs_cryp_plot

# Calculate empirical, two-sided p-value comparing % support at each threshold for cryptics vs annotated
# (two-sided = abs difference from the mean in annotated, count number of times abs annotated diff from mean is greater)

# QC: does mean reflect the 'central tendency' of the empirical distribution at each threshold?

# Calculate mean and median for each threshold
annot_frac_summary <- annot_patr_counts_long %>%
  group_by(distancethresh) %>%
  summarize(
    mean_frac = mean(frac, na.rm = TRUE),
    median_frac = median(frac, na.rm = TRUE)
  )

# Make plot across window sizes
annot_frac_hists <- annot_patr_counts_long %>%
  ggplot(aes(x = frac)) +
  facet_wrap(~ distancethresh, scales = "free_y") +
  geom_histogram(binwidth = 0.005) +
  geom_vline(data = annot_frac_summary, aes(xintercept = mean_frac), color = "#1b9e77", linetype = "dashed",) +
  geom_vline(data = annot_frac_summary, aes(xintercept = median_frac), color = "#1f78b4", linetype = "longdash",) +
  labs(
    title = "PATR support in annotated PAS samples",
    subtitle = "Dashed lines = mean (green) and median (blue)",
    x = "Events with PATR Support (fraction)",
    y = "Frequency"
  ) + 
  theme_bw(base_size = 10)

annot_frac_hists

# Mean and median are m/l indistinguishable - calculate two-sided pvalue
cryptic_annot_p_value <- calculate_p_value(cryptic_patr_counts_long, annot_patr_counts_long, group_col = "distancethresh", stat_col = "frac") %>%
  mutate(distancethresh = fct_inseq(distancethresh),
         adj_p_value = p.adjust(p_value, method = "BH")
         )

cryptic_annot_p_value

# Make plot with p-value added
annot_vs_cryp_plot_p <- annot_vs_cryp_plot + add_pvalue(cryptic_annot_p_value, x = "distancethresh", xmin = "distancethresh", y.position = 70,
                                label = "p = {p_value}", remove.bracket = T,label.size = 2.5
                              )
annot_vs_cryp_plot_p



if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

ggsave(filename = file.path(outdir, "2024-11-04_patr_support.cryptic_vs_annotated.all_thresholds.box_plot.png"),
       plot = annot_vs_cryp_plot,
       dpi = "retina",
       units = "mm", width = 125, height = 100)

ggsave(filename = file.path(outdir, "2024-11-04_patr_support.cryptic_vs_annotated.all_thresholds.box_plot.pvals.png"),
       plot = annot_vs_cryp_plot_p,
       dpi = "retina",
       units = "mm", width = 125, height = 100)

ggsave(filename = file.path(outdir, "2024-11-04_patr_support.annotated_empirical_dbrn.all_thresholds.histogram.png"),
       plot = annot_frac_hists,
       dpi = "retina",
       units = "mm", width = 125, height = 100)

# Save as PDFs
ggsave(filename = file.path(outdir, "2024-11-04_patr_support.cryptic_vs_annotated.all_thresholds.box_plot.pdf"),
       plot = annot_vs_cryp_plot,
       dpi = "retina",
       units = "mm", width = 125, height = 100)

ggsave(filename = file.path(outdir, "2024-11-04_patr_support.cryptic_vs_annotated.all_thresholds.box_plot.pvals.pdf"),
       plot = annot_vs_cryp_plot_p,
       dpi = "retina",
       units = "mm", width = 125, height = 100)

ggsave(filename = file.path(outdir, "2024-11-04_patr_support.annotated_empirical_dbrn.all_thresholds.histogram.pdf"),
       plot = annot_frac_hists,
       dpi = "retina",
       units = "mm", width = 125, height = 100)



