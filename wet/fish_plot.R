library(tidyverse)
library(rstatix)
library(ggpubr)

norm_control_by_batch <- function(df, value_col, conds = c("CTRL", "TDP")) {
  
  # group by batch
  rep_grpd <- fish_counts %>%
    group_by(replicate) 
  
  rep_grpd %>%
    #list of dfs
    group_split() %>%
    #set names to batch values for later tracking
    set_names(pull(group_keys(rep_grpd))) %>%
    # convert count col rows to columns for each condition
    map(~ pivot_wider(.x, id_cols = probe, names_from = condition, values_from = !!sym(value_col))
    ) %>%
    # normalise each condition with respect to control
    map(~ mutate(.x, across(all_of(conds),
                            ~ .x / !!sym(conds[1])
                            )
                 )
    ) %>%
    # combine batch-wise control normalised values to single df
    bind_rows(.id = "replicate")
  
  
}


fish_counts <- read_tsv("processed/2023-10-13_fish_counts_processed.tsv")

# Plots/comparisons to make
# 1. Ratio of foci counts per cell (normalised with respect to control sample in each replicate), for both probes
# 2. extra-nuclear:nuclear ratio of foci counts per cell (normalised with respect to control sample in each replicate), for green probe only

# normalise foci counts and ratios with respect to control samples in each batch
# (eacj element in list is count/ratio col - all tables contain each probe)
counts_norm <- c("mean_cell", "extranuc_nuc_ratio", "nuc_extranuc_ratio") %>%
  set_names() %>%
  map(~ norm_control_by_batch(fish_counts, .x))
  
# 1 sample t.test for each probe, log transforming input values with null hypotheis that difference is not eaual to log(1) (0)

counts_norm_ttest <- counts_norm %>%
  map(~ mutate(.x, TDP = log(TDP)) %>%
        group_by(probe) %>%
        t_test(TDP ~ 0) %>% # log(1) = 0 
        ungroup()
      ) %>%
  bind_rows(.id = "metric")

# adjust p-values for multiple testing
# since the ratios give equivalent p-values, drop the nuc_extranuc ratios
# also drop the ratio for distal probe (as don't expect much in CTRL, and proximal probe can represent either isoform)
counts_norm_ttest_adj <- counts_norm_ttest %>%
  filter(metric != "nuc_extranuc_ratio",
         !(metric == "extranuc_nuc_ratio" & probe == "distal")
         ) %>%
  adjust_pvalue(p.col = "p", output.col = "p.adj", method = "BH") %>%
  add_significance(p.col = "p.adj", output.col = "p.adj.signif")

counts_norm_ttest_adj  
