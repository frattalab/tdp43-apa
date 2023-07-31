library(tidyverse)
library(broom)
library(readxl)
library(rstatix)


#
n_df <- read_csv("data/number_foci_ratio.csv")
n_df

#
ratio_df <- read_csv("data/nuclear_cyto_ratio.csv")
ratio_df

n_green_ttest <- t.test(n_df$Green, mu = 1) %>%
  broom::tidy()

n_red_ttest <- t.test(n_df$Red, mu = 1) %>%
  broom::tidy()

ratio_green_ttest <- t.test(ratio_df$Green, mu = 1) %>%
  broom::tidy()

# combine dfs
comb_df <- bind_rows(count_green = n_green_ttest,
          count_red = n_red_ttest,
          ratio_green = ratio_green_ttest,
          .id = "test_type") 

# adjust p-values for multiple testing
comb_df <- comb_df %>%
  adjust_pvalue(p.col = "p.value", output.col = "p.adj", method = "bonferroni") 


write_csv(comb_df, "fish_ttest_results.csv", col_names = T)


n_green_ttest_log <- t.test(log(n_df$Green), mu = log(1)) %>%
  broom::tidy()

n_red_ttest_log <- t.test(log(n_df$Red), mu = log(1)) %>%
  broom::tidy()

ratio_green_ttest_log <- t.test(log(ratio_df$Green), mu = log(1)) %>%
  broom::tidy()

# combine dfs
comb_df_log <- bind_rows(count_green = n_green_ttest_log,
                     count_red = n_red_ttest_log,
                     ratio_green = ratio_green_ttest_log,
                     .id = "test_type") 

# adjust p-values for multiple testing
comb_df_log <- comb_df_log %>%
  adjust_pvalue(p.col = "p.value", output.col = "p.adj", method = "BH") %>%
  add_significance(p.col = "p.adj", output.col = "p.adj.signif")


write_csv(comb_df_log, "fish_ttest_log_transform_results.csv", col_names = T)
