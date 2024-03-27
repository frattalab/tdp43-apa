library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggprism)
set.seed(123)

norm_control_by_batch <- function(df, value_col, conds = c("CTRL", "TDP43KD")) {
  
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
# update condition labels
fish_counts <- mutate(fish_counts,
                      condition = if_else(condition == "TDP", "TDP43KD", condition))

# Plots/comparisons to make
# 1. Ratio of foci counts per cell (normalised with respect to control sample in each replicate), for both probes
# 2. extra-nuclear:nuclear ratio of foci counts per cell (normalised with respect to control sample in each replicate), for green probe only

# normalise foci counts and ratios with respect to control samples in each batch
# (eacj element in list is count/ratio col - all tables contain each probe)
counts_norm <- c("mean_cell", "extranuc_nuc_ratio", "nuc_extranuc_ratio") %>%
  set_names() %>%
  map(~ norm_control_by_batch(fish_counts, .x)
      )

  
# 1 sample t.test for each probe, log transforming input values with null hypotheis that difference is not eaual to log(1) (0)

counts_norm_ttest <- counts_norm %>%
  map(~ mutate(.x, TDP43KD = log(TDP43KD)) %>%
        group_by(probe) %>%
        t_test(TDP43KD ~ 0) %>% # log(1) = 0 
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

# Plot the count ratios with respective pvalues

# Make pvalue df reading for plotting
plot_df_counts_all_pval <- counts_norm_ttest_adj %>%
  filter(metric == "mean_cell") %>%
  mutate(plot_probe = factor(if_else(probe == "proximal", "Total", "Cryptic"),
                             levels = c("Total", "Cryptic")),
         group1 = "CTRL", group2 = "TDP43KD")
  
# prep coutns for plotting
plot_df_counts_all <- counts_norm$mean_cell %>%
  pivot_longer(cols = all_of(c("CTRL", "TDP43KD")),
               names_to = "condition",
               values_to = "mean_count"
  ) %>%
  mutate(plot_probe = factor(if_else(probe == "proximal", "Total", "Cryptic"),
                             levels = c("Total", "Cryptic"))) 


plot_counts_all <- plot_df_counts_all %>%
  ggplot(aes(x = condition, y = mean_count)) +
  facet_wrap("~ plot_probe") +
  scale_y_continuous(limits = c(0,10),
                     breaks = seq(0,10, 2)) +
  # pval df doesn't have replicate aesthetic, so add to plot first before plotting the values
  ggprism::add_pvalue(plot_df_counts_all_pval, y.position = 9.5,tip.length = 0,label.size = 5) +
  geom_point(aes(x = condition, y= mean_count, shape = replicate),
             data = plot_df_counts_all,
             position = position_dodge(width = 0.4),
             size = 3) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "",
       y = "Mean foci per cell ratio",
       shape = "Replicate")

# same plot, jsut y-axis starts at 1
plot_counts_all_1 <- plot_counts_all +
  scale_y_continuous(limits = c(1,10),
                     breaks = seq(1,9, 1)) +
  labs(y = "Mean foci per cell\n(normalised to CTRL)")

plot_counts_all_1


# Repeat for ratio
# Make pvalue df reading for plotting
plot_df_ratio_prox_pval <- counts_norm_ttest_adj %>%
  filter(metric == "extranuc_nuc_ratio") %>%
  mutate(group1 = "CTRL", group2 = "TDP43KD")

# prep coutns for plotting
plot_df_ratio_prox <- counts_norm$extranuc_nuc_ratio %>%
  pivot_longer(cols = all_of(c("CTRL", "TDP43KD")),
               names_to = "condition",
               values_to = "extranuc_nuc_ratio"
  ) %>% filter(probe == "proximal")

# plot
plot_ratio_prox <- plot_df_ratio_prox %>%
  ggplot(aes(x = condition, y = extranuc_nuc_ratio)) +
  scale_y_continuous(limits = c(0,2.25),
                     breaks = seq(0,2, 0.5)) +
  # pval df doesn't have replicate aesthetic, so add to plot first before plotting the values
  ggprism::add_pvalue(plot_df_ratio_prox_pval, label = "p.adj", y.position = 2.1, tip.length = 0,label.size = 5) +
  geom_point(aes(x = condition, y= extranuc_nuc_ratio, shape = replicate),
             data = plot_df_ratio_prox,
             position = position_dodge(width = 0.3),
             size = 3) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "",
       y = "Extranuclear:nuclear foci ratio\n(normalised to CTRL)",
       shape = "Replicate")

plot_ratio_prox

# same plot, just shrink the range so y axis starts at 1
plot_ratio_prox_1 <- plot_ratio_prox + 
  scale_y_continuous(limits = c(1,2.25),
                     breaks = seq(1,2,0.5))

plot_ratio_prox_1



### Supplementaty figure updates


# 1. How variable are foci & nuclear counts replicate to replicate?

plot_nuclei_sums <- fish_counts %>%
  mutate(plot_probe = factor(if_else(probe == "proximal", "Total", "Cryptic"),
                             levels = c("Total", "Cryptic"))
         ) %>%
  ggplot(aes(x = condition, y = sum_nuclei, shape = as.factor(replicate))) +
  facet_wrap("~ plot_probe") +
  scale_y_continuous(limits = c(0,1000), breaks = seq(0,1000, 200)) +
  geom_point(size = 3) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "",
       y = "Total nuclei",
       shape = "Replicate")

plot_foci_sums <- fish_counts %>%
  mutate(plot_probe = factor(if_else(probe == "proximal", "Total", "Cryptic"),
                             levels = c("Total", "Cryptic"))
  ) %>%
  ggplot(aes(x = condition, y = sum_foci, shape = as.factor(replicate))) +
  facet_wrap("~ plot_probe") +
  scale_y_continuous(limits = c(0,1000), breaks = seq(0,1000, 200)) +
  geom_point(size = 3) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "",
       y = "Total foci",
       shape = "Replicate")
  
plot_nuclei_sums
plot_foci_sums


# Nuclear/extranuclear foci counts
plot_subcellfoci_sums <- fish_counts %>%
  mutate(plot_probe = factor(if_else(probe == "proximal", "Total", "Cryptic"),
                             levels = c("Total", "Cryptic"))
  ) %>%
  pivot_longer(cols = starts_with("sum_foci_"), names_to = "subcellular_fraction", values_to = "sum_subcell_foci", names_prefix = "sum_foci_") %>%
  mutate(plot_subcell_fraction = str_to_sentence(subcellular_fraction)) %>%
  ggplot(aes(x = condition, y = sum_subcell_foci, colour = plot_subcell_fraction, shape = as.factor(replicate))) +
  facet_wrap("~ plot_probe") +
  # scale_y_continuous(limits = c(0,1000), breaks = seq(0,1000, 200)) +
  geom_point(size = 3,
             position = position_dodge(width = 0.5)
             ) +
  scale_colour_manual(values = c("#b2df8a", "#1f78b4")) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "",
       y = "Total foci",
       shape = "Replicate",
       colour = "Location")


### Output plots


if (!dir.exists("processed/")) {dir.create("processed")}

ggsave("2024-01-09_fish_probe_count_ratio_all_cell_facet.png",
       plot = plot_counts_all,
       path = "processed/",
       device = "png",
       units = "in",
       height = 5,
       width = 5,
       dpi = "retina")

ggsave("2024-01-09_fish_probe_count_ratio_all_cell_facet.svg",
       plot = plot_counts_all,
       path = "processed/",
       device = svg,
       units = "in",
       height = 5,
       width = 5,
       dpi = "retina")

ggsave("2024-01-09_fish_prox_subcell_ratio.png",
       plot = plot_ratio_prox,
       path = "processed/",
       device = "png",
       units = "in",
       height = 5,
       width = 5,
       dpi = "retina")

ggsave("2024-01-09_fish_prox_subcell_ratio.svg",
       plot = plot_ratio_prox,
       path = "processed/",
       device = svg,
       units = "in",
       height = 5,
       width = 5,
       dpi = "retina")

ggsave("2024-01-09_fish_probe_count_ratio_all_cell_facet_start1.png",
       plot = plot_counts_all_1,
       path = "processed/",
       device = "png",
       units = "in",
       height = 5,
       width = 5,
       dpi = "retina")

ggsave("2024-01-09_fish_probe_count_ratio_all_cell_facet_start1.svg",
       plot = plot_counts_all_1,
       path = "processed/",
       device = svg,
       units = "in",
       height = 5,
       width = 5,
       dpi = "retina")


ggsave("2024-01-09_fish_prox_subcell_ratio_start1.png",
       plot = plot_ratio_prox_1,
       path = "processed/",
       device = "png",
       units = "in",
       height = 5,
       width = 5,
       dpi = "retina")

ggsave("2024-01-09_fish_prox_subcell_ratio_start1.svg",
       plot = plot_ratio_prox_1,
       path = "processed/",
       device = svg,
       units = "in",
       height = 5,
       width = 5,
       dpi = "retina")
