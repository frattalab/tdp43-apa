library(tidyverse)
library(ggupset)

facet_heatmap <- function(df,
                          plot_title = "Cryptic event deltas across datasets",
                          plot_subtitle = "** = cryptic criteria, * = padj < 0.05, blank = padj > 0.05",
                          plot_x = "Dataset",
                          plot_y = "Last exon ID",
                          plot_fill = "PolyA site usage % (TDP43KD - CTRL)") {
  df %>%
  ggplot(aes(x = experiment_name_simple, y = plot_le_id, fill = delta_PPAU_treatment_control, label = plot_label)) + 
    facet_wrap("~ plot_event_type", ncol = 3, scales = "free_y") +
    geom_tile() +
    scale_fill_gradient2(low ="#998ec3", mid = "#f7f7f7", high = "#f1a340",
                         breaks = seq(-1,1,0.2)) + # "#fee8c8"
    geom_text(size = rel(3), nudge_y = 0 ) + # -0.25
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(size = rel(0.75), angle = 90),
          axis.text.y = element_text(size = rel(0.55)),
          legend.position = "top",
          legend.key.width = unit(2, "cm")
    ) +
    labs(title = plot_title,
         subtitle = plot_subtitle,
         x = plot_x,
         y = plot_y,
         fill = plot_fill)
}


df <- read_tsv("data/2023-12-10_i3_cortical_zanovello.all_datasets.dexseq_apa.results.processed.cleaned.tsv")
mv_df <- read_tsv("data/bleedthrough_manual_validation.tsv")
cryptics_df <- read_tsv("data/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation_complex.tsv")


# remove manually validated isoforms
mv_fail_ids <- filter(mv_df, event_manual_validation != "yes") %>% pull(le_id)
df <- filter(df, !le_id %in% mv_fail_ids)

# annotate cryptic events (any dataset)
df <- mutate(df, cryptic_any = le_id %in% unique(cryptics_df$le_id))

# Define a cleaned df mapping experiment names to plot names
experiment_names <- tibble(experiment_name = c("humphrey_i3_cortical",
                                          "brown_i3_cortical",
                                          "seddighi_i3_cortical",
                                          "zanovello_i3_cortical_upf1_tdp_tdpkd_upf1ctl_vs_tdpctl_upf1ctl",
                                          "klim_i3_motor",
                                          "zanovello_shsy5y_curve_0075",
                                          "zanovello_shsy5y_chx_kd_ctl_vs_ctl_ctl",
                                          "brown_shsy5y",
                                          "zanovello_skndz_curve_1",
                                          "brown_skndz"),
                      plot_experiment_name = c("Humphrey i3 cortical",
                                               "Brown i3 cortical",
                                               "Seddighi i3 cortical",
                                               "Zanovello i3 cortical",
                                               "Klim i3 motor",
                                               "Zanovello SH-SY-5Y curve",
                                               "Zanovello SH-SY-5Y CHX",
                                               "Brown SH-SY-5Y",
                                               "Zanovello SK-N-BE(2) curve",
                                               "Brown SK-N-BE(2)")
                      ) %>%
  mutate(plot_experiment_name = factor(plot_experiment_name, levels = c("Humphrey i3 cortical",
                                                                        "Brown i3 cortical",
                                                                        "Seddighi i3 cortical",
                                                                        "Zanovello i3 cortical",
                                                                        "Klim i3 motor",
                                                                        "Zanovello SH-SY-5Y curve",
                                                                        "Zanovello SH-SY-5Y CHX",
                                                                        "Brown SH-SY-5Y",
                                                                        "Zanovello SK-N-BE(2) curve",
                                                                        "Brown SK-N-BE(2)"))
         )




#1. number of datasets cryptic, regulated, expressed/evaluated (and fractions)

# annotate regulated, cryptic & expressed for each event + dataset
df <- df %>%
  mutate(cryptic = padj < 0.05 & mean_PPAU_base < 0.1 & delta_PPAU_treatment_control > 0.1,
         regulated = !cryptic & padj < 0.05,
         evaluated = T)

# get all possible combos of le_id and experiment name, fill in missing rows with False for 'evaluated' column
df_all_combos <- expand(df, le_id, experiment_name) %>%
  left_join(select(df, le_id, experiment_name, gene_name, pvalue, padj, mean_PPAU_base, mean_PPAU_treatment, delta_PPAU_treatment_control, cryptic, regulated, evaluated),
            by = c("le_id", "experiment_name")) %>%
  replace_na(list("cryptic" = F, "regulated" = F, evaluated = F)) 

# now generate a count of datasets where event cryptic, regulated, evaluated
# TODO: add gene_name to this df so can output
le_exper_summ_counts <- df_all_combos %>%
  distinct(le_id, experiment_name, cryptic, regulated, evaluated) %>%
  group_by(le_id) %>%
  summarise(n_cryptic = sum(cryptic),
            n_regulated = sum(regulated),
            n_evaluated = sum(evaluated),
            n_experiments = n()) %>%
  ungroup()

# double check number of cryptic events is the same
le_exper_summ_counts %>%
  filter(n_cryptic > 0) %>%
  nrow()
# [1] 227
  
plot_cryptic_vs_expressed_counts <- le_exper_summ_counts %>%
  filter(n_cryptic > 0 & n_evaluated > 0) %>%
  count(n_evaluated, n_cryptic) %>%
  ggplot(aes(x = n_cryptic, y = n_evaluated, label = n, fill = n)) +
  geom_tile() +
  scale_fill_gradient(low = "#fee8c8", high = "#e34a33") +
  geom_text(size = rel(6)) +
  scale_x_continuous(labels = seq(1,10,1),
                     breaks = seq(1,10,1)) +
  scale_y_continuous(labels = seq(1,10,1),
                     breaks = seq(1,10,1)) +
  labs(title = "Cryptics are rarely called in >1 dataset",
       x = "Number of datasets cryptic",
       y = "Number of datasets expressed",
       fill = "Count") +
  theme_classic(base_size = 20)

plot_cryptic_vs_expressed_counts


plot_cryptic_vs_expressed_counts_notitle <- plot_cryptic_vs_expressed_counts + 
  labs(title = "")

# Plot number of cryptic events per dataset

# get df of cryptic counts per dataset
dataset_cryptic_counts <-  df %>%
  filter(cryptic) %>%
  distinct(experiment_name, le_id) %>%
  count(experiment_name) %>%
  # add in tidy dataset name for plotting
  left_join(experiment_names, by = "experiment_name")

barplot_cryptic_counts <- dataset_cryptic_counts %>%
  ggplot(aes(x = n,
             y = plot_experiment_name)) +
  geom_col() +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(0,125)) +
  geom_text(aes(y = plot_experiment_name, label = n), hjust = -0.25, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Number of Cryptics",
       y = "Dataset")

barplot_cryptic_counts


# Upset plot of overlap of cryptic events between datasets
upsetplot_cryptic_counts <- df %>%
  filter(cryptic) %>%
  distinct(experiment_name, le_id) %>%
  left_join(experiment_names, by = "experiment_name") %>%
  group_by(le_id) %>%
  summarise(datasets = list(plot_experiment_name)) %>%
  ggplot(aes(x = datasets)) +
  geom_bar() +
  scale_x_upset() +
  scale_y_continuous(limits = c(0,65)) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.25, size = 2) +
  labs(x = "Dataset",
       y = "Count") +
  theme_bw()

upsetplot_cryptic_counts


#2. Calculate regulation consistency score (basically sum of absolute/ signed-ranks of -log10 pvalues
# Way to rank/visualise events according to their consistency score

# 1. Score last exons by -log10(pvalue) * delta usage (upweights significant pvalues, then scales by magnitude of change)
# 2a) sum scaled score across datasets, perhaps also add normalised to number of datasets
# 2b) normalise within experiment - min-max scale then multiple by sign of delta, 
# OR calculate fractional rank of score, sign by delta & sum

working_df <- df %>%
  select(experiment_name, le_id, gene_name, pvalue, padj, contains("PPAU")) %>%
  distinct(.keep_all = T)

# some cases where p-values are 0 - need to replace these so don't get infinite -log10(pvalue)
working_df %>%
  filter(pvalue == 0)

# replace experiment_wise with the minimum pvalue from the experiment (shifted down by 1 % to make smaller pvalues)
exper_min_pval <- working_df %>%
  filter(pvalue != 0) %>%
  group_by(experiment_name) %>%
  summarise(min_pvalue = min(pvalue))

working_df <- working_df %>%
  left_join(exper_min_pval, by = "experiment_name") %>%
  # if not 0 keep existing pvalue, otherwise make 1* smaller than minimium non-zero pvalue
  mutate(score_pvalue = if_else(pvalue == 0, 
                                min_pvalue - (min_pvalue * 0.01),
                                pvalue))

# calculate regulation score
working_df <- working_df %>%
  mutate(regn_score = -log10(score_pvalue) * delta_PPAU_treatment_control)


# sum regulation score across experiments
regn_score_sum <- working_df %>%
  group_by(le_id, gene_name) %>%
  summarise(sum_regn_score = sum(regn_score)) %>%
  arrange(desc(sum_regn_score))
  
# add number of datasets summary counts & normalise to number of datasets evaluated
regn_score_sum <- regn_score_sum %>%
  left_join(le_exper_summ_counts, by = c("le_id")) %>%
  mutate(sum_regn_score_norm = sum_regn_score / n_evaluated) %>%
  arrange(desc(sum_regn_score_norm))

# TODO: to get 'crypticy-ness', downweight each regulation score based on usage in control cells/'closeness' to threshold?
# e.g. further weight score by 1 - control usage (so bigger usage in control cells = smaller value)
# also will penalise further events that are cryptic in one dataset but physiologically expressed in another dataset
## (if that doesn't appear in ranks, should be able to see by difference in unadjusted and adjusted score)

# Make cryptic heatmap plot, ordering within event type by regulation score, fill scaled by change in usage
# Group datasets by event type

# get to plot df
# generate readable le_id - gene_name + _<1> isoform suffix
# add label of significance threshold
# text label for whether classed as cryptic, regulated or not

# define plotting order of datasets
unique(df$experiment_name)
plot_exper_name_order <- c("humphrey_i3_cortical",
                           "brown_i3_cortical",
                           "seddighi_i3_cortical",
                           "zanovello_i3_cortical_upf1_tdp_tdpkd_upf1ctl_vs_tdpctl_upf1ctl",
                           "klim_i3_motor",
                           "zanovello_shsy5y_curve_0075",
                           "zanovello_shsy5y_chx_kd_ctl_vs_ctl_ctl",
                           "brown_shsy5y",
                           "zanovello_skndz_curve_1",
                           "brown_skndz"
)

plot_exper_name_simple <- c("Humphrey i3 cortical",
                            "Brown i3 cortical",
                            "Seddighi i3 cortical",
                            "Zanovello i3 cortical",
                            "Klim i3 motor",
                            "Zanovello SH-SY-5Y curve",
                            "Zanovello SH-SY-5Y CHX",
                            "Brown SH-SY-5Y",
                            "Zanovello SK-N-BE(2) curve",
                            "Brown SK-N-BE(2)")

# make a named vector where names are the existing values
names(plot_exper_name_simple) <- plot_exper_name_order
names(plot_exper_name_simple)


# order first by per-experiment average regulation score
plot_df <- df %>%
  mutate(plot_le_id = paste(gene_name, str_split_i(le_id, "_", 2), sep = "_"), 
       # plot_shape = if_else(padj < 0.05, T, F),
       plot_label = case_when(cryptic ~ "**",
                              regulated ~ "*",
                              TRUE ~ "")
       ) %>%
  # add in scores
  left_join(select(regn_score_sum, -gene_name),
            by = c("le_id")) %>%
  # Arrange le_ids in order descending order of regulation_score
  # Arange datasets in chosen order order
  mutate(plot_le_id = fct_reorder(plot_le_id, sum_regn_score_norm),
         experiment_name_simple = plot_exper_name_simple[experiment_name],
         experiment_name_simple = factor(experiment_name_simple, levels = plot_exper_name_simple),
         plot_event_type = case_when(simple_event_type == "bleedthrough" ~ "IPA",
                                     simple_event_type == "distal_3utr_extension" ~ "3'Ext",
                                     simple_event_type == "spliced" ~ "ALE",
                                     T ~ ""),
         plot_event_type = factor(plot_event_type, levels = c("ALE", "IPA", "3'Ext"))
         )


# normalised score plot df
cryptic_norm_score_plot_df <- plot_df %>%
  # first subset for events that are cryptic in at least one dataset
  group_by(le_id) %>%
  filter(any(cryptic)) %>%
  ungroup()


# Remove redundant/unnecessary suffixes from name labels (i.e. keep only when multiple isoforms for same gene)
cryptic_norm_score_plot_df <- cryptic_norm_score_plot_df %>%
  group_by(gene_name) %>%
  mutate(gene_ncryp = n_distinct(le_id),
         plot_le_id = if_else(gene_ncryp == 1, gene_name, plot_le_id)) %>%
  ungroup() %>%
  select(-gene_ncryp) %>%
  # Re-arrange le_ids in order descending order of regulation_score
  mutate(plot_le_id = fct_reorder(plot_le_id, sum_regn_score_norm))


cryptics_norm_score_sort_heatmap <- cryptic_norm_score_plot_df %>%
  facet_heatmap(plot_title = "Cryptic event deltas across datasets - sorted by normalised sum of regulation score")


# repeat - this time just using sum of regulation score across datasets
cryptic_sum_score_plot_df <- plot_df %>%
  # first subset for events that are cryptic in at least one dataset
  group_by(le_id) %>%
  filter(any(cryptic)) %>%
  ungroup()

# Remove redundant/unnecessary suffixes from name labels
cryptic_sum_score_plot_df <- cryptic_sum_score_plot_df %>%
  group_by(gene_name) %>%
  mutate(gene_ncryp = n_distinct(le_id),
         plot_le_id = if_else(gene_ncryp == 1, gene_name, plot_le_id)) %>%
  ungroup() %>%
  select(-gene_ncryp) %>%
  # Re-order using sum of regulation score across datasets
  mutate(plot_le_id = fct_reorder(plot_le_id, sum_regn_score))

cryptics_sum_score_sort_heatmap <- cryptic_sum_score_plot_df %>%
  facet_heatmap(plot_title = "Cryptic event deltas across datasets - sorted by sum of regulation score")

# make both plots without axis labels
cryptics_norm_score_sort_heatmap_nolab <- facet_heatmap(cryptic_norm_score_plot_df,
                                                       plot_title = "",
                                                       plot_subtitle = "",
                                                       plot_x = "",
                                                       plot_y = "",
)


cryptics_sum_score_sort_heatmap_nolab <- facet_heatmap(cryptic_sum_score_plot_df,
                                                       plot_title = "",
                                                       plot_subtitle = "",
                                                       plot_x = "",
                                                       plot_y = "",
                                                       )

cryptics_norm_score_sort_heatmap
cryptics_norm_score_sort_heatmap_nolab
cryptics_sum_score_sort_heatmap
cryptics_sum_score_sort_heatmap_nolab

if (!dir.exists("processed")) {dir.create("processed", recursive = T)}


ggsave(filename = "2024-11-15_ndatasets_cryptic_vs_expressed_binplot.png",
       plot = plot_cryptic_vs_expressed_counts,
              path = "processed",
              width = 10,
              height = 10,
              units = "in",
              dpi = "retina")

ggsave(filename = "2024-11-15_ndatasets_cryptic_vs_expressed_binplot.svg",
       plot = plot_cryptic_vs_expressed_counts,
       device = svg,
       path = "processed",
       width = 10,
       height = 10,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-15_ndatasets_cryptic_vs_expressed_binplot_notitle.png",
       plot = plot_cryptic_vs_expressed_counts_notitle,
       path = "processed",
       width = 10,
       height = 10,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-15_ndatasets_cryptic_vs_expressed_binplot_notitle.svg",
       plot = plot_cryptic_vs_expressed_counts_notitle,
       device = svg,
       path = "processed",
       width = 10,
       height = 10,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-15_cryptics_normed_regn_score_sort_facet_heatmap.png",
       plot = cryptics_norm_score_sort_heatmap,
       path = "processed",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-15_cryptics_sum_regn_score_sort_facet_heatmap.png",
       plot = cryptics_sum_score_sort_heatmap,
       path = "processed",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-15_cryptics_normed_regn_score_sort_facet_heatmap.svg",
       plot = cryptics_norm_score_sort_heatmap,
       path = "processed",
       device = svg,
       width = 12,
       height = 18,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-15_cryptics_sum_regn_score_sort_facet_heatmap.svg",
       plot = cryptics_sum_score_sort_heatmap,
       path = "processed",
       device = svg,
       width = 12,
       height = 18,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-15_cryptics_sum_regn_score_sort_facet_heatmap_nolabs.svg",
       plot = cryptics_sum_score_sort_heatmap_nolab,
       path = "processed",
       device = svg,
       width = 12,
       height = 18,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-15_cryptics_sum_regn_score_sort_facet_heatmap_nolabs.pdf",
       plot = cryptics_sum_score_sort_heatmap_nolab,
       path = "processed",
       device = "pdf",
       width = 12,
       height = 18,
       units = "in",
       dpi = "retina")

# cryptic counts and overlap plots

ggsave(filename = "2024-11-14_cryptics_count_barplot.png",
       plot = barplot_cryptic_counts,
       path = "processed",
       width = 100,
       height = 100,
       units = "mm",
       dpi = "retina")

ggsave(filename = "2024-11-14_cryptics_count_barplot.pdf",
       plot = barplot_cryptic_counts,
       path = "processed",
       width = 100,
       height = 100,
       units = "mm",
       dpi = "retina")


ggsave(filename = "2024-11-14_cryptics_upsetplot.png",
       plot = upsetplot_cryptic_counts,
       path = "processed",
       width = 200,
       height = 100,
       units = "mm",
       dpi = "retina")

ggsave(filename = "2024-11-14_cryptics_upsetplot.pdf",
       plot = upsetplot_cryptic_counts,
       path = "processed",
       width = 200,
       height = 100,
       units = "mm",
       dpi = "retina")

ggsave(filename = "2024-11-14_cryptics_upsetplot.svg",
       plot = upsetplot_cryptic_counts,
       path = "processed",
       device = svg,
       width = 200,
       height = 100,
       units = "mm",
       dpi = "retina")


# write enrichment scores to tsv (summarised)
write_tsv(regn_score_sum, "processed/2024-01-10_last_exon_regulation_scores.tsv", col_names = T)

# write df with cryptic, detected and regulated counts across datasets for each event
le_exper_summ_counts <- le_exper_summ_counts %>%
  left_join(distinct(df, le_id, gene_name)) %>%
  relocate(gene_name, .after = le_id) %>%
  arrange(desc(n_cryptic), desc(n_regulated), desc(n_evaluated))

write_tsv(le_exper_summ_counts, "processed/2024-01-10_last_exon_dataset_regulation_summary_counts.tsv", col_names = T)

# save to Rdata
regn_score_all <- working_df

save(plot_cryptic_vs_expressed_counts, cryptics_norm_score_sort_heatmap, cryptics_sum_score_sort_heatmap, regn_score_all, regn_score_sum, plot_df, cryptic_norm_score_plot_df, cryptic_sum_score_plot_df, file = "processed/supplementary_fig1_cryptics_consistency.Rdata")