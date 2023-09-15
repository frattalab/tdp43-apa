library(tidyverse)

df <- read_tsv("data/2023-05-24_i3_cortical_zanovello.all_datasets.dexseq_apa.results.processed.cleaned.tsv")
mv_df <- read_tsv("data/bleedthrough_manual_validation.tsv")

# remove some of the intermediate depletion curve datasets (i.e. keep highest KD only)
# zanovello_skndz_curve_1
# zanovello_shsy5y_curve_0075

exp_to_keep <- unique(df$experiment_name)[str_detect(unique(df$experiment_name), "_curve_",negate = T) | 
                                            unique(df$experiment_name) %in% c("zanovello_skndz_curve_1", "zanovello_shsy5y_curve_0075")]

df <- filter(df, experiment_name %in% exp_to_keep)

# remove manually validated isoforms
mv_fail_ids <- filter(mv_df, event_manual_validation != "yes") %>% pull(le_id)
df <- filter(df, !le_id %in% mv_fail_ids)

# annotate cryptic events (any dataset)
df <- mutate(df, cryptic_any = padj < 0.05 & mean_PPAU_base < 0.1 & delta_PPAU_treatment_control > 0.1)

#1. number of datasets cryptic, regulated, expressed/evaluated (and fractions)

# annotate regulated, cryptic & expressed for each event
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
            n_experiments = n())
  
le_exper_summ_counts %>%
  filter(n_cryptic > 0 & n_evaluated > 0) %>%
  ggplot(aes(x = n_cryptic, y = n_evaluated)) +
  geom_bin2d(binwidth = c(1,1)) +
  stat_bin2d(geom = "text", aes(label = after_stat(count)), size = rel(3), binwidth = c(1,1)) +
  scale_fill_gradient(low = "#fee8c8", high = "#e34a33") +
  scale_x_continuous(breaks = seq(0,10,1)) +
  scale_y_continuous(breaks = seq(0,10,1)) +
  labs(title = "Cryptics are rarely called in >1 dataset",
       x = "N datasets cryptic",
       y = "N datasets tested") +
  theme_classic(base_size = 14)


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





