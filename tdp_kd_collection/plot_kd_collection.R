library(tidyverse)
library(ggrepel)
library(ggupset)
set.seed(123)

#' Convert dataframe of median delta usages for signficiant events into df reading for scatter plot
#' Can label genes of interest with vector of gene_names (pass as 2nd argument)
#' Optionally labelling only cryptic events (label_all = F, default) or all events
get_scatter_df <- function(df, genes_to_label, label_all = F) {
  
  if (label_all) {
    
    df_out <- mutate(df, plot_name = if_else(gene_name %in% genes_to_label,
                                                 gene_name,
                                                 "")
    )
    
  } else {
    # only label cryptics of provided genes
    df_out <- mutate(df, plot_name = if_else(abs(median_delta) > 0.1 & median_ctl < 0.1 & gene_name %in% genes_to_label,
                                                 gene_name,
                                                 "")
                     )
    
  }
  
  
  df_out <- df_out %>%
    mutate(plot_alpha = case_when(plot_name != "" ~ 1,
                                  abs(median_delta) > 0.1 & median_ctl < 0.1 ~ 0.5,
                                  abs(median_delta) > 0.1 & median_ctl > 0.1 ~ 0.2,
                                  TRUE ~ 0.01),
           plot_colour = if_else(plot_alpha == 1 | plot_alpha == 0.5,
                                 "orange", "grey"))
  

  
  df_out

}



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

sig_df <- filter(df, padj < 0.05)

# annotate cryptic events (any dataset)
sig_df <- mutate(sig_df, cryptic_any = padj < 0.05 & mean_PPAU_base < 0.1 & delta_PPAU_treatment_control > 0.1)

# how many events are cryptic in at least 1 dataset?
cryp_any_df <- sig_df %>%
  arrange(le_id, desc(cryptic_any)) %>%
  # keep first row, favouring if cryptic in at least 1 dataset
  distinct(le_id, .keep_all = T) 

count(cryp_any_df, cryptic_any)
# A tibble: 2 × 2
# cryptic_any     n
# <lgl>       <int>
#   1 FALSE        5321
# 2 TRUE          227

# Calculate median base usage & change in usage for each le_id (i.e. each dataset in which a significant change)
sig_med_df <- sig_df %>% 
  group_by(le_id) %>%
  summarise(median_ctl = median(mean_PPAU_base),
            median_delta = median(delta_PPAU_treatment_control)) %>%
  ungroup() %>%
  left_join(select(sig_df, le_id, gene_name, event_type), by = "le_id") %>%
  distinct(le_id, median_ctl, median_delta, gene_name, event_type, .keep_all = T)

# with median criteria, how many isoforms are cryptic in one dataset but not in terms of median usage
cryp_any_med_df <- sig_med_df %>%
  mutate(cryptic_med = median_ctl < 0.1 & median_delta > 0.1) %>%
  distinct(le_id, cryptic_med, median_ctl, median_delta) %>%
  left_join(cryp_any_df,., by = "le_id") 

cryp_any_med_df %>%
  count(cryptic_any, cryptic_med)
# A tibble: 4 × 3
# cryptic_any cryptic_med     n
# <lgl>       <lgl>       <int>
#   1 FALSE       FALSE        5319
# 2 FALSE       TRUE            2
# 3 TRUE        FALSE          91
# 4 TRUE        TRUE          136

cryp_any_med_df %>%
  count(simple_event_type, cryptic_any, cryptic_med,)
# A tibble: 10 × 4
# simple_event_type     cryptic_any cryptic_med     n
# <chr>                 <lgl>       <lgl>       <int>
#   1 bleedthrough          FALSE       FALSE        1082
# 2 bleedthrough          TRUE        FALSE          21
# 3 bleedthrough          TRUE        TRUE           41
# 4 distal_3utr_extension FALSE       FALSE        1272
# 5 distal_3utr_extension FALSE       TRUE            1
# 6 distal_3utr_extension TRUE        FALSE          44
# 7 distal_3utr_extension TRUE        TRUE           61
# 8 spliced               FALSE       FALSE        3337
# 9 spliced               TRUE        FALSE          48
# 10 spliced               TRUE        TRUE           78

# seems to affect event types quite generally

# get a table of genes which miss cryptic criteria by median delta
cryp_any_not_med_df <- cryp_any_med_df %>%
  filter(cryptic_any & !cryptic_med) %>%
  distinct(le_id, gene_name, simple_event_type, median_ctl, median_delta) %>%
  arrange(desc(median_ctl))

# Categorise cryptics according to their expression criteria for medians of significant datasets
cryp_type_med_df <- cryp_any_med_df %>%
  filter(cryptic_any | cryptic_med) %>%
  distinct(le_id, gene_name, simple_event_type, median_ctl, median_delta, cryptic_any, cryptic_med) %>%
  mutate(median_cryp_status = case_when(cryptic_med ~ "median_cryptic",
                                        median_ctl < 0.1 & median_delta < 0 ~ "median_low_exprn_negative_delta",
                                        median_ctl > 0.1 & median_delta < 0 ~ "median_high_exprn_negative_delta",
                                        median_ctl < 0.1 & median_delta < 0.1 ~ "median_low_exprn_low_delta",
                                        median_ctl >= 0.1 & median_delta > 0.1 ~ "median_high_exprn_high_delta",
                                        T ~ "median_high_exprn_low_delta"
                                        )) 

cryp_type_med_df_counts <- cryp_type_med_df %>%
  count(median_cryp_status) %>%
  mutate(frac = n / sum(n))



# prepare base df for plotting
highlight_genes <- c("ELK1", "ARHGAP32", "STMN2", "CNPY3", "TLX1", "ANKRD27")
plot_med_df <- get_scatter_df(cryp_any_med_df, highlight_genes)

med_scatter <- plot_med_df %>%
  ggplot(aes(x = median_ctl*100,
             y = median_delta*100, 
             alpha = plot_alpha,
             colour = plot_colour, label = plot_name)) + 
  geom_point() +
  geom_hline(yintercept = -10, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  scale_colour_manual(values = c("#000000", "#d95f02")) +
  geom_text_repel(max.overlaps = 700, 
                  size = rel(8),
                  force = 60,
                  force_pull = 0.75,
                  direction = "both",
                  min.segment.length = 0,
                  seed = 123,
                  xlim = c(0,70),
                  ylim = c(-40, 100)
  ) +
  scale_x_continuous(breaks = seq(0,100,10)) + 
  scale_y_continuous(limits = c(-100,100),
                     breaks = seq(-100,100,10)) +
  labs(x = "Median of CTL mean PAS usage %",
       y = "Median of change in usage (TDP-43 KD - CTL)") +
  theme_bw(base_size = 16) + 
  guides(alpha = "none",
         colour = "none") +
  theme(axis.title.x = element_text(size = rel(1.75)),
        axis.title.y = element_text(size = rel(1.75)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5))
  )

med_scatter

# Alternative - plot gene names for all significant genes (to highlight mirroring)
# highlight the top n cryptic genes where separation is more obvious
highlight_genes_mirror <- cryp_any_med_df %>%
  filter(median_ctl < 0.1) %>%
  slice_max(median_delta, n = 10) %>%
  pull(gene_name)

plot_med_df_all_gn <- get_scatter_df(cryp_any_med_df, highlight_genes_mirror, label_all = T)
# plot_med_df_all_gn

med_scatter_all_gn <- plot_med_df_all_gn %>%
  ggplot(aes(x = median_ctl*100,
             y = median_delta*100, 
             alpha = plot_alpha,
             colour = plot_colour, label = plot_name)) + 
  geom_point() +
  geom_hline(yintercept = -10, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  scale_colour_manual(values = c("#000000", "#d95f02")) +
  geom_text_repel(max.overlaps = 1000, 
                  size = rel(6),
                  force = 80,
                  force_pull = 0.5,
                  direction = "both",
                  min.segment.length = 0,
                  seed = 123,
                  xlim = c(0,100)
                  ) +
  scale_x_continuous(breaks = seq(0,100,10)) + 
  scale_y_continuous(limits = c(-100,100),
                     breaks = seq(-100,100,10)) +
  labs(x = "Median of CTL mean PAS usage %",
       y = "Median of change in usage (TDP-43 KD - CTL)") +
  theme_bw(base_size = 16) + 
  guides(alpha = "none",
         colour = "none") +
  theme(axis.title.x = element_text(size = rel(1.75)),
        axis.title.y = element_text(size = rel(1.75)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5))
  )

med_scatter_all_gn


# first, get counts of unique IDs per gene
cryp_id_counts <- cryp_any_med_df %>%
  distinct(gene_id, gene_name, le_id) %>%
  count(gene_name, sort = T)

# visualise isoform counts across genes
plot_cryp_id_counts <- cryp_id_counts %>%
  ggplot(aes(x = n)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  scale_x_continuous(breaks = seq(1,10,1)) +
  scale_y_continuous(limits = c(0,2250)) +
  theme_bw(base_size = 14) +
  labs(x = "Unique Isoform Count",
       y = "Gene Count")

plot_cryp_id_counts


# Alternative plot - highlight events that are cryptic in at least one dataset (so can still highlight these, useful for later discussion)
x <- cryp_any_med_df %>%
  mutate(plot_name = if_else(cryptic_any, gene_name, ""),
         # 
         plot_alpha = case_when(plot_name != "" ~ 1, # keep all cryptics as 1 alpha
                                abs(median_delta) > 0.1 & median_ctl < 0.1 ~ 0.5,
                                abs(median_delta) > 0.1 & median_ctl > 0.1 ~ 0.2,
                                TRUE ~ 0.01),
         plot_colour = case_when(cryptic_any & cryptic_med ~ "orange",
                                 cryptic_any & !cryptic_med ~ "purple",
                                 T ~ "grey"
                                 )
         )

# 
# x %>%
#   ggplot(aes(x = median_ctl*100,
#              y = median_delta*100, 
#              alpha = plot_alpha,
#              colour = plot_colour, label = plot_name)) + 
#   geom_point() +
#   geom_hline(yintercept = -10, linetype = "dashed") +
#   geom_hline(yintercept = 10, linetype = "dashed") +
#   geom_vline(xintercept = 10, linetype = "dashed") +
#   scale_colour_manual(values = c("#000000", "#d95f02", "#7570b3"))


# old purple - "#7570b3"
colours_3grp <- c("Not cryptic" = "darkgrey", "Cryptic by medians" = "#d95f02", "Cryptic in >= 1 dataset" = "#1b9e77")

med_scatter_lab_1dataset <- ggplot(filter(x, !cryptic_any & !cryptic_med), aes(x = median_ctl*100,
                                                   y = median_delta*100,
                                                   alpha = plot_alpha,
                                                   colour = "Not cryptic"
                                                   )
       ) +
  geom_point() +
  geom_point(data = filter(x, cryptic_any & !cryptic_med),
             aes(x = median_ctl*100,
                 y = median_delta*100,
                 alpha = plot_alpha,
                 colour = "Cryptic in >= 1 dataset")) +
  geom_point(data = filter(x, cryptic_any & cryptic_med),
             aes(x = median_ctl*100,
                 y = median_delta*100,
                 alpha = plot_alpha,
             colour = "Cryptic by medians")
             ) + 
  scale_color_manual(values = colours_3grp) +
  geom_hline(yintercept = -10, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,100,10)) + 
  scale_y_continuous(limits = c(-100,100),
                     breaks = seq(-100,100,10)) +
  labs(x = "Median of CTL mean PAS usage %",
       y = "Median of change in usage (TDP-43 KD - CTL)",
       colour = "") +
  theme_bw(base_size = 16) + 
  guides(alpha = "none"
         ) +
  theme(axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.2)),
        legend.position = "top",
        legend.text = element_text(size = rel(1.25))
  )


med_scatter_lab_1dataset

if (!dir.exists("processed")) {dir.create("processed", recursive = T)}

ggsave(filename = "2024-11-17_tdp_kd_collection_cryptics_scatter_colour_medians_only_gene_name_mirror_eg.png",
       plot = med_scatter_all_gn,
       path = "processed",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-17_tdp_kd_collection_cryptics_scatter_colour_medians_only_gene_name_mirror_eg.pdf",
       plot = med_scatter_all_gn,
       path = "processed",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2024-11-17_tdp_kd_collection_sig_isoform_gene_counts_bar.png",
       plot = plot_cryp_id_counts,
       path = "processed",
       width = 125,
       height = 75,
       units = "mm",
       dpi = "retina")

ggsave(filename = "2024-11-17_tdp_kd_collection_sig_isoform_gene_counts_bar.pdf",
       plot = plot_cryp_id_counts,
       path = "processed",
       width = 125,
       height = 75,
       units = "mm",
       dpi = "retina")


ggsave(filename = "2023-10-02_tdp_kd_collection_cryptics_scatter_colour_medians_only_gene_name.png",
       plot = med_scatter,
       path = "processed",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-10-02_tdp_kd_collection_cryptics_scatter_colour_medians_only_gene_name.svg",
       plot = med_scatter,
       path = "processed",
       device = svg,
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-10-02_tdp_kd_collection_cryptics_scatter_colour_medians_only_gene_name.pdf",
       plot = med_scatter,
       path = "processed",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")



ggsave(filename = "2023-10-06_tdp_kd_collection_cryptics_scatter_colour_any_cryptic_no_gene_name.png",
       plot = med_scatter_lab_1dataset,
       path = "processed",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-10-06_tdp_kd_collection_cryptics_scatter_colour_any_cryptic_no_gene_name.svg",
       plot = med_scatter_lab_1dataset,
       path = "processed",
       device = svg,
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")


# write tsv of median values for each cryptic/regulated evetn
plot_med_df %>%
  select(le_id, groupID, gene_name, contains("median"), contains("cryptic"), simple_event_type, contains("plot")) %>%
  write_tsv("processed/2023-09-15_cryptics_scatter_standard_plot_tbl.tsv", col_names = T)

cryp_any_not_med_df %>%
  write_tsv("processed/2023-09-15_cryptics_1dataset_not_median_base_delta_tbl.tsv", col_names = T)

cryp_type_med_df %>%
  write_tsv("processed/2023-10-12_cryptic_median_exprn_category_tbl.tsv", col_names = T)

cryp_type_med_df_counts %>%
  write_tsv("processed/2023-10-12_cryptic_median_exprn_category_counts.tsv", col_names = T)

# save to Rdata
save(med_scatter_lab_1dataset, med_scatter, plot_med_df, cryp_any_not_med_df, sig_med_df, cryp_type_med_df_counts, cryp_type_med_df, file = "processed/2024-11-14_fig1_cryptics_scatter.Rdata")



