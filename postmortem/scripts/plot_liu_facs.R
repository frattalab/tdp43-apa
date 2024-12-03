library(tidyverse)
source("scripts/utils_liu_facs.R")

ppau_delta_paired_median_all <- read_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_delta_ppau.all_samples.all_ales.tsv.gz")
ppau_delta_paired_median_disease <- read_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_delta_ppau.subtype_split.all_ales.tsv.gz")
ppau_delta_paired <- read_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_per_sample_delta_ppau.all_ales.tsv.gz")
ppau_delta_mean_cryp <- read_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_popn_delta_ppau.cryptics.tsv")


# get an order of le_ids according to ranks of median delta between negative and positive
# (this puts events with most consistent enrichment at the top of the heatmap)
ppau_order_cryp <- ppau_delta_paired_median_all %>%
  filter(cryptic_status) %>%
  # rank deltas from largest to smallest, giving same rank if identical delta
  mutate(rank_mean_up = min_rank(desc(mean_paired_delta_ppau)),
         rank_median_up = min_rank(desc(median_paired_delta_ppau))) %>%
  arrange(rank_median_up, rank_mean_up)

# get ranked order of IDs from most enriched to least enriched 
ppau_order_cryp_median_gn <- ppau_order_cryp %>%
  arrange(rank_median_up) %>%
  pull(plot_le_id)

# get a df of events passing the median threshold
ppau_order_cryp_min5 <- ppau_order_cryp %>%
  filter(median_paired_delta_ppau > 0.05)

# 
delta_paired_summary <- ppau_delta_paired_median_all %>%
  filter(cryptic_status) %>%
  summary(paired_delta_ppau_neg_pos)

# List of IDs with at least 5 % median delta
ppau_order_cryp_median_gn_5 <- ppau_order_cryp %>%
  filter(median_paired_delta_ppau > 0.05) %>%
  pull(plot_le_id)

# Decay of number of enriched events as increase the delta cut-off
cryp_median_min_delta_n <- seq(0,0.25,0.05) %>%
  set_names() %>%
  map(~ ppau_order_cryp %>%
        filter(median_paired_delta_ppau > .x) %>%
        summarise(n_events = n_distinct(plot_le_id))
  ) %>%
  bind_rows(.id = "min_median_delta_pau")

# Decay of number of enriched events as increase the delta cut-off with counts split by event type
cryp_median_min_delta_nevent_by_type <- seq(0,0.25,0.05) %>%
  set_names() %>%
  map(~ ppau_order_cryp %>%
        mutate(enriched = median_paired_delta_ppau > .x) %>%
        group_by(simple_event_type) %>%
        summarise(n_enriched = sum(enriched),
                  n = n_distinct(plot_le_id),
        )
  ) %>%
  bind_rows(.id = "median_delta_cutoff")

# Generate plot-ready df - all events with at least +5 % delta median paired PPAU
# Order IDs in descending order of enrichment in each category
ppau_delta_paired_cryp <- ppau_delta_paired %>%
  filter(cryptic_status)

plot_df_median_5_delta <- ppau_delta_paired_cryp %>%
  filter(plot_le_id %in% ppau_order_cryp_median_gn_5) %>%
  prep_heatmap_df(le_id_order = rev(ppau_order_cryp_median_gn_5))

median5_std_events_heatmap <- plot_df_median_5_delta %>%
  filter(plot_event_type %in% c("ALE", "IPA", "3'Ext")) %>%
  plot_heatmap(fill_name ="Delta PAS usage %\n(TDPnegative - TDPpositive)",
               plot_labs = labs(x = "",
                                y = "",
                                ),
               plot_theme = theme_bw(base_size = 14),
               theme_args = list(
                 legend.position = "top",
                 legend.key.size = unit(20*0.75, "mm"),
                 legend.key.width = unit(25*0.75, "mm"),
                 legend.key.height = unit(10*0.75, "mm")
               ))

median5_std_events_heatmap

# heatmaps for other categories passing enrichment criteria
median5_other_events_heatmap <- plot_df_median_5_delta %>%
  filter(!plot_event_type %in% c("ALE", "IPA", "3'Ext")) %>%
  plot_heatmap()

median5_other_events_heatmap


# # Heatmap for NYGC selective events that aren't found in FACS-seq (manually curated)
# nygc_sel_only <- c("SYNJ2", "PHF2", "SERGEF", "PATJ", "DLGAP1") %>%
#   fct_inorder()
# 
# plot_df_nygc_sel_only <- ppau_delta_paired_cryp %>%
#   filter(gene_name %in% nygc_sel_only) %>%
#   prep_heatmap_df(le_id_order = rev(levels(nygc_sel_only)))
# 
# nygc_sel_only_heatmap <- plot_heatmap(plot_df_nygc_sel_only) +
#   geom_text(aes(label = round(paired_delta_ppau_neg_pos, 2))) +
#   scale_fill_gradient2(name = "Delta polyA usage % (TDPnegative - TDPpositive)",
#                        low = "#998ec3",
#                        mid = "#f7f7f7",
#                        high = "#f1a340",
#                        midpoint = 0,
#                        # limits = c(-1, 1),
#                        breaks = seq(-20, 20, 5)
#   )


if (!dir.exists("processed/liu_facs")) {dir.create("processed/liu_facs", recursive = T)}

# standard heatmaps
ggsave("2024-11-26_liu_facs_cryptic_median_delta_05_facet.std_event_types.png",
       plot = median5_std_events_heatmap,
       path = "processed/liu_facs",
       height = 150,
       width = 200,
       units = "mm",
       dpi = "retina")
ggsave("2024-11-26_liu_facs_cryptic_median_delta_05_facet.std_event_types.pdf",
       plot = median5_std_events_heatmap,
        path = "processed/liu_facs",
        height = 150,
        width = 200,
        units = "mm",
        dpi = "retina")

# other event types
ggsave("2024-11-21_liu_facs_cryptic_median_delta_05_facet.other_event_types.png",
       plot = median5_other_events_heatmap,
       path = "processed/liu_facs",
       height = 8,
       width = 12,
       units = "in",
       dpi = "retina")

ggsave("2024-11-21_liu_facs_cryptic_median_delta_05_facet.other_event_types.svg",
       plot = median5_other_events_heatmap,
       path = "processed/liu_facs",
       height = 10,
       width = 15,
       device = svg,
       units = "in",
       dpi = "retina")

ggsave("2024-11-21_liu_facs_cryptic_median_delta_05_facet.other_event_types.pdf",
       plot = median5_other_events_heatmap,
       path = "processed/liu_facs",
       height = 10,
       width = 15,
       units = "in",
       dpi = "retina")

# NYGC selective only
# ggsave("2024-11-21_liu_facs_nygc_selective_only.ales.png",
#        plot = nygc_sel_only_heatmap,
#        path = "processed/liu_facs",
#        height = 8,
#        width = 12,
#        units = "in",
#        dpi = "retina")
# 
# ggsave("2024-11-21_liu_facs_nygc_selective_only.ales.pdf",
#        plot = nygc_sel_only_heatmap,
#        path = "processed/liu_facs",
#        height = 8,
#        width = 12,
#        units = "in",
#        dpi = "retina")

# enriched count tables
write_tsv(cryp_median_min_delta_nevent_by_type, "processed/liu_facs/2024-11-21_liu_facs_min_delta_range_numevents_eventtype.tsv", col_names = T)
write_tsv(cryp_median_min_delta_n, "processed/liu_facs/2024-11-21_liu_facs_min_delta_cutoffs_numevents_all.tsv", col_names = T)
# enriched table (with ranks and summarised ppau)
write_tsv(ppau_order_cryp_min5, "processed/liu_facs/2024-11-22_liu_facs_median_delta_5.delta_ppau.all_samples.cryptics.tsv", col_names = T)
          