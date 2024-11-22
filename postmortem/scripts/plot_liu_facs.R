library(tidyverse)

# means and median of pairwise PPAUs across all samples and split by subtype
ppau_delta_paired_median_all <- read_tsv("processed/2023-09-11_liu_facs_decoys_delta_ppau.all_samples.all_ales.tsv.gz")
ppau_delta_paired_median_disease <- read_tsv("processed/2023-09-11_liu_facs_decoys_delta_ppau.subtype_split.all_ales.tsv.gz")
ppau_delta_paired <- read_tsv("processed/2023-09-11_liu_facs_decoys_per_sample_delta_ppau.all_ales.tsv.gz")

# get an order of le_ids according to ranks of median delta between negative and positive
# (this puts events with most consistent enrichment at the top of the heatmap)
ppau_order_cryp <- ppau_delta_paired_median_all %>%
  filter(cryptic_status) %>%
  # rank deltas from largest to smallest, giving same rank if identical delta
mutate(rank_mean_up = min_rank(desc(mean_paired_delta_ppau)),
       rank_median_up = min_rank(desc(median_paired_delta_ppau))) %>%
  arrange(rank_median_up, rank_mean_up)

# get a df of events passing the median threshold
ppau_order_cryp_min5 <- ppau_order_cryp %>%
  filter(median_paired_delta_ppau > 0.05)



# ppau_order_cryp_median <- ppau_order_cryp %>%
#   arrange(rank_median_up) %>%
#   pull(le_id)

ppau_order_cryp_median_gn <- ppau_order_cryp %>%
  arrange(rank_median_up) %>%
  pull(plot_le_id)

ppau_delta_paired_median_all %>%
  filter(cryptic_status) %>%
  summary(paired_delta_ppau_neg_pos)

# rejoin gene information 
# ppau_delta_paired_median_all <- ppau_delta_paired_median_all %>%
#   left_join(distinct(ppau_delta_paired, le_id, gene_id, gene_name), by = "le_id") %>%
#   relocate(all_of(c("mean_paired_delta_ppau", "median_paired_delta_ppau",
#                     "rank_mean_up", "rank_median_up")),
#            .after = last_col()
#            )


ppau_delta_paired_cryp <- ppau_delta_paired %>%
  filter(cryptic_status)

# add event
# at least 5 % median delta

ppau_order_cryp_median_gn_5 <- ppau_order_cryp %>%
  filter(median_paired_delta_ppau > 0.05) %>%
  pull(plot_le_id)

# decay of number of enriched events as increase the delta cut-off
cryp_median_min_delta_n <- seq(0,0.25,0.05) %>%
  set_names() %>%
  map(~ ppau_order_cryp %>%
                               filter(median_paired_delta_ppau > .x) %>%
        summarise(n_events = n_distinct(plot_le_id))
  ) %>%
  bind_rows(.id = "min_median_delta_pau")


# Get the number of enriched events by event type at delta cut-off
cryp_median_min_5_nevent_by_type <- ppau_order_cryp %>%
  left_join(distinct(ppau_delta_paired_cryp, plot_le_id, simple_event_type), by = "plot_le_id") %>%
  mutate(enriched = median_paired_delta_ppau > 0.05) %>%
  group_by(simple_event_type) %>%
  summarise(n_enriched = sum(enriched),
            n = n_distinct(plot_le_id),
            )

cryp_median_min_delta_nevent_by_type <- seq(0,0.25,0.05) %>%
  set_names() %>%
  map(~ ppau_order_cryp %>% # lazy
        left_join(distinct(ppau_delta_paired_cryp, plot_le_id, simple_event_type), by = "plot_le_id") %>%
        mutate(enriched = median_paired_delta_ppau > .x) %>%
        group_by(simple_event_type) %>%
        summarise(n_enriched = sum(enriched),
                  n = n_distinct(plot_le_id),
        )
        ) %>%
  bind_rows(.id = "median_delta_cutoff")
  

  
####
plot_df_median_5_delta <- ppau_delta_paired_cryp %>%
  filter(plot_le_id %in% ppau_order_cryp_median_gn_5) %>%
  mutate(
    plot_le_id = factor(plot_le_id, levels = rev(ppau_order_cryp_median_gn_5)),
    # simple_event_type = factor(simple_event_type, levels = c("spliced", "bleedthrough", "distal_3utr_extension",
    #                                                          "bleedthrough,spliced",
    #                                                          "bleedthrough,distal_3utr_extension",
    #                                                          NA
    #                                                          )
    #                            ),
    plot_event_type = case_when(simple_event_type == "spliced" ~ "AS-ALE",
                                simple_event_type == "bleedthrough" ~ "Bleedthrough-ALE",
                                simple_event_type == "distal_3utr_extension" ~ "3'UTR-ALE",
                                str_detect(simple_event_type, ",") ~ "Complex"),
    plot_event_type = factor(plot_event_type, levels = c("AS-ALE", "Bleedthrough-ALE", "3'UTR-ALE", "Complex", NA)),
    plot_patient_id = factor(paste("Sample_", str_extract(patient_id, "\\d$"), sep = ""), levels = c("Sample_2", "Sample_6", "Sample_7", "Sample_1", "Sample_3", "Sample_4", "Sample_5"))
  )

plot_df_median_5_delta %>%
  ggplot(aes(x = plot_patient_id,
             y = plot_le_id,
             fill = paired_delta_ppau_neg_pos)) +
  facet_wrap("~ plot_event_type", scales = "free_y") +
  geom_tile() +
  scale_fill_gradientn(name = "Delta polyA usage (TDPnegative - TDPpositive)",
                       colours = c("#998ec3", "#f7f7f7", "#f1a340"),
                       limits = c(-1, 1),
                       breaks = seq(-1, 1, 0.2)) +
  theme_bw(base_size = 20) +
  labs(x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "top",
        legend.key.size = unit(20, "mm"),
        legend.key.width = unit(25, "mm"),
        legend.key.height = unit(10, "mm") 
        )

if (!dir.exists("processed/liu_facs")) {dir.create("processed/liu_facs", recursive = T)}

ggsave("2023-10-18_liu_facs_cryptic_median_delta_05_event_type_facet.png",
       path = "processed/liu_facs",
       height = 8,
       width = 12,
       units = "in",
       dpi = "retina")


ggsave("2023-10-18_liu_facs_cryptic_median_delta_05_event_type_facet.svg",
       path = "processed/liu_facs",
       height = 10,
       width = 15,
       device = svg,
       units = "in",
       dpi = "retina")

write_tsv(cryp_median_min_5_nevent_by_type, "processed/liu_facs/2023-10-18_liu_facs_min_delta_5_numevents_eventtype.tsv", col_names = T)
write_tsv(cryp_median_min_delta_nevent_by_type, "processed/liu_facs/2023-10-18_liu_facs_min_delta_range_numevents_eventtype.tsv", col_names = T)
write_tsv(cryp_median_min_delta_n, "processed/liu_facs/2023-10-18_liu_facs_min_delta_cutoffs_numevents_all.tsv", col_names = T)
write_tsv(ppau_order_cryp_min5, "processed/liu_facs/2023-10-18_liu_facs_median_delta_5.delta_ppau.all_samples.cryptics.tsv", col_names = T)
# ppau_delta_paired_cryp %>%
#   filter(plot_le_id %in% ppau_order_cryp_median_gn_5) %>%
#   mutate(
#     plot_le_id = factor(plot_le_id, levels = rev(ppau_order_cryp_median_gn_5))
#   ) %>%
#   ggplot(aes(x = patient_id,
#              y = plot_le_id,
#              fill = paired_delta_ppau_neg_pos)) +
#   geom_tile() +
#   scale_fill_gradientn(name = "Sample-wise dPPAU (TDPneg - TDPpos)",
#                        colours = c("#998ec3", "#f7f7f7", "#f1a340"),
#                        limits = c(-1, 1)) +
#   theme_bw() +
#   labs(title = "Cell line cryptic last exons - relative usage change in FACS nuclei",
#        x = "Sample ID",
#        y = "Gene name") +
#   theme(axis.text.x = element_text(angle = 45,hjust = 1),
#         axis.text.y = element_text(size=rel(0.82)))


# ppau_delta_paired_cryp %>%
#   mutate(le_id = factor(le_id, levels = rev(ppau_order_cryp_median)),
#          plot_le_id = factor(plot_le_id, levels = rev(ppau_order_cryp_median_gn))
#          ) %>%
#   ggplot(aes(x = patient_id,
#              y = le_id,
#              fill = paired_delta_ppau_neg_pos)) +
#   geom_tile() +
#   scale_fill_gradientn(name = "Sample-wise dPPAU (TDPneg - TDPpos)",
#                        colours = c("#998ec3", "#f7f7f7", "#f1a340"),
#                        limits = c(-1, 1)) +
#   theme_bw() +
#   labs(title = "Cell line cryptic last exons - relative usage change in FACS nuclei",
#        x = "Sample ID",
#        y = "Gene name") +
#   theme(axis.text.x = element_text(angle = 45,hjust = 1),
#         axis.text.y = element_text(size=rel(0.82)))
# 
# ppau_delta_paired_cryp %>%
#   mutate(le_id = factor(le_id, levels = rev(ppau_order_cryp_median)),
#          plot_le_id = factor(plot_le_id, levels = rev(ppau_order_cryp_median_gn))
#   ) %>%
#   ggplot(aes(x = patient_id,
#              y = plot_le_id,
#              fill = paired_delta_ppau_neg_pos)) +
#   geom_tile() +
#   scale_fill_gradientn(name = "Sample-wise dPPAU (TDPneg - TDPpos)",
#                        colours = c("#998ec3", "#f7f7f7", "#f1a340"),
#                        limits = c(-1, 1)) +
#   theme_bw() +
#   labs(title = "Cell line cryptic last exons - relative usage change in FACS nuclei",
#        x = "Sample ID",
#        y = "Gene name") +
#   theme(axis.text.x = element_text(angle = 45,hjust = 1),
#         axis.text.y = element_text(size=rel(0.82)))
# 
# 






# use ggside to add gender (& possibly gene expression) as side tiles to heatmap, event type
# https://cran.r-project.org/web/packages/ggside/vignettes/ggside_basic_usage.html

