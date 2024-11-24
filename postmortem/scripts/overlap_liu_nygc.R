library(tidyverse)
library(ggVennDiagram)
source("scripts/utils_liu_facs.R")

nygc_summ_ale <- read_tsv("processed/nygc/expression_by_pathology_ale_all.tsv")
nygc_sel_ale <- nygc_summ_ale %>%
  filter(selective)

liu_min5_all <- read_tsv("processed/liu_facs/2024-11-22_liu_facs_median_delta_5.delta_ppau.all_samples.cryptics.tsv")
liu_min5_ale <- liu_min5_all %>%
  filter(simple_event_type == "spliced")

# Liu per-sample PPAU
liu_ppau_delta_paired <- read_tsv("processed/liu_facs/2024-11-20_liu_facs_decoys_per_sample_delta_ppau.all_ales.tsv.gz") %>%
  filter(cryptic_status)

outdir <- "processed/liu_facs_nygc_overlap"


# map le_ids to gene names for later annotation
le2name <- bind_rows(distinct(liu_min5_ale, le_id, gene_name),
          distinct(nygc_sel_ale, le_id, gene_name))%>%
  distinct(.keep_all = T)

# list of IDs for overlap
le_id_hits_list <- list("NYGC Selective" = unique(nygc_sel_ale$le_id),
     "Liu Enriched" = unique(liu_min5_ale$le_id)
     )

# try ggVennDiagram
hits_venn_df <- process_region_data(Venn(le_id_hits_list))
hits_venn_plot <- ggVennDiagram(le_id_hits_list,
                                label = "count",  
                                label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "darkgrey", name = "Count") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(legend.position = "bottom")

hits_venn_plot

# Clean up and annotate interssection df
# row per ID & add gene name/plot id
hits_venn_df <- hits_venn_df %>%
  unnest(item) %>%
  rename(le_id = item) %>%
  left_join(le2name, by = "le_id")

hits_venn_df


# Liu specific - plot junction detection across path status in NYGC
liu_specific <- filter(hits_venn_df, name == "Liu Enriched") %>% pull(le_id)
nygc_specific <- filter(hits_venn_df, name == "NYGC Selective") %>% pull(le_id)

nygc_liu_ale <- filter(nygc_summ_ale, le_id %in% liu_specific)

# prepare df for plotting
# Calculate path fold enrichment (for ordering the plot)
nygc_liu_ale <-  nygc_liu_ale %>%
  mutate(fold_enrichment = (fraction_path + 0.001) / (fraction_not_path + 0.001))

# make an additional, even simpler IDs just gene name suffixed with an event number (based on genomic start-end)
plot_nygc_liu_ale <- nygc_liu_ale %>%
  mutate(coords = str_split_i(paste_into_igv_junction, ":", 2),
  ) %>% 
  separate(coords, into = c("start", "end"), sep = "-",convert = T) %>%
  # for each gene, assign a number based on left-right genomic position order in gene
  group_by(gene_name) %>%
  arrange(start, end, .by_group = T) %>%
  mutate(event_number = row_number(),
         n_events = n_distinct(paste_into_igv_junction)
  ) %>%
  ungroup() %>%
  # now generate simple ID & sort according to detection enrichment
  mutate(plot_name_simple = if_else(n_events == 1, gene_name, paste(gene_name, event_number, sep = "_")),
         plot_name_simple = fct_reorder(plot_name_simple, fold_enrichment)
  ) %>%
  select(-event_number, -n_events)

# longer format for plotting
plot_nygc_liu_ale <- plot_nygc_liu_ale %>%
  pivot_longer(cols = all_of(c("fraction_not_path", "fraction_path")),
               names_to = "path_status",
               values_to = "fraction_detected",
               names_prefix = "^fraction_") %>%
  mutate(plot_path = if_else(path_status == "path", "TDP-43 pathology", "No pathology"))

# bar plot of fraction detected in each group, with sorting by enrichment ratio
liu_specific_bar <- plot_nygc_liu_ale %>%
  ggplot(aes(x = fraction_detected * 100, y = plot_name_simple, fill = plot_path)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#999999", "#E69F00"), labels = c("False", "True")) +
  theme_bw(base_size = 10) +
  labs(x = "Detected tissues (%)",
       y = "",
       fill = "TDP-43 Pathology") +
  theme(legend.position = "bottom")

liu_specific_bar


# NYGC specific in Liu FACS
liu_nygc_ppau_delta_paired <- liu_ppau_delta_paired %>%
  filter(le_id %in% nygc_specific)

# calculate median pairwise delta PPAU (to define plotting order)
liu_nygc_median_ppau_delta_paired <- liu_nygc_ppau_delta_paired %>%
  group_by(le_id, plot_le_id) %>%
  summarise(median_paired_delta_ppau = median(paired_delta_ppau_neg_pos)) %>%
  ungroup() %>%
  mutate(rank_median_up = min_rank(desc(median_paired_delta_ppau))) %>%
  arrange(rank_median_up)

# prepare df for generating heatmap
plot_df_liu_nygc_ppau <- prep_heatmap_df(liu_nygc_ppau_delta_paired,
                le_id_order = rev(pull(liu_nygc_median_ppau_delta_paired,
                                       plot_le_id))
                )

# plot heatmap, modifying defaults to fit delta range of these events
nygc_sel_only_heatmap <- plot_heatmap(plot_df_liu_nygc_ppau,
             plot_theme = theme_bw(base_size = 10),
             theme_args = list(legend.position = "top")
             ) +
  geom_text(aes(label = round(paired_delta_ppau_neg_pos, 2)), size = 3) +
  scale_fill_gradient2(name = "Delta polyA usage %\n(TDPnegative - TDPpositive)",
                       low = "#998ec3",
                       mid = "#f7f7f7",
                       high = "#f1a340",
                       midpoint = 0,
                       # limits = c(-1, 1),
                       breaks = seq(-5, 15, 5)
)

nygc_sel_only_heatmap

# save plots and dfs to file

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

# Venn diagram
ggsave(filename = "2024-11-24_nygc_liu_venn.png",
       plot = hits_venn_plot,
       path = outdir,
       height = 125,
       width = 75,
       units = "mm",
       dpi = "retina"
)

ggsave(filename = "2024-11-24_nygc_liu_venn.pdf",
       plot = hits_venn_plot,
       path = outdir,
       height = 125,
       width = 75,
       units = "mm",
       dpi = "retina"
)

# Liu FACS specific
ggsave(filename = "2024-11-24_liu_enriched_specific.detection_bar.png",
       plot = liu_specific_bar,
       path = outdir,
       height = 75,
       width = 125,
       units = "mm",
       dpi = "retina"
)

ggsave(filename = "2024-11-24_liu_enriched_specific.detection_bar.pdf",
       plot = liu_specific_bar,
       path = outdir,
       height = 75,
       width = 125,
       units = "mm",
       dpi = "retina"
)

# NYGC specific
ggsave("2024-11-21_liu_facs_nygc_selective_only.ales.png",
       plot = nygc_sel_only_heatmap,
       path =outdir,
       height = 75,
       width = 125,
       units = "mm",
       dpi = "retina")

ggsave("2024-11-21_liu_facs_nygc_selective_only.ales.pdf",
       plot = nygc_sel_only_heatmap,
       path = outdir,
       height = 75,
       width = 125,
       units = "mm",
       dpi = "retina")




