library(tidyverse)


path_summ_ale <- read_tsv("processed/nygc/expression_by_pathology_ale_all.tsv")
outdir <- "processed/nygc"

#
liu_specific_gn <- c("PCBP3", "ONECUT1", "PRKAR1B",
                     "RNF144A", "ZNF771", "ADARB2", "KIF26B",
                     "ACOT11", "DNM1P35", "HS6ST3", "SHLD1",
                     "GTF2IRD2", "PTPRN2", "NEK4", "GPR146", "PLXDC2", "ARHGAP32")

path_summ_ale_liu <- filter(path_summ_ale, gene_name %in% liu_specific_gn)


# prepare df for plotting
# Calculate path fold enrichment (for ordering the plot)
path_summ_ale_liu <-  path_summ_ale_liu %>%
  mutate(fold_enrichment = (fraction_path + 0.001) / (fraction_not_path + 0.001))

# make an additional, even simpler IDs just gene name suffixed with an event number (based on genomic start-end)
plot_path_summ_ale_liu <- path_summ_ale_liu %>%
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
  # now generate simple ID
  mutate(plot_name_simple = if_else(n_events == 1, gene_name, paste(gene_name, event_number, sep = "_")),
         plot_name_simple = fct_reorder(plot_name_simple, fold_enrichment)
         ) %>%
  select(-event_number, -n_events)

# longer format for plotting
plot_path_summ_ale_liu <- plot_path_summ_ale_liu %>%
  pivot_longer(cols = all_of(c("fraction_not_path", "fraction_path")),
               names_to = "path_status",
               values_to = "fraction_detected",
               names_prefix = "^fraction_") %>%
  mutate(plot_path = if_else(path_status == "path", "TDP-43 pathology", "No pathology"))

# bar plot of fraction detected in each group, with sorting by enrichment ratio
liu_specific_bar <- plot_path_summ_ale_liu %>%
  ggplot(aes(x = fraction_detected * 100, y = plot_name_simple, fill = plot_path)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#999999", "#E69F00"), labels = c("False", "True")) +
  theme_bw(base_size = 10) +
  labs(x = "Detected tissues (%)",
       y = "",
       fill = "TDP-43 Pathology") +
  theme(legend.position = "bottom")


liu_specific_bar

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

ggsave(filename = "2024-11-22_liu_enriched_specific.detection_bar.png",
       plot = liu_specific_bar,
       path = outdir,
       height = 75,
       width = 125,
       units = "mm",
       dpi = "retina"
       )

ggsave(filename = "2024-11-22_liu_enriched_specific.detection_bar.pdf",
       plot = liu_specific_bar,
       path = outdir,
       height = 75,
       width = 125,
       units = "mm",
       dpi = "retina"
)

