library(tidyverse)
library(ggrepel)
set.seed(123)


corticali3_kinetic <- read_tsv("processed/2023-08-22_i3cortical_slamseq_grandr_kinetics.tsv")

# ranks GW
corticali3_kinetic %>%
  mutate(hl_rank = rank(desc(log2FoldHalfLife)),
         synth_rank = rank(desc(log2FoldSynthesis))) %>%
  arrange(synth_rank) %>% View()

# calculate ranks specifically for up and down, or normalise to distance from middle possibly? 
# i.e. 1 = furthest above middle rank, -1 = furthest below middle rank?


# what are the relative differences 
corticali3_kinetic %>%
  filter(gene_name %in% c("ELK1", "SIX3", "TLX1")) %>%
  pivot_longer(cols = starts_with("log2Fold"),
               names_prefix = "^log2Fold",
               names_to = "metric",
               values_to = "log2FoldKDvsWT") %>%
  ggplot(aes(x = gene_name, y = log2FoldKDvsWT, colour = metric)) + 
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(limits = c(-1.5, 1.5),
                     breaks = seq(-1.5, 1.5, 0.25)
                     ) +
  scale_color_manual(values = c("#d95f02", "#1b9e77")) +
  theme_bw(base_size = 16) +
  labs(title = "Cryptic 3'UTRs - Stability vs synthesis",
       x = "Gene Name",
       )

# only ELK1 has a dramatic difference in stability and synthesis rate (opposing direction - preventing too high levels?)
# SIX3 & TLX1 have comparable increases in stabiluty and translation - can't say stability is driver of increased RNA & protein levels
ggsave(filename = "2023-08-22_cryptic_3utr_delta_stability_synthesis_point.png",
       path = "processed/",
       device = "png",
       dpi = "retina",
       height = 8,
       width = 8,
       units = "in")


# ELK1's dramatic difference, how does this look transcriptome wide?
plot_corticali3_synth_hl <- corticali3_kinetic %>%
  mutate(plot_colour = if_else(gene_name %in% c("ELK1", "SIX3", "TLX1"),
                               T,
                               F),
         plot_label = if_else(plot_colour, gene_name, "")
         )

plot_corticali3_synth_hl %>%
  filter(!plot_colour) %>%
  ggplot(aes(x = log2FoldHalfLife,
             y = log2FoldSynthesis,
             colour = plot_colour,
             label = plot_label
             ))+ 
  geom_point(alpha = 0.5) +
  geom_point(data = filter(plot_corticali3_synth_hl, plot_colour)) +
  scale_colour_manual(values = c("#bdbdbd", "#d95f02")) +
  geom_text_repel(data = filter(plot_corticali3_synth_hl, plot_colour),
                  force = 50,max.overlaps = 100, size = rel(5)) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(-8, 8)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw(base_size = 14) +
  guides(alpha = "none", colour = "none")

# SIX3 & TLX1 certainly represent an uncommon category of changes - increased stability and synthesis are rare
# most common is decreased stability & increased synthesis (likely a compensatory mechanism, possibly also true for the likes of ELK1?)
ggsave(filename = "2023-08-22_txome_wide_stability_synthesis_scatter.png",
       path = "processed/",
       device = "png",
       dpi = "retina",
       height = 8,
       width = 8,
       units = "in")


