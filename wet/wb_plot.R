library(tidyverse)
# library(ggpubr)

# read in Halo-TDP quantification
df <- read_csv("data/2023_12_17_elk1_quantification_relative_normalised.csv")

# convert conditions in to single column
df <- df %>%
  pivot_longer(cols = -sample, names_to = "condition", values_to = "rel_expr") %>%
  # rename TDPKD to TDP43KD for consistency with current FISH plots
  mutate(plot_condition = if_else(condition == "TDPKD", "TDP43KD", "Control"))

rel_exp_scatter <- df %>%
  ggplot(aes(x = plot_condition, y = rel_expr)) +
  geom_point(position = position_jitter(width = 0.2, seed = 123), size = 2) +
  scale_y_continuous(breaks = seq(0,8,1)) +
  theme_bw(base_size = 20) +
  labs(x = "",
       y = "Relative Expression")
  
if (!dir.exists("processed/")) {dir.create("processed")}

ggsave("2023-12-19_i3n_tdp_halo_elk1_wb_scatter.png",
       plot = rel_exp_scatter,
       path = "processed/",
       device = "png",
       units = "in",
       height = 8,
       width = 8,
       dpi = "retina")

ggsave("2023-12-19_i3n_tdp_halo_elk1_wb_scatter.svg",
       plot = rel_exp_scatter,
       path = "processed/",
       device = svg,
       units = "in",
       height = 4,
       width = 4,
       dpi = "retina")



