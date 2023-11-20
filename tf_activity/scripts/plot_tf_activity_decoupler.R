library(tidyverse)


#####

bind_rows(dorothea_all = ferguson_ulm_dorothea_elks,
          dorothea_nospl = ferguson_nospl_ulm_dorothea_elks,
          collectri_all = ferguson_ulm_collectri_elks,
          collectri_nospl = ferguson_nospl_ulm_collectri_elks,
          .id = "targets_source") %>%
  filter(source == "ELK1") %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  ggplot(aes(x = targets_source, y = score, label = round(p_value, 4))) +
  geom_col() +
  geom_text() + 
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs("univariate linear model Ferguson HeLa")