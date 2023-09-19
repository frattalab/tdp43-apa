library(tidyverse)


path_summ_ale <- read_tsv("processed/nygc/expression_by_pathology_ale_all.tsv")

# Classify as selective or enriched according to AL's criteria
path_summ_ale <- path_summ_ale %>%
  mutate(selective = fraction_not_path <= 0.005 & fraction_path >= 0.01,
         enriched = !selective & (fraction_path > fraction_not_path * 2.5))

count(path_summ_ale, selective, enriched)
# A tibble: 3 × 3
# selective enriched     n
# <lgl>     <lgl>    <int>
#   1 FALSE     FALSE       48
# 2 FALSE     TRUE        60
# 3 TRUE      FALSE        7


path_summ_sel <- filter(path_summ_ale, selective)
path_summ_enr <- filter(path_summ_ale, enriched)



# create name as junction coords - gene_name
selective_jnc <- path_summ_sel %>% pull(paste_into_igv_junction) %>%
  set_names(paste(path_summ_sel$paste_into_igv_junction, path_summ_sel$gene_name, sep = "_") %>% str_replace_all(":|-", "_"))

enriched_jnc <- path_summ_enr %>% pull(paste_into_igv_junction) %>%
  set_names(paste(path_summ_enr$paste_into_igv_junction, path_summ_enr$gene_name, sep = "_") %>% str_replace_all(":|-", "_"))

enriched_rpm <- ale_enriched_by_rpm %>%
  pull(paste_into_igv_junction) %>%
  set_names(paste(ale_enriched_by_rpm$paste_into_igv_junction, ale_enriched_by_rpm$gene_name, sep = "_") %>% str_replace_all(":|-", "_"))

sel_jnc_plots <- map(selective_jnc,
                     ~ plot_junction(.x))

enr_jnc_plots <- map(enriched_jnc,
                     ~ plot_junction(.x))

enr_rpm_jnc_plots <- map(enriched_rpm,
                         ~ plot_junction(.x))

dir.create("processed/nygc/selective_jncs/", recursive = T)
dir.create("processed/nygc/enriched_jncs/", recursive = T)
dir.create("processed/nygc/enriched_rpm_jncs/", recursive = T)
#
walk2(sel_jnc_plots,
      names(sel_jnc_plots),
      ~ ggsave(filename = glue::glue("processed/nygc/selective_jncs/2023-09-13_nygc_papa_as_ale.selective.spliced_reads.{.y}.png"),
               plot = .x,
               height = 14,
               width = 14,
               units = "in",
               dpi = "retina"
      )
)

walk2(enr_jnc_plots,
      names(enr_jnc_plots),
      ~ ggsave(filename = glue::glue("processed/nygc/enriched_jncs/2023-09-13_nygc_papa_as_ale.enriched.spliced_reads.{.y}.png"),
               plot = .x,
               height = 14,
               width = 14,
               units = "in",
               dpi = "retina"
      )
)

walk2(enr_rpm_jnc_plots,
      names(enr_rpm_jnc_plots),
      ~ ggsave(filename = glue::glue("processed/nygc/enriched_rpm_jncs/2023-09-13_nygc_papa_as_ale.enriched_rpm.spliced_reads.{.y}.png"),
               plot = .x,
               height = 14,
               width = 14,
               units = "in",
               dpi = "retina"
      )
)

# Make bar plot of %/fraction of path tissues detected
# do separately for specific and enriched
# Also a plot of combined


# just selective, ordered by most tissues obs
path_summ_sel %>%
  mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
         plot_name = fct_reorder(plot_name, fraction_path)) %>%
  ggplot(aes(x = fraction_path * 100, y = plot_name)) +
  geom_col() +
  theme_bw(base_size = 14) +
  labs(title = "Selective AS-ALE junctions",
       x = "% of TDP-43 pathological tissues cryptic detected",
       y = "Junction|Gene name")

# just selective, include fraction non_path
path_summ_sel %>%
  mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
         plot_name = fct_reorder(plot_name, fraction_path)) %>%
  pivot_longer(starts_with("fraction"), names_to = "tissue_group", values_to = "fraction") %>%
  ggplot(aes(x = fraction * 100, y = plot_name, fill = tissue_group)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
  theme_bw(base_size = 14) +
  labs(title = "Selective AS-ALE junctions",
       x = "% of TDP-43 pathological tissues cryptic detected",
       y = "Junction|Gene name")


# just enriched, include frac non-path
# (potentially a bit misleading if expressed in a lot of tissues, could be an expression change rather than strong cryp expression - check PSI)
path_summ_enr %>%
  mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
         plot_name = fct_reorder(plot_name, fraction_path)) %>%
  pivot_longer(starts_with("fraction"), names_to = "tissue_group", values_to = "fraction") %>%
  ggplot(aes(x = fraction * 100, y = plot_name, fill = tissue_group)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
  theme_bw(base_size = 14) +
  labs(title = "Enriched AS-ALE junctions",
       subtitle = "Enriched = fraction path > 2.5* fraction non-path",
       x = "% of TDP-43 pathological tissues cryptic detected",
       y = "Junction|Gene name")

bind_rows(Enriched = path_summ_enr, Selective = path_summ_sel, .id = "sel_group") %>%
  mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
         plot_name = fct_reorder(plot_name, fraction_path),
         sel_group = factor(sel_group, levels = c("Selective", "Enriched"))) %>%
  pivot_longer(starts_with("fraction"), names_to = "tissue_group", values_to = "fraction") %>%
  ggplot(aes(x = fraction * 100, y = plot_name, fill = tissue_group)) +
  facet_wrap("~ sel_group", ncol = 1, scales = "free_y") +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
  theme_bw(base_size = 10) +
  labs(title = "Selective & Enriched AS-ALE junctions",
       x = "% of TDP-43 pathological tissues cryptic detected",
       y = "Junction|Gene name") +
  theme(axis.text.y = element_text(size = rel(0.65)))


#enriched by rpm
ale_enriched_by_rpm %>%
  filter(!paste_into_igv_junction %in% path_summ_sel$paste_into_igv_junction) %>%
  mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
         plot_name = fct_reorder(plot_name, fraction_path)) %>%
  pivot_longer(starts_with("fraction"), names_to = "tissue_group", values_to = "fraction") %>%
  ggplot(aes(x = fraction * 100, y = plot_name, fill = tissue_group)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
  theme_bw(base_size = 14) +
  labs(title = "Enriched AS-ALE junctions",
       subtitle = "Enriched by wilcoxon test of spliced RPM",
       x = "% of TDP-43 pathological tissues cryptic detected",
       y = "Junction|Gene name")