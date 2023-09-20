library(tidyverse)
library(data.table)

plot_junction = function(junc,plotin_table = spliced_counts_ale){
  
  vals = c(`ALS-TDP` = "#E1BE6A", Control = "#40B0A6", `FTD-TDP` = "#E1BE6A", 
           `ALS\nnon-TDP` = "#408A3E", `FTD\nnon-TDP` = "#408A3E")
  
  gene_name = plotin_table[paste_into_igv_junction == junc,unique(gene_name)]
  
  plot_title = glue::glue("{gene_name} - {junc}")
  
  plt = plotin_table[paste_into_igv_junction == junc] |>
    filter(!tissue_clean %in% c("Other","Choroid","Hippocampus","Liver","Sensory_Cortex","Occipital_Cortex")) |> 
    mutate(disease =  gsub("-n","\nn",disease)) |> 
    mutate(disease = fct_relevel(disease,"Control", "ALS\nnon-TDP","ALS-TDP", "FTD\nnon-TDP", "FTD-TDP")) |>
    ggplot(aes(x = disease, y = spliced_reads, fill = disease)) +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(height = 0,alpha = 0.7, pch = 21) +
    facet_wrap(~tissue_clean, scales = "free_y") +
    scale_fill_manual(values = vals)  +
    ylab("N spliced reads") +
    xlab("") +
    ggtitle(plot_title) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 16)) 
  
  
  
  return(plt)
  
}



path_summ_ale <- read_tsv("processed/nygc/expression_by_pathology_ale_all.tsv")
path_summ_seddighi <- read_tsv("processed/nygc/expression_by_pathology_seddighi_all.tsv")
spliced_counts_ale <- fread("processed/nygc/spliced_counts_ale_all.tsv")
path_summ_seddighi_counts <- read_tsv("processed/nygc/expression_by_pathology_seddighi_counts.tsv")
path_summ_ale_counts <- read_tsv("processed/nygc/expression_by_pathology_ale_counts.tsv")


# Make Bar plot of frac detected
path_summ_comb <- bind_rows("AS-ALE" = path_summ_ale,
          "Other SJs" = path_summ_seddighi, .id = "event_type") 


# How many ALE junctions are identical to MAJIQ?
shared <- path_summ_comb %>%
  group_by(paste_into_igv_junction) %>%
  filter(n_distinct(event_type) > 1) %>%
  ungroup() %>%
  arrange(paste_into_igv_junction, event_type)

shared_ids <- unique(shared$paste_into_igv_junction)

# Where shared, prefer the AS-ALE event type
path_summ_comb <- path_summ_comb %>%
  group_by(paste_into_igv_junction) %>%
  arrange(event_type) %>%
  slice_head(n=1) %>%
  ungroup()

# assign a neater ID for plotting
path_summ_comb <- path_summ_comb %>%
  mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
         plot_name = fct_reorder(plot_name, fraction_path)
         )

# make an additional, even simpler IDs just gene name suffixed with an event number (based on genomic start-end)
path_summ_comb <- path_summ_comb %>%
  # "chr10:133166939-133167381" to start, end coords
  mutate(coords = str_split_i(paste_into_igv_junction, ":", 2)) %>% 
  separate(coords, into = c("start", "end"), sep = "-",convert = T) %>%
  # for each gene, assign a number based on left-right genomic position order in gene
  group_by(gene_name) %>%
  arrange(start, end, .by_group = T) %>%
  mutate(event_number = row_number()) %>%
  ungroup() %>%
  # now generate simple ID
  mutate(plot_name_simple = paste(gene_name, event_number, sep = "_"),
         plot_name_simple = fct_reorder(plot_name_simple, fraction_path))



# path_summ_comb %>%
#   mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
#          plot_name = fct_reorder(plot_name, fraction_path)
#   )



# bar plot of fraction observed coloured by SJ
sel_bar_jnc_gn <- path_summ_comb %>%
  filter(selective) %>%
  ggplot(aes(x = fraction_path * 100, y = plot_name, fill = event_type)) +
  geom_col() +
  scale_fill_manual(values = c("#d95f02", "#7570b3")) +
  theme_bw(base_size = 14) +
    labs(title = "Selective AS-ALE junctions",
         subtitle = "Non path detected fraction < 0.005, path detected > 0.01. Min spliced reads = 2",
         x = "% of TDP-43 pathological tissues cryptic detected",
         y = "Junction|Gene name",
         fill = "Event type")

sel_bar_gn <- path_summ_comb %>%
  filter(selective) %>%
  ggplot(aes(x = fraction_path * 100, y = plot_name_simple, fill = event_type)) +
  geom_col() +
  scale_fill_manual(values = c("#d95f02", "#7570b3")) +
  theme_bw(base_size = 14) +
  labs(title = "Selective AS-ALE junctions",
       subtitle = "Non path detected fraction < 0.005, path detected > 0.01. Min spliced reads = 2",
       x = "% of TDP-43 pathological tissues cryptic detected",
       y = "Event ID",
       fill = "Event type")
  
dir.create("processed/nygc/selective_jncs/ale/svg", recursive = T)
dir.create("processed/nygc/selective_jncs/ale/png", recursive = T)
dir.create("processed/nygc/enriched_jncs/ale/svg", recursive = T)
dir.create("processed/nygc/enriched_jncs/ale/png", recursive = T)

ggsave(filename = "2023-09-20_nygc_papa_seddighi_selective_jnc_gn_bar.png",
       plot = sel_bar_jnc_gn,
       path = "processed/nygc/",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-09-20_nygc_papa_seddighi_selective_gn_bar.png",
       plot = sel_bar_gn,
       path = "processed/nygc/",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-09-20_nygc_papa_seddighi_selective_jnc_gn_bar.svg",
       device = svg,
       plot = sel_bar_jnc_gn,
       path = "processed/nygc/",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

ggsave(filename = "2023-09-20_nygc_papa_seddighi_selective_gn_bar.svg",
       plot = sel_bar_gn,
       device = svg,
       path = "processed/nygc/",
       width = 12,
       height = 12,
       units = "in",
       dpi = "retina")

#
path_summ_ale_counts %>%
  mutate(frac = n / sum(n))
#
path_summ_seddighi_counts %>%
  mutate(frac = n / sum(n))


#### - Make box plots for each junction


sel_path_sum_ale <- path_summ_comb %>%
  filter(selective & event_type == "AS-ALE") 

enr_path_sum_ale <- path_summ_comb %>%
  filter(enriched & event_type == "AS-ALE") 

# create named lists of selective & enriched splice junctions
selective_jnc_ale <- sel_path_sum_ale %>%
  pull(paste_into_igv_junction) %>%
  set_names(paste(sel_path_sum_ale$paste_into_igv_junction, sel_path_sum_ale$gene_name, sep = "_") %>% 
              str_replace_all(":|-", "_")
            )

enriched_jnc_ale <- enr_path_sum_ale %>%
  pull(paste_into_igv_junction) %>%
  set_names(paste(enr_path_sum_ale$paste_into_igv_junction, enr_path_sum_ale$gene_name, sep = "_") %>% 
              str_replace_all(":|-", "_")
  )

# Make box plots of spliced reads for each junction 
sel_jnc_ale_plots <- map(selective_jnc_ale,
                     ~ plot_junction(.x))

enr_jnc_ale_plots <- map(enriched_jnc_ale,
                          ~ plot_junction(.x))





# Save splice read box plots to file

walk2(sel_jnc_ale_plots,
      names(sel_jnc_ale_plots),
      ~ ggsave(filename = glue::glue("2023-09-20_nygc_papa_as_ale.selective.spliced_reads.{.y}.png"),
               path = "processed/nygc/selective_jncs/ale/png/",
               plot = .x,
               height = 14,
               width = 14,
               units = "in",
               dpi = "retina"
      ),
      .progress = T
)

walk2(enr_jnc_ale_plots,
      names(enr_jnc_ale_plots),
      ~ ggsave(filename = glue::glue("2023-09-20_nygc_papa_as_ale.enriched.spliced_reads.{.y}.png"),
               path = "processed/nygc/enriched_jncs/ale/png/",
               plot = .x,
               height = 14,
               width = 14,
               units = "in",
               dpi = "retina"
      ),
      .progress = T
)

# save as SVGs
walk2(sel_jnc_ale_plots,
      names(sel_jnc_ale_plots),
      ~ ggsave(filename = glue::glue("2023-09-20_nygc_papa_as_ale.selective.spliced_reads.{.y}.svg"),
               path = "processed/nygc/selective_jncs/ale/svg/",
               device = svg,
               plot = .x,
               height = 14,
               width = 14,
               units = "in",
               dpi = "retina"
      ),
      .progress = T
)

walk2(enr_jnc_ale_plots,
      names(enr_jnc_ale_plots),
      ~ ggsave(filename = glue::glue("2023-09-20_nygc_papa_as_ale.enriched.spliced_reads.{.y}.svg"),
               path = "processed/nygc/enriched_jncs/ale/svg/",
               plot = .x,
               device = svg,
               height = 14,
               width = 14,
               units = "in",
               dpi = "retina"
      ),
      .progress = T
)


#



# # create name as junction coords - gene_name
# selective_jnc <- path_summ_sel %>% pull(paste_into_igv_junction) %>%
#   set_names(paste(path_summ_sel$paste_into_igv_junction, path_summ_sel$gene_name, sep = "_") %>% str_replace_all(":|-", "_"))
# 
# enriched_jnc <- path_summ_enr %>% pull(paste_into_igv_junction) %>%
#   set_names(paste(path_summ_enr$paste_into_igv_junction, path_summ_enr$gene_name, sep = "_") %>% str_replace_all(":|-", "_"))
# 
# enriched_rpm <- ale_enriched_by_rpm %>%
#   pull(paste_into_igv_junction) %>%
#   set_names(paste(ale_enriched_by_rpm$paste_into_igv_junction, ale_enriched_by_rpm$gene_name, sep = "_") %>% str_replace_all(":|-", "_"))
# 




# path_summ_ale <- read_tsv("processed/nygc/expression_by_pathology_ale_all.tsv")
# 
# # Classify as selective or enriched according to AL's criteria
# path_summ_ale <- path_summ_ale %>%
#   mutate(selective = fraction_not_path <= 0.005 & fraction_path >= 0.01,
#          enriched = !selective & (fraction_path > fraction_not_path * 2.5))
# 
# count(path_summ_ale, selective, enriched)
# # A tibble: 3 × 3
# # selective enriched     n
# # <lgl>     <lgl>    <int>
# #   1 FALSE     FALSE       48
# # 2 FALSE     TRUE        60
# # 3 TRUE      FALSE        7
# 
# 
# path_summ_sel <- filter(path_summ_ale, selective)
# path_summ_enr <- filter(path_summ_ale, enriched)
# 
# 
# 
# # create name as junction coords - gene_name
# selective_jnc <- path_summ_sel %>% pull(paste_into_igv_junction) %>%
#   set_names(paste(path_summ_sel$paste_into_igv_junction, path_summ_sel$gene_name, sep = "_") %>% str_replace_all(":|-", "_"))
# 
# enriched_jnc <- path_summ_enr %>% pull(paste_into_igv_junction) %>%
#   set_names(paste(path_summ_enr$paste_into_igv_junction, path_summ_enr$gene_name, sep = "_") %>% str_replace_all(":|-", "_"))
# 
# enriched_rpm <- ale_enriched_by_rpm %>%
#   pull(paste_into_igv_junction) %>%
#   set_names(paste(ale_enriched_by_rpm$paste_into_igv_junction, ale_enriched_by_rpm$gene_name, sep = "_") %>% str_replace_all(":|-", "_"))
# 
# sel_jnc_plots <- map(selective_jnc,
#                      ~ plot_junction(.x))
# 
# enr_jnc_plots <- map(enriched_jnc,
#                      ~ plot_junction(.x))
# 
# enr_rpm_jnc_plots <- map(enriched_rpm,
#                          ~ plot_junction(.x))
# 
# dir.create("processed/nygc/selective_jncs/", recursive = T)
# dir.create("processed/nygc/enriched_jncs/", recursive = T)
# dir.create("processed/nygc/enriched_rpm_jncs/", recursive = T)
# #
# walk2(sel_jnc_plots,
#       names(sel_jnc_plots),
#       ~ ggsave(filename = glue::glue("processed/nygc/selective_jncs/2023-09-13_nygc_papa_as_ale.selective.spliced_reads.{.y}.png"),
#                plot = .x,
#                height = 14,
#                width = 14,
#                units = "in",
#                dpi = "retina"
#       )
# )
# 
# walk2(enr_jnc_plots,
#       names(enr_jnc_plots),
#       ~ ggsave(filename = glue::glue("processed/nygc/enriched_jncs/2023-09-13_nygc_papa_as_ale.enriched.spliced_reads.{.y}.png"),
#                plot = .x,
#                height = 14,
#                width = 14,
#                units = "in",
#                dpi = "retina"
#       )
# )
# 
# walk2(enr_rpm_jnc_plots,
#       names(enr_rpm_jnc_plots),
#       ~ ggsave(filename = glue::glue("processed/nygc/enriched_rpm_jncs/2023-09-13_nygc_papa_as_ale.enriched_rpm.spliced_reads.{.y}.png"),
#                plot = .x,
#                height = 14,
#                width = 14,
#                units = "in",
#                dpi = "retina"
#       )
# )
# 
# # Make bar plot of %/fraction of path tissues detected
# # do separately for specific and enriched
# # Also a plot of combined
# 
# 
# # just selective, ordered by most tissues obs
# path_summ_sel %>%
#   mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
#          plot_name = fct_reorder(plot_name, fraction_path)) %>%
#   ggplot(aes(x = fraction_path * 100, y = plot_name)) +
#   geom_col() +
#   theme_bw(base_size = 14) +
#   labs(title = "Selective AS-ALE junctions",
#        x = "% of TDP-43 pathological tissues cryptic detected",
#        y = "Junction|Gene name")
# 
# # just selective, include fraction non_path
# path_summ_sel %>%
#   mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
#          plot_name = fct_reorder(plot_name, fraction_path)) %>%
#   pivot_longer(starts_with("fraction"), names_to = "tissue_group", values_to = "fraction") %>%
#   ggplot(aes(x = fraction * 100, y = plot_name, fill = tissue_group)) +
#   geom_col(position = "dodge") +
#   scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
#   theme_bw(base_size = 14) +
#   labs(title = "Selective AS-ALE junctions",
#        x = "% of TDP-43 pathological tissues cryptic detected",
#        y = "Junction|Gene name")
# 
# 
# # just enriched, include frac non-path
# # (potentially a bit misleading if expressed in a lot of tissues, could be an expression change rather than strong cryp expression - check PSI)
# path_summ_enr %>%
#   mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
#          plot_name = fct_reorder(plot_name, fraction_path)) %>%
#   pivot_longer(starts_with("fraction"), names_to = "tissue_group", values_to = "fraction") %>%
#   ggplot(aes(x = fraction * 100, y = plot_name, fill = tissue_group)) +
#   geom_col(position = "dodge") +
#   scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
#   theme_bw(base_size = 14) +
#   labs(title = "Enriched AS-ALE junctions",
#        subtitle = "Enriched = fraction path > 2.5* fraction non-path",
#        x = "% of TDP-43 pathological tissues cryptic detected",
#        y = "Junction|Gene name")
# 
# bind_rows(Enriched = path_summ_enr, Selective = path_summ_sel, .id = "sel_group") %>%
#   mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
#          plot_name = fct_reorder(plot_name, fraction_path),
#          sel_group = factor(sel_group, levels = c("Selective", "Enriched"))) %>%
#   pivot_longer(starts_with("fraction"), names_to = "tissue_group", values_to = "fraction") %>%
#   ggplot(aes(x = fraction * 100, y = plot_name, fill = tissue_group)) +
#   facet_wrap("~ sel_group", ncol = 1, scales = "free_y") +
#   geom_col(position = "dodge") +
#   scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
#   theme_bw(base_size = 10) +
#   labs(title = "Selective & Enriched AS-ALE junctions",
#        x = "% of TDP-43 pathological tissues cryptic detected",
#        y = "Junction|Gene name") +
#   theme(axis.text.y = element_text(size = rel(0.65)))
# 
# 
# #enriched by rpm
# ale_enriched_by_rpm %>%
#   filter(!paste_into_igv_junction %in% path_summ_sel$paste_into_igv_junction) %>%
#   mutate(plot_name = paste(paste_into_igv_junction, gene_name, sep = "|"),
#          plot_name = fct_reorder(plot_name, fraction_path)) %>%
#   pivot_longer(starts_with("fraction"), names_to = "tissue_group", values_to = "fraction") %>%
#   ggplot(aes(x = fraction * 100, y = plot_name, fill = tissue_group)) +
#   geom_col(position = "dodge") +
#   scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
#   theme_bw(base_size = 14) +
#   labs(title = "Enriched AS-ALE junctions",
#        subtitle = "Enriched by wilcoxon test of spliced RPM",
#        x = "% of TDP-43 pathological tissues cryptic detected",
#        y = "Junction|Gene name")