library(tidyverse)
library(data.table)
library(arrow)
library(here)

create_formated_metadata = function(dt){
    dt[,disease := ifelse(disease == "FTD" & pathology %in% c("FTD-TAU","FTD-FUS"),"FTD-non-TDP",disease)]
    dt[,disease := ifelse(disease == "FTD","FTD-TDP",disease)]
    
    dt[,disease := ifelse(grepl("ALS",disease), "ALS-TDP",disease)]
    dt[,disease := ifelse(grepl("ALS",disease) & mutations %in% c("FUS","SOD1"),"ALS-non-TDP",disease)]
    
    
    dt = dt |>  mutate(disease_tissue = case_when(grepl("FTD",disease_full) & grepl("Front|Temp",tissue_clean)  ~ T,
                                                  (grepl("ALS",disease) & grepl("Cord|Motor|Front|Temp",tissue_clean))  ~ T,
                                                  (grepl("Occipital|Sensory",tissue_clean)) ~ F,
                                                  (grepl("Control",disease) & grepl("Cord|Cortex",tissue_clean)) ~ T,
                                                  TRUE ~ F)) 
    dt = dt |> 
        mutate(tdp_path = case_when((disease %in% c("ALS-TDP","FTD-TDP") & disease_tissue == T) ~ 'path',
                                    T ~ "not_path"))
    dt = unique(dt[disease != "Other"])
    return(dt)
}


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

#
path_summ <- read_csv("processed/nygc/expression_by_pathology_ale (1).csv")


library_depth = arrow::read_parquet('processed/nygc/nygc_library_sizes.parquet') |> 
    dplyr::rename(sample = sample_id)
meta_data_full = fread("data/nygc/NYGC_all_RNA_samples_support.tsv")

meta_data = meta_data_full[,.SD,.SDcols = c("sample", "individual", "region", "tissue", "tissue_clean",
                                       "disease", "disease_full", "age",
                                       "onset", "mutations", "pathology")]


raw_spliced_ale = fread('data/nygc/2023-09-11_papa_as_ale_cryptics.aggregated.clean.annotated.bed')
raw_spliced_ale = raw_spliced_ale |> 
    mutate(sample = str_extract(V4,"CGND-HRA-\\d{5}"))

raw_spliced_ale[, gene_id := tstrsplit(V7,"\\|")[[1]]]
raw_spliced_ale[, gene_name := tstrsplit(V7,"\\|")[[2]]]
raw_spliced_ale[,paste_into_igv_junction := paste0(V1,":",V2,"-",V3)]
# 
junction_info = unique(raw_spliced_ale[,.(paste_into_igv_junction,gene_id,gene_name)])
# setnames(junction_info, "V6","strand")
temp = complete(raw_spliced_ale[,.(paste_into_igv_junction,sample,V5)],
         paste_into_igv_junction, sample,fill = list(V5 = 0)) |> as.data.table()
# 
setnames(temp,"V5","spliced_reads")
# 
temp = temp |> left_join(junction_info)
# 
temp = temp |> left_join(meta_data)
spliced_counts_ale = create_formated_metadata(temp) |> unique()
rm(temp)

total_path = spliced_counts_ale |> 
    group_by(tdp_path) |> 
    summarize(n_samp = n_distinct(sample)) |> 
    filter(tdp_path == "path") |> pull()

total_notpath = spliced_counts_ale |> 
    group_by(tdp_path) |> 
    summarize(n_samp = n_distinct(sample)) |> 
    filter(tdp_path == "not_path") |> pull()

expression_by_pathology_ale = spliced_counts_ale |> 
    dplyr::select(disease_tissue,spliced_reads,tdp_path,sample,paste_into_igv_junction,gene_id,gene_name) |> 
    unique() |> 
    filter(disease_tissue == TRUE) |> 
    mutate(observed = spliced_reads >= 2) |> 
    group_by(tdp_path,paste_into_igv_junction) |> 
    summarise(n_obs = sum(observed)) |> 
    ungroup() |> 
    pivot_wider(values_from = 'n_obs',
                names_from = 'tdp_path') |> 
    left_join(junction_info)  |> 
    as.data.table() |> 
    mutate(fraction_not_path = not_path / total_notpath) |> 
    mutate(fraction_path = path / total_path)



ale_selective = expression_by_pathology_ale |>
    filter(fraction_not_path <= 0.005) |>
    filter(fraction_path >= 0.01)

ale_enriched = expression_by_pathology_ale[fraction_path > fraction_not_path * 2.5]


# another way to check enrichment -----------------------------------------

obvserved_enough = spliced_counts_ale[spliced_reads >=2] %>% 
    group_by(paste_into_igv_junction) %>% 
    mutate(n_obs = n_distinct(sample)) %>% 
    filter(n_obs > 10) %>% 
    pull(paste_into_igv_junction) |> unique()

long_rpm = spliced_counts_ale |> 
    filter(paste_into_igv_junction %in% obvserved_enough) |> 
    left_join(library_depth) |> 
    mutate(rpm = 10^6 * (spliced_reads / library_size)) |> 
    filter(disease_tissue == TRUE) |> 
    mutate(predict_path = ifelse(disease %in% c('ALS-TDP','FTD-TDP'),
                                 "TDP-path",
                                 "non-TDP path")) 

wilcox_rpm = long_rpm |> 
    group_by(paste_into_igv_junction) |> 
    nest() %>%
    mutate(lm_obj = map(data, ~wilcox.test(rpm~predict_path,data = .))) %>%
    mutate(lm_tidy = map(lm_obj, broom::tidy)) %>%
    ungroup() |> 
    unnest(cols = c(lm_tidy)) |> 
    select(-data)

sig_wil = wilcox_rpm |> filter(p.value < 0.01) |> pull(paste_into_igv_junction)

ale_enriched_by_rpm = expression_by_pathology_ale[paste_into_igv_junction %in% sig_wil]



### Some quick plots

# Classify as selective or enriched according to AL's criteria
path_summ <- path_summ %>%
  mutate(selective = fraction_not_path <= 0.005 & fraction_path >= 0.01,
         enriched = !selective & (fraction_path > fraction_not_path * 2.5))

count(path_summ, selective, enriched)
# A tibble: 3 × 3
# selective enriched     n
# <lgl>     <lgl>    <int>
#   1 FALSE     FALSE       48
# 2 FALSE     TRUE        60
# 3 TRUE      FALSE        7


path_summ_sel <- filter(path_summ, selective)
path_summ_enr <- filter(path_summ, enriched)



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
