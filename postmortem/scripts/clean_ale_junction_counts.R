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


#' Add nygc metadata to an SJ count bed
#' dt - data.table of BED file of SJs, where score column (5th) is number of spliced reads. must also have 'sample' and 'paste_into_igv_junction' cols
#' metadata - patient metadata to join to junctions by 'sample' column
#' id_cols - additional cols other than paste_into_igv_junction to retain as information for the junction
bed_add_metadata <- function(dt, metadata, id_cols) {
  
  info_cols <- c("paste_into_igv_junction", id_cols)
  
  # get a dataframe of unique SJs searched (+ associated gene IDs)
  junction_info = distinct(dt, !!!syms(info_cols)) %>% as.data.table()
  # junction_info = unique(raw_spliced_ale[,.(paste_into_igv_junction,gene_id,gene_name)])
  
  # create a table of all combinations of junctions, sample_ids and counts (filling missing counts as 0 if applicable)
  temp = complete(dt[,.(paste_into_igv_junction, sample, V5)],
                  paste_into_igv_junction, sample,fill = list(V5 = 0)) |> 
    as.data.table()
  
  # rename counts column
  setnames(temp,"V5","spliced_reads")
  # Add in sample ID & junction annotation information (gene id etc)
  temp = temp |> left_join(junction_info)
  # add tissue, disease classification & other metadata 
  temp = temp |> left_join(metadata)
  
  # add columns specifying whether tissue with expected TDP path (disease_tissue) & whether has TDP path (tdp_path)
  create_formated_metadata(temp) |> unique()
}

summarise_path <- function(dt, id_cols, observed_minreads = 2,
                           selective_minfrac_path = 0.01,
                           selective_maxfrac_nonpath = 0.005, enriched_minratio = 2.5) {
  
  info_cols <- c("paste_into_igv_junction", id_cols)
  
  # get a dataframe of unique SJs searched (+ associated gene IDs)
  junction_info = distinct(dt, !!!syms(info_cols)) %>% as.data.table()
  
  summ_df <- dt |> 
    dplyr::select(disease_tissue,spliced_reads,tdp_path,sample,paste_into_igv_junction, all_of(id_cols)) |> 
    unique() |> 
    filter(disease_tissue == TRUE) |> 
    mutate(observed = spliced_reads >= observed_minreads) |> 
    group_by(tdp_path,paste_into_igv_junction) |> 
    summarise(n_obs = sum(observed)) |> 
    ungroup() |> 
    pivot_wider(values_from = 'n_obs',
                names_from = 'tdp_path') |> 
    left_join(junction_info)  |> 
    as.data.table() |> 
    mutate(fraction_not_path = not_path / total_notpath) |> 
    mutate(fraction_path = path / total_path)
  
  # add whether passes selective/enriched criteria
  summ_df %>%
  mutate(selective = fraction_path >= selective_minfrac_path & fraction_not_path <= selective_maxfrac_nonpath,
         enriched = !selective & fraction_path > (fraction_not_path * enriched_minratio)) %>%
    arrange(desc(selective), desc(enriched))
  
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


library_depth = arrow::read_parquet('processed/nygc/nygc_library_sizes.parquet') |> 
    dplyr::rename(sample = sample_id)
meta_data_full = fread("data/nygc/NYGC_all_RNA_samples_support.tsv")

meta_data = meta_data_full[,.SD,.SDcols = c("sample", "individual", "region", "tissue", "tissue_clean",
                                       "disease", "disease_full", "age",
                                       "onset", "mutations", "pathology")]


# read in BED file of ALE SJ counts (output of bedops_parse_star_junctions)
raw_spliced_ale = fread('data/nygc/2023-09-11_papa_as_ale_cryptics.aggregated.clean.annotated.bed')
raw_spliced_ale = raw_spliced_ale |> 
    mutate(sample = str_extract(V4,"CGND-HRA-\\d{5}"))

# extract metadata from 'name' field of ALE BED file
raw_spliced_ale[, le_id := tstrsplit(V7,"\\|")[[1]]]
raw_spliced_ale[, gene_name := tstrsplit(V7,"\\|")[[2]]]
raw_spliced_ale[,paste_into_igv_junction := paste0(V1,":",V2,"-",V3)]

spliced_counts_ale <- bed_add_metadata(raw_spliced_ale, meta_data, c("le_id", "gene_name"))

# number of path tissue samples
total_path = spliced_counts_ale |> 
    group_by(tdp_path) |> 
    summarize(n_samp = n_distinct(sample)) |> 
    filter(tdp_path == "path") |> pull()

# number of non-path tissue samples
total_notpath = spliced_counts_ale |> 
    group_by(tdp_path) |> 
    summarize(n_samp = n_distinct(sample)) |> 
    filter(tdp_path == "not_path") |> pull()

# for disease tissues, identify whether junctions are observed (>= 2 spliced reads)
# and calculate fraction of tissues where observed (relative to own group)
# Also add labels whether passes path specificity/enrichment criteria
expression_by_pathology_ale <- summarise_path(spliced_counts_ale, c("le_id", "gene_name")) 

counts_enriched_selective_ale <- count(expression_by_pathology_ale, selective, enriched)

write_tsv(expression_by_pathology_ale, "processed/nygc/expression_by_pathology_ale_all.tsv")
write_tsv(filter(expression_by_pathology_ale, selective), "processed/nygc/expression_by_pathology_ale_selective.tsv")
write_tsv(filter(expression_by_pathology_ale, enriched), "processed/nygc/expression_by_pathology_ale_enrichedfrac.tsv")
write_tsv(counts_enriched_selective_ale, "processed/nygc/expression_by_pathology_ale_counts.tsv")

## Process Seddighi df

raw_spliced_seddighi <- fread("data/nygc/2023-09-19_seddighi_all_cryptics.aggregated.clean.annotated.bed")

# add junction ID
raw_spliced_seddighi[,paste_into_igv_junction := paste0(V1,":",V2,"-",V3)]

# extract sample from Name field of BED file
raw_spliced_seddighi = raw_spliced_seddighi |> 
  mutate(sample = str_extract(V4,"CGND-HRA-\\d{5}"))

# extract ID columns
raw_spliced_seddighi <- raw_spliced_seddighi %>%
  rename(Symbol = V7) %>%
  mutate(gene_name = str_split_i(Symbol, "_", 1)) %>%
  as.data.table()

spliced_counts_seddighi <- bed_add_metadata(raw_spliced_seddighi, meta_data, c("Symbol", "gene_name"))

# number of path tissue samples
total_path_seddighi = spliced_counts_seddighi |> 
  group_by(tdp_path) |> 
  summarize(n_samp = n_distinct(sample)) |> 
  filter(tdp_path == "path") |> pull()

# number of non-path tissue samples
total_notpath_seddighi = spliced_counts_seddighi |> 
  group_by(tdp_path) |> 
  summarize(n_samp = n_distinct(sample)) |> 
  filter(tdp_path == "not_path") |> pull()

# summarise 
expression_by_pathology_seddighi <- summarise_path(spliced_counts_seddighi, c("Symbol", "gene_name"))
counts_enriched_selective_seddighi <- count(expression_by_pathology_seddighi, selective, enriched)

# write to file
write_tsv(expression_by_pathology_seddighi, "processed/nygc/expression_by_pathology_seddighi_all.tsv")
write_tsv(filter(expression_by_pathology_seddighi, selective), "processed/nygc/expression_by_pathology_seddighi_selective.tsv")
write_tsv(filter(expression_by_pathology_seddighi, enriched), "processed/nygc/expression_by_pathology_seddighi_enrichedfrac.tsv")
write_tsv(counts_enriched_selective_seddighi, "processed/nygc/expression_by_pathology_seddighi_counts.tsv")

# another way to check enrichment -----------------------------------------

# obvserved_enough = spliced_counts_ale[spliced_reads >=2] %>% 
#     group_by(paste_into_igv_junction) %>% 
#     mutate(n_obs = n_distinct(sample)) %>% 
#     filter(n_obs > 10) %>% 
#     pull(paste_into_igv_junction) |> unique()
# 
# long_rpm = spliced_counts_ale |> 
#     filter(paste_into_igv_junction %in% obvserved_enough) |> 
#     left_join(library_depth) |> 
#     mutate(rpm = 10^6 * (spliced_reads / library_size)) |> 
#     filter(disease_tissue == TRUE) |> 
#     mutate(predict_path = ifelse(disease %in% c('ALS-TDP','FTD-TDP'),
#                                  "TDP-path",
#                                  "non-TDP path")) 
# 
# # spliced reads per million - sig different between TDP-path and non-TDP path?
# wilcox_rpm = long_rpm |> 
#     group_by(paste_into_igv_junction) |> 
#     nest() %>%
#     mutate(lm_obj = map(data, ~wilcox.test(rpm~predict_path,data = .))) %>%
#     mutate(lm_tidy = map(lm_obj, broom::tidy)) %>%
#     ungroup() |> 
#     unnest(cols = c(lm_tidy)) |> 
#     select(-data) |>
#   mutate(p.value.adj = p.adjust(p.value, method = "BH")) %>%
#   arrange(p.value.adj)
#   
# 
# sig_wil = wilcox_rpm |> filter(p.value < 0.01) |> pull(paste_into_igv_junction)
# 
# ale_enriched_by_rpm = expression_by_pathology_ale[paste_into_igv_junction %in% sig_wil]


