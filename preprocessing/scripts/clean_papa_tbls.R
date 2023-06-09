library(tidyverse)

# Script to re-process PAPA results tables and add additional output
# 1. Combined df containing results for all experiments (alogn with per-experiment)
# 2. De-duplicating gene name values (due to a stupid design decision on my part)
# 3. Re-annotating bleedthrough events if called as novel extension of annotated bleedthrough event

# Script assumes you have different comparisons under different directories. 
# The directory structure should also match that of a PAPA run e.g.
# <main_output_dir>/differential_apa/dexseq_apa.results.processed.tsv
# The script grabs the experiment_name from 'main_output_dir', so make sure it is named appropriately


### Functions

#' convert dataset specific PPAU columns to standardised base/treatment columns
generalise_ppau_cols <- function(df, delta_col_prefix = "delta_PPAU", mean_col_prefix = "mean_PPAU") {
  
  delta_col <- colnames(df)[str_detect(colnames(df), paste("^", delta_col_prefix, sep = ""))]
  # base condition column always comes leftmost, followed by treatment/contrast column
  base_col <- colnames(df)[str_detect(colnames(df), paste("^", mean_col_prefix, sep = ""))][1]
  treat_col <- colnames(df)[str_detect(colnames(df), paste("^", mean_col_prefix, sep = ""))][2]
  
  df %>%
    dplyr::rename(delta_PPAU_treatment_control = !!sym(delta_col),
                  mean_PPAU_base = !!sym(base_col),
                  mean_PPAU_treatment = !!sym(treat_col))
}


#' generalise the exon usage coefficient columns to base and contrast for a given comparison
generalise_coeff_cols <- function(dapa_df, contrast_sep = "vs", out_prefix = "UsageCoefficient") {
  
  # contrast name always in format of KD vs CTL
  contrast_name <- unique(dapa_df$contrast_name)
  stopifnot(length(contrast_name) == 1)
  
  keys <- str_split_1(contrast_name, contrast_sep)
  cont_key <- keys[1]
  base_key <- keys[2]
  
  out_base <- paste(out_prefix, "base", sep = "_")
  out_cont <- paste(out_prefix, "contrast", sep = "_")
  
  rename(dapa_df,
         !!out_base := !!sym(base_key),
         !!out_cont := !!sym(cont_key)
  )
  
}


#' Standardise DEXSeq's contrast-specific log2FoldChange column
generalise_log2fold_col <- function(dapa_df, contrast_sep = "vs", in_prefix = "log2fold", out_prefix = "log2fold") {
  
  # contrast name always in format of KD vs CTL
  contrast_name <- unique(dapa_df$contrast_name)
  stopifnot(length(contrast_name) == 1)
  
  keys <- str_split_1(contrast_name, contrast_sep)
  cont_key <- keys[1]
  base_key <- keys[2]
  
  exp_col <- paste(in_prefix, cont_key, base_key, sep = "_")
  
  stopifnot(exp_col %in% colnames(dapa_df))
  
  out_col <- paste(out_prefix, "contrast", "base", sep = "_")
  
  rename(dapa_df, !!out_col := !!sym(exp_col))
  
}

#' helper function to collapse duplicated values in a delimited string to non-redundant values
#' I wasn't very smart with gene_name column in PAPA output - often get duplicated gene name values...
collapse_names <- function(col, split=",") {
  
  apply(str_split(col,  split, simplify = T), # generate a matrix of strings separated by split 
        MARGIN = 1, # loop over rows
        FUN = function(x) {paste(unique(x)[unique(x) != ""],
                                 collapse = split)
        }, # collapse split gene names to unique values, then recombine if there are multiple IDs remaining
        simplify = T)
}


## code


# find dexseq results tables recursively
dapa_paths <- list.files(path = "data/gencode_v40_reruns/i3_cortical_upf1_zanovello_novel",
                         pattern = "dexseq_apa.results.processed.tsv",
                         full.names = T,
                         recursive = T,
                         include.dirs = T)

# set name of each path to the 'main_output_dir' defined by PAPA
dapa_paths <- set_names(dapa_paths, basename(dirname(dirname(dapa_paths))))

# Read in and combine dfs, standardising column names so can be combined without column duplication
dapa_comb <- dapa_paths %>%
  map(~ read_tsv(.x, show_col_types = F)) %>%
  # make sure PPAU columns have standardised names
  map(~ generalise_ppau_cols(.x)) %>%
  # now standardise dexseq results columns
  map(~ generalise_log2fold_col(.x) %>%
        generalise_coeff_cols()
  ) %>%
  bind_rows(.id = "experiment_name")


# Sometimes gene names have duplicate values separated by commas (WHY SAM)
# Will drop to unique values and recombine if necessary
dapa_comb <- mutate(dapa_comb, gene_name = collapse_names(gene_name))

dapa_comb

# Simple event type annotation (spliced, bleedthrough, distal_3utr_extension) - properly identifying bleedthrough events
# Few bleedthroughs get called as novel last_exon_extensions of annotated bleedthrough last exons 
# Re-classify based on relative position in gene. If last_exon_extension event_types: 
# le_id has expressed further downstream last exons = bleedthrough
# otherwise = distal_3utr_extension

# internal_exon_extension = bleedthrough
# any _spliced = spliced

dapa_comb <- dapa_comb %>%
  group_by(experiment_name, gene_id) %>%
  # le_ids are annotated with numeric suffix to denote position of le in gene (1 = most proximal). Also strand aware
  arrange(le_id, .by_group = T) %>%
  mutate(posn = row_number(),
         simple_event_type = case_when(event_type == "last_exon_extension" & posn == max(posn) ~ "distal_3utr_extension",
                                       event_type == "last_exon_extension" & posn != max(posn) ~ "bleedthrough",
                                       str_detect(event_type, "spliced$") ~ "spliced",
                                       event_type %in% c("first_exon_extension", "internal_exon_extension") ~ "bleedthrough"
         )
  ) %>%
  ungroup() %>%
  # remove intermediate columns
  select(-posn) %>%
  # sort columns as in input
  arrange(experiment_name, gene.qvalue)

# dapa_comb


if (!dir.exists("processed/PAPA/")) { dir.create("processed/PAPA/", recursive = T) }

# write combined TSV
write_tsv(dapa_comb,
          file = "processed/PAPA/2023-05-24_i3_cortical_zanovello.all_datasets.dexseq_apa.results.processed.cleaned.tsv",
          col_names = T )


# Write per experiment TSVs
for (experiment in unique(dapa_comb$experiment_name)) {
  
  dapa_comb %>%
    filter(experiment_name == experiment) %>%
    write_tsv(.,
              file = paste("processed/PAPA/2023-05-24_i3_cortical_zanovello.",
                           experiment,
                           ".dexseq_apa.results.processed.cleaned.tsv",
                           sep = ""),
              col_names = T
    )
}
