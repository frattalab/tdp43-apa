library(tidyverse)

# Script to get a table of non-cryptic APAs produced by the pipeline 

# Combined differential analysis
dexseq_all <- read_tsv("data/2023-05-24_i3_cortical_zanovello.all_datasets.dexseq_apa.results.processed.cleaned.tsv.gz")
# df containing cleaned cooridnate columns
le_id_coords <- read_tsv("processed/le_id_collapsed_coords.quant.last_exons.tsv")
outdir <- "processed"


dexseq_sig <- filter(dexseq_all, padj < 0.05)

# keep all but the intermediate curve datasets
exp_to_keep <- unique(dexseq_sig$experiment_name)[str_detect(unique(dexseq_sig$experiment_name), "_curve_",negate = T) | 
                                            unique(dexseq_sig$experiment_name) %in% c("zanovello_skndz_curve_1", "zanovello_shsy5y_curve_0075")]

dexseq_sig <- filter(dexseq_sig, experiment_name %in% exp_to_keep)

# exclude cryptics from the output
dexseq_sig_noncryp <- dexseq_sig %>%
  mutate(cryptic = padj < 0.05 & mean_PPAU_base < 0.1 & delta_PPAU_treatment_control > 0.1) %>%
  group_by(experiment_name, gene_id) %>%
  # annotate gene as cryptic-containing per experiment
  mutate(cryptic_gene = if_else(any(cryptic), T, F)) %>%
  ungroup() %>%
  # now, never report a gene if ever cryptic
  group_by(gene_id) %>%
  filter(!any(cryptic_gene)) %>%
  ungroup() %>%
  select(-cryptic, -cryptic_gene)
  

# get number of unique experiments, the name, median deltas and control usages, annotation status
summ_dexseq_sig_noncryp <- dexseq_sig_noncryp %>%
  group_by(le_id) %>%
  summarise(gene_name = paste(unique(gene_name), collapse = ","), # note that order is not guaranteed for le_ids with multiple gene_names
            simple_event_type = paste(sort(unique(simple_event_type)), collapse = ","),
            n_exper_sig = n_distinct(experiment_name),
            experiment_name = paste(sort(unique(experiment_name)), collapse = ","),
            annot_status = paste(sort(unique(annot_status)), collapse = ","),
            chromosome = unique(chromosome),
            start = paste(sort(unique(start)), collapse = ","),
            end = paste(sort(unique(end)), collapse = ","),
            strand = unique(strand)
            ) %>%
  arrange(desc(n_exper_sig))

dexseq_sig_noncryp 

# summarise effect sizes (for sorting)
dexseq_sig_eff_sizes <- dexseq_sig_noncryp %>%
  distinct(le_id, experiment_name, mean_PPAU_base, mean_PPAU_treatment, delta_PPAU_treatment_control) %>%
  group_by(le_id) %>%
  summarise(median_delta = median(delta_PPAU_treatment_control)
  )


# Collapse to minimal statistics per le_id and experiment
stat_cols <- c("exonBaseMean", "pvalue", "padj", "mean_PPAU_base", "mean_PPAU_treatment", "delta_PPAU_treatment_control")
dexseq_sig_noncryp_stats <- dexseq_sig_noncryp %>%
  select(le_id, experiment_name, all_of(stat_cols)) %>%
  distinct(.keep_all = T)


# Combine with annotation columns
dexseq_sig_noncryp_out <- dexseq_sig_noncryp_stats %>%
  left_join(select(summ_dexseq_sig_noncryp, -experiment_name), by = "le_id") %>%
  left_join(dexseq_sig_eff_sizes, by = "le_id") %>%
  arrange(desc(n_exper_sig),desc(abs(median_delta)), experiment_name, le_id)
  
# Re-annotate event types 

# Re-annotated 'spliced' as 3'Shortening events i.e. are the proximal 3'UTR of a locus with a novel 3'UTR extension

# vector of le_ids that are novel 3'UTR extensions
d3utr_le_ids <- dexseq_all %>%
  filter(simple_event_type == "distal_3utr_extension") %>%
  distinct(le_id) %>%
  pull()

# construct a 'putative le_id' for sig events - i.e. if it had a distal 3'UTR, this would be it's le_id
# gene_id prefix is same, but number suffix = le_id + 1
dexseq_sig_noncryp_out <- dexseq_sig_noncryp_out %>%
  mutate(distal_le_id = paste(str_remove_all(le_id, "_\\d+$"), # remove suffixed number
                              as.character(as.integer(str_split_i(le_id, "_", 2)) + 1),
                              sep = "_"))


# now check if putative distal_le_id is in fact a novel 3'UTR extension
# if so, update the simple event type annotation to proximal_3utr_extension
dexseq_sig_noncryp_out <- dexseq_sig_noncryp_out %>%
  mutate(has_distal_3utr = simple_event_type == "spliced" & distal_le_id %in% d3utr_le_ids) %>%
  mutate(simple_event_type = if_else(has_distal_3utr,
                                     "proximal_3utr_extension", simple_event_type)) %>%
  # remove intermediate columns
  select(-distal_le_id, -has_distal_3utr)

# Update event type labels for consistency with paper
dexseq_sig_noncryp_out <- dexseq_sig_noncryp_out %>%
  mutate(event_type = case_when(simple_event_type == "spliced" ~ "ALE",
                              simple_event_type == "bleedthrough" ~ "IPA",
                              simple_event_type == "distal_3utr_extension" ~ "3'Ext",
                              simple_event_type == "proximal_3utr_extension" ~ "3'Shortening",
                              TRUE ~ "Complex")) %>%
  select(-simple_event_type) %>%
  relocate(all_of(c("gene_name", "event_type")), .after = le_id) 


dexseq_sig_noncryp_out

# TODO: clean up coordinates
dexseq_sig_noncryp_out <- dexseq_sig_noncryp_out %>%
  select(-start, -end) %>%
  left_join(le_id_coords, by = "le_id") %>%
  relocate(all_of(c("start", "end")), .after = chromosome)
  
# check that no coord cols have weird characters (i.e. x +e power y)
dexseq_sig_noncryp_out %>%
  filter(str_detect(start, "e") | str_detect(end, "e")) %>%
  nrow()
# [1] 0

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}
write_tsv(dexseq_sig_noncryp_out, file.path(outdir, "2023-05-24_i3_cortical_zanovello.all_datasets.non_cryptic_sig_apas.summary.tsv"))

