library(tidyverse)
library(writexl)

datasets <- read_tsv("data/2023-11-22_paper_tdp43_collection_library_statistics.tsv")
# df containing cleaned cooridnate columns extracted from quant GTF
le_id_coords <- read_tsv("processed/le_id_collapsed_coords.quant.last_exons.tsv")
# yes/no binding within plotting windows at representative coordinates/boundaries
iclip_summary <- read_tsv(file.path(outdir, "2024-11-26_cryptic_iclip_summary_cleaned_combined.tsv"))
outdir <- "processed"

# replace unpublished with this study
datasets <- mutate(datasets, paper_doi = if_else(paper_doi == "unpublished",
                                                 "this-study",
                                                 paper_doi)
                   )


# replace Klim doi with URL for consistency - https://doi.org/10.1038/s41593-018-0300-4
datasets <- mutate(datasets, paper_doi = if_else(paper_doi == "10.1038/s41593-018-0300-4",
                                                 "https://doi.org/10.1038/s41593-018-0300-4",
                                                 paper_doi))


# selective ALE junctions - include all searched but highlight selective at top
nygc_ale_all <- read_tsv("data/expression_by_pathology_ale_all.tsv")
nygc_ale_all


# Since enriched not really referred to in text (and a rough metric), instead sort remaining by fold enrichment and remove column
nygc_ale_all_cleaned <- nygc_ale_all %>%
  mutate(fold_enrichment = (fraction_path + 0.001) / (fraction_not_path + 0.001)) %>%
  arrange(desc(selective), desc(fold_enrichment)) %>%
  select(-enriched) %>%
  rename(count_not_path = not_path, count_path = path) %>%
  relocate(starts_with("count_"), .after = gene_name)


# ribo-seq for APA-containing genes
riboseq_all <- read_tsv("data/2023-09-26_riboseq_volcano_plot_df.tsv")

# filter for ALEs only
riboseq_ale <- filter(riboseq_all, !is.na(le_id))

riboseq_ale_clean <- riboseq_ale %>%
  select(gene_id, gene_name, le_id, event_type = plot_event_type, baseMean, log2FoldChange, lfcSE, log2FoldChangeShrink, lfcSEShrink,
         pvalue, padj) %>%
  distinct(.keep_all = T) %>%
  arrange(pvalue)

# some genes are duplicated (multiple le_ids) - collapse le_ids and report single row per APA gene
riboseq_ale_le_id_clpsd <- riboseq_ale_clean %>%
  group_by(gene_id, gene_name) %>%
  summarise(le_id_clpsd = paste(le_id, collapse = ";"))

riboseq_ale_clean <- riboseq_ale_clean %>%
  left_join(riboseq_ale_le_id_clpsd, by = c("gene_id", "gene_name")) %>%
  distinct(gene_id, gene_name, le_id_clpsd, .keep_all = T) %>%
  select(-le_id, le_id = le_id_clpsd) %>%
  relocate(le_id, .after = gene_name)


# RNA-seq for APA-containing genes
## TODO: output from script

# Cryptic ALE event summary df
cryptics_summary <- read_tsv("data/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv")


unique(cryptics_summary$simple_event_type)
# [1] "spliced"                 "distal_3utr_extension"   "bleedthrough"            "bleedthrough,spliced"    "proximal_3utr_extension"

# update event type label for consistency with paper
cryptics_summary_clean <- cryptics_summary %>%
  mutate(event_type = case_when(simple_event_type == "spliced" ~ "ALE",
                   simple_event_type == "bleedthrough" ~ "IPA",
                   simple_event_type == "distal_3utr_extension" ~ "3'Ext",
                   simple_event_type == "proximal_3utr_extension" ~ "3'Shortening",
                   TRUE ~ "Complex")) %>%
  select(-simple_event_type) %>%
  relocate(event_type, .after = gene_name) 

# replace commas with semi colons in experiment_name
cryptics_summary_clean <- mutate(cryptics_summary_clean, experiment_name = str_replace_all(experiment_name, ",", ";"))

# replace skndz with sknbe2 to reflect new genotyping
cryptics_summary_clean <- mutate(cryptics_summary_clean, experiment_name = str_replace_all(experiment_name, "skndz", "sknbe2"))


# Unsure how coordinates selected, so use representative coordinates from last exon BEDs (iCLIP, motifs) to annotate
rep_coords_paths <- c("ALE" = "data/2023-12-14_papa_cryptic_spliced.background_shsy5y.last_exons.bed",
  "IPA" = "data/2023-12-14_papa_cryptic_bleedthrough_uniq.background_shsy5y.last_exons.bed",
  "3'Ext" = "data/2023-12-15_papa_cryptic_d3utr.background_shsy5y.last_exons.distal.bed",
  "3'Shortening" = "data/2023-12-15_papa_cryptic_d3utr_proximal.background_shsy5y.last_exons.proximal.bed")

#
rep_coords_df <-  rep_coords_paths %>%
  map(~ read_tsv(.x, col_names = c("chromosome", "start", "end", "name", "score", "strand"), show_col_types = F) %>%
        # extract info out of name column (e.g. ENSG00000120942.14_3|UBIAD1|distal|background)
        separate(name, into = c("le_id", "gene_name", NA, "regn_group"), sep = "\\|") %>%
        # get cryptics only and tidy up
        filter(regn_group == "cryptic") %>%
        select(-all_of(c("regn_group", "score", "gene_name"))) %>%
        mutate(across(all_of(c("start", "end")), as.character))
      ) %>%
  bind_rows()

# won't be represented in these BEDs - keep original coordinates
complex_coords <- cryptics_summary_clean %>%
  filter(event_type == "Complex") %>%
  select(le_id, chromosome, start, end, strand) %>%
  mutate(across(all_of(c("start", "end")), as.character))
  

# join into summary df and replace coordinates
cryptics_summary_clean <- select(cryptics_summary_clean, -all_of(c("chromosome", "start", "end", "strand"))) %>%
  left_join(bind_rows(rep_coords_df, complex_coords), by = "le_id")

# Complex coords not used for iCLIP so do not assign start/end coords
# Add in all coordinates from quant GTF and differentiate coords
cryptics_summary_clean <- cryptics_summary_clean %>%
  mutate(start = if_else(event_type == "Complex", NA, start),
         end = if_else(event_type == "Complex", NA, end)) %>%
  left_join(select(le_id_coords, le_id, start, end),
            by = "le_id", suffix = c("_rep", "_all")) %>%
  relocate(strand, .after = everything())
  
# Add in binding information for each event

cryptics_summary_clean <- cryptics_summary_clean %>%
  left_join(select(iclip_summary, le_id, binding_start, binding_end), by = c("le_id"))

cryptics_summary_clean

# summmary df containing non-cryptic significant hits from the analysis
noncryptics_summary_clean <- read_tsv(file.path(outdir, "2023-05-24_i3_cortical_zanovello.all_datasets.non_cryptic_sig_apas.summary.tsv"))

# NYGC metadata summary
nygc_metadata <- read_tsv("../postmortem/processed/nygc/2024-11-27_nygc_metadata_summary.tdp_subtypes.tsv")

# HeLa target genes
hela_chipseq_targets <- read_tsv("data/2024-01-09_ferguson_hela_chipseq_target_gene_lists.tsv")

# table of primer sequences
primers <- read_tsv("data/Oligonucleotides_table_final.tsv")

# table mapping cell type/model to the experiment(s) in which they were used
cell2exper <- read_csv("data/celltype-to-experiment-mapping.csv")
cell2exper <- drop_na(cell2exper)

suppl_list <- list("Supplementary_Table_1" = datasets,
                   "Supplementary_Table_2" = cryptics_summary_clean,
                   "Supplementary_Table_3" = noncryptics_summary_clean,
                   "Supplementary_Table_4" = nygc_metadata,
                   "Supplementary_Table_5" = nygc_ale_all_cleaned,
                   "Supplementary_Table_6" = riboseq_ale_clean,
                   "Supplementary_Table_7" = hela_chipseq_targets,
                   "Supplementary_Table_8" = cell2exper,
                   "Supplementary_Table_9" = primers)

if (!dir.exists(outdir)) {dir.create(outdir)}

write_xlsx(suppl_list,
           path = file.path(outdir, "2024-11-28_supplementary_tables.xlsx"), col_names = T, format_headers = T)

walk2(.x = suppl_list,
      .y = names(suppl_list),
      ~ write_tsv(.x, 
                  file = file.path(outdir, paste("2024-11-28_", .y, ".tsv", sep = "")),
                  col_names = T,
                  na = ""
                  )
      )
