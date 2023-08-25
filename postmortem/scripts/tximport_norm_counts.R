library(tximport)
library(DESeq2)
library(tidyverse)



salmon_quant_dir <- "../data/liu_facs/2023-08-23_salmon_reference_filtered_cryptic_le_combined/salmon_quant/"
tx2gene2name <- read_tsv("../processed/2023-08-23_reference_filtered_cryptic_le_combined.tx2gene2name.tsv")

# tximport needs transcript_id | gene_id in that order
tx2gene <- tx2gene2name %>%
  select(transcript_id, gene_id)

# get a named vector of paths to quant.sf files
# name of final directory = sample name
quant_files <- list.files(path = "../data/liu_facs/2023-08-23_salmon_reference_filtered_cryptic_le_combined/salmon_quant",
           pattern = "quant.sf", 
           full.names = T, 
           recursive = T) %>%
  set_names(basename(dirname(.)))

# run tximport - 1st to get matrix for summarised counts
liu_facs_gene_counts <- tximport(files = quant_files,
                                 type = "salmon",
                                 txOut = FALSE, # summarise to gene
                                 countsFromAbundance = "lengthScaledTPM",
                                 tx2gene = tx2gene)


# for TPMs and counts, return gene_ids to column & add gene name information
liu_facs_dfs <-  c("abundance", "counts") %>%
  set_names() %>%
  map(~ liu_facs_gene_counts[[.x]] %>%
        as_tibble(rownames = "gene_id") %>%
        left_join(distinct(select(tx2gene2name, gene_id, gene_name)),
                  by = "gene_id") %>%
        relocate(gene_name, .after = gene_id)
        )

if (!dir.exists("../processed/gene_exprn/")) {dir.create("../processed/gene_exprn", recursive = T)}

write_tsv(liu_facs_dfs$abundance, "../processed/gene_exprn/2023-08-25_liu_facs_salmon_summarised_abundance.tsv", col_names = T)
write_tsv(liu_facs_dfs$counts, "../processed/gene_exprn/2023-08-25_liu_facs_salmon_summarised_counts.tsv", col_names = T)
