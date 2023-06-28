set.seed(123)
library(DESeq2)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)
library(tibble)
library(ggrepel)
library(apeglm)

cds_counts <- list.files("data/gene_cds_counts", pattern = "_results.txt$",full.names = T) %>%
  set_names(str_remove(basename(.), "_featureCounts_results.txt")) %>%
  map(~ read_tsv(.,comment = "#", show_col_types = F)) %>%
  # counts colname = path to provided BAM file
  map(~ rename(., counts = ends_with(".bam"))) %>%
  bind_rows(.id = "sample_id")

slice_sample(cds_counts, n = 100)

# pivot to wide count matrix
counts_mtx <- cds_counts %>%
  select(sample_id, gene_id = Geneid, counts) %>%
  pivot_wider(names_from = sample_id, values_from = counts) %>%
  column_to_rownames("gene_id")

counts_mtx


# generate a normalised count matrix (for later output)
# named vector of size factors for each sample in matrix
sfs <-DESeq2::estimateSizeFactorsForMatrix(counts_mtx)

# Divide each sample's counts by corresponding size factor
norm_counts_mtx <- sweep(counts_mtx, MARGIN = 2, STATS = sfs, FUN = '/')

norm_counts_mtx

# add gene names back to normed count matrix
norm_counts_mtx_gn <- norm_counts_mtx %>%
  rownames_to_column("Geneid") %>%
  left_join(distinct(select(cds_counts, Geneid, gene_name), Geneid, .keep_all=T), by = "Geneid")


### RUN DESeq2 analysis

# make a metadata table - rownames = sample_names, condition = CTL/KD 
coldata <- colnames(counts_mtx) %>% 
  as_tibble_col("sample_name") %>%
  mutate(condition = as.factor(if_else(str_detect(sample_name,"tdp"),
                                       "KD", "CTL")
  )
  ) %>%
  column_to_rownames("sample_name")

# coldata
# doublecheck sample names
all(rownames(coldata) %in% colnames(counts_mtx))
all(rownames(coldata) == colnames(counts_mtx))


dds <- DESeqDataSetFromMatrix(countData = counts_mtx,
                              colData = coldata,
                              design = ~ condition)

# pre-filter genes with low counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set denominator to CTL just to be safe
dds$condition <- relevel(dds$condition, ref = "CTL")

dds <- DESeq(dds,parallel = T)

# results(dds)

deseq_res_df_shrink <- lfcShrink(dds,coef = "condition_KD_vs_CTL", type = "apeglm")
#deseq_res_df_shrink


# convert to df and add gene name info
deseq_res_df <- results(dds) %>%
  as_tibble(rownames = "gene_id") %>%
  left_join(select(norm_counts_mtx_gn, Geneid, gene_name),
            by = c("gene_id" = "Geneid"))

# add lfcShrink and lfcSE to df 
deseq_res_df <- deseq_res_df_shrink %>%
  as_tibble(rownames = "gene_id") %>%
  select(gene_id,
         log2FoldChangeShrink = log2FoldChange,
         lfcSEShrink = lfcSE) %>%
  left_join(deseq_res_df, ., by = "gene_id") %>%
  relocate(ends_with("Shrink"), .before = stat)

# deseq_res_df

if (!dir.exists("processed")) {dir.create("processed", recursive = T)}

write_tsv(deseq_res_df, "processed/2023-05-08_i3_cortical_riboseq_deseq2_results.tsv",col_names = T)
write_tsv(norm_counts_mtx_gn, "processed/2023-05-08_i3_cortical_riboseq_sf_normed_count_matrix.tsv", col_names = T)
