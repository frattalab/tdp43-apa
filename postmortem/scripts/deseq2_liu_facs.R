library(tidyverse)
library(DESeq2)

sample_tbl <- read_csv("data/liu_facs/liu_facs_sample_sheet_dexseq_meta.csv")
counts <- read_tsv("processed/gene_exprn/2023-09-07_liu_facs_salmon_summarised_counts.tsv")
base_key="TDPpos"
contrast_key="TDPneg"
contrast_name="TDPneg_vs_TDPpos"
covs <- c("gender", "patient")

# convert condition col to factor, setting base_key as reference level
sample_tbl <- mutate(sample_tbl, condition = factor(condition, levels = c(base_key, contrast_key)))
sample_tbl <- mutate(sample_tbl, across(all_of(covs), as.factor))

# make column names in counts matrix and sample table align

# columns = only sample names
counts_mtx <- column_to_rownames(select(counts, -gene_name), "gene_id")
# subset for sampels in sample table
counts_mtx <- counts_mtx[, colnames(counts_mtx) %in% sample_tbl$sample_name]

# check the column oreder matches the row order in sample table
all(sample_tbl$sample_name == colnames(counts_mtx))
# TRUE


# construct deseq object
dds <- DESeqDataSetFromMatrix(countData = round(counts_mtx),
                              colData = sample_tbl,
                              design = ~ patient + condition)

# low level filter for min counts in min number of samples
smallestGroupSize <- 7 # number of samples in one condition
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# run deseq
dds <- DESeq(dds,parallel = T)

# 
res <- results(dds, name = "condition_TDPneg_vs_TDPpos") %>%
  as_tibble(rownames = "gene_id") %>%
  left_join(distinct(counts, gene_id, gene_name)) %>%
  arrange(padj)

# shrink the FCs...
fc_shrink <- lfcShrink(dds, coef="condition_TDPneg_vs_TDPpos", type="apeglm")

# add to results df
res <- fc_shrink %>% 
  as_tibble(rownames = "gene_id") %>%
  dplyr::rename(log2FoldChangeShrink = log2FoldChange,
                lfcSEShrink = lfcSE
                ) %>%
  dplyr::select(gene_id, log2FoldChangeShrink, lfcSEShrink) %>%
  left_join(res, .)


# transformation whilst blind to the design formula
vsd_blind <- vst(dds, blind=TRUE)

# transformation whilst blind to the design formula
vsd_notblind <- vst(dds, blind=FALSE)

plotPCA(vsd_blind, intgroup=c("condition", "patient","gender"))

plotPCA(vsd_blind, intgroup=c("condition", "patient","gender"), returnData = TRUE) %>%
  ggplot(aes(x = PC1, y = PC2, shape = gender, colour = condition, label = patient)) +
  geom_point(size = 4) +
  ggrepel::geom_label_repel()

# do full transofrmation in awareness of design matrix
vsd_full_notblind <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "parametric")

###
if (!dir.exists("processed/gene_exprn/deseq2/")) {dir.create("processed/gene_exprn/deseq2/", recursive = T)}
write_tsv(res, "processed/gene_exprn/deseq2/2023-20-11_deseq2_liu_facs_results.tsv")
assay(vsd_full_notblind) %>%
  as_tibble(rownames = "gene_id") %>%
  left_join(distinct(counts, gene_id, gene_name)) %>%
  write_tsv("processed/gene_exprn/deseq2/2023-20-11_deseq2_liu_facs_vst_notblind_counts_mtx.tsv")
