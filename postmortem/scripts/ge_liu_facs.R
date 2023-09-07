library(tidyverse)


calc_means <- function(df, sample_cols, out_col) {
  
  df[, out_col] <- rowMeans(df[, sample_cols])
  df

}


# gene expression matrices
liu_abundance <- read_tsv("processed/gene_exprn/2023-09-07_liu_facs_salmon_summarised_abundance.tsv")
liu_counts <- read_tsv("processed/gene_exprn/2023-09-07_liu_facs_salmon_summarised_counts.tsv")
liu_counts_norm <- read_tsv("processed/gene_exprn/2023-09-07_liu_facs_salmon_summarised_counts_deseq_normalised.tsv")

# cryptic-containing genes
cryptics_summary <- read_tsv("data/cryptics_summary_all_events.tsv")
cryptic_gn <- distinct(cryptics_summary, gene_name) %>% pull()


# all cols minus gene_id and gene_name
sample_cols <- colnames(liu_abundance)[3:ncol(liu_abundance)]
# all TDP-negative samples
neg_sample_cols <- colnames(liu_abundance)[str_detect(colnames(liu_abundance), "_negative$")]
# all TDP-positive samples
pos_sample_cols <- colnames(liu_abundance)[str_detect(colnames(liu_abundance), "_positive$")]

# calculate condition-wise & overall means for all three categories
liu_abundance <- calc_means(liu_abundance, sample_cols, "mean_all") 
liu_abundance <- calc_means(liu_abundance, neg_sample_cols, "mean_neg") 
liu_abundance <- calc_means(liu_abundance, pos_sample_cols, "mean_pos") 

liu_counts <- calc_means(liu_counts, sample_cols, "mean_all") 
liu_counts <- calc_means(liu_counts, neg_sample_cols, "mean_neg") 
liu_counts <- calc_means(liu_counts, pos_sample_cols, "mean_pos") 

liu_counts_norm <- calc_means(liu_counts_norm, sample_cols, "mean_all") 
liu_counts_norm <- calc_means(liu_counts_norm, neg_sample_cols, "mean_neg") 
liu_counts_norm <- calc_means(liu_counts_norm, pos_sample_cols, "mean_pos") 

# filter for cryptics
liu_abundance_cryp <- filter(liu_abundance, gene_name %in% cryptic_gn)
liu_counts_cryp <- filter(liu_counts, gene_name %in% cryptic_gn)
liu_counts_norm_cryp <- filter(liu_counts_norm, gene_name %in% cryptic_gn)

# combine into single df of means
liu_all_means <-  bind_rows(counts_norm = select(liu_counts_norm_cryp, gene_id, gene_name, starts_with("mean_")),
          counts = select(liu_counts_cryp, gene_id, gene_name, starts_with("mean_")),
          tpm = select(liu_abundance_cryp, gene_id, gene_name, starts_with("mean_")),
          .id = "abundance_type")


# standard counts filter - try a few different mean cutoffs (how many excluded?)

seq(0,100, 10) %>%
  set_names() %>%
  map(~ liu_all_means %>%
        filter(abundance_type == "counts_norm",
               mean_all <= .x) %>%
        summarise(n_genes = n_distinct(gene_name))
      )  %>%
  bind_rows(.id = "max_counts")

# max_counts n_genes
# <chr>        <int>
#   1 0                1
# 2 10              11
# 3 20              14
# 4 30              15
# 5 40              20
# 6 50              22
# 7 60              28
# 8 70              28
# 9 80              31
# 10 90              36
# 11 100             39

# ~ 40 / 283 have GENE-level counts < 100. I.e. unique fragments across the whole gene, so could be very few reads over individual isoforms/exons. Hard to say what's enough coverage to say difficult to detect an individual exon
# abvoe is conventional filter for gene-level stats analysis, not necessarily event-level
# 11 genes have GE below conventional filter

# Try some plots to visualise
liu_all_means %>%
  filter(abundance_type == "counts_norm") %>%
  mutate(gene_name = fct_reorder(gene_name, mean_all)) %>%
  ggplot(aes(x = log2(mean_all + 1), y = gene_name)) +
  geom_point(size = rel(0.1)) +
  geom_vline(xintercept = log2(11), linetype = "dashed") +
  geom_vline(xintercept = log2(26), linetype = "dashed") +
  geom_vline(xintercept = log2(50), linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(title = "Cryptic-containing gene expression") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
  
liu_all_means %>%
    filter(abundance_type == "counts_norm") %>%
    mutate(gene_name = fct_reorder(gene_name, mean_all)) %>%
    ggplot(aes(x = log2(mean_all + 1))) +
  stat_ecdf() +
    geom_vline(xintercept = log2(11), linetype = "dashed") +
    geom_vline(xintercept = log2(26), linetype = "dashed") +
    geom_vline(xintercept = log2(51), linetype = "dashed") +
    theme_bw(base_size = 14) +
    labs(title = "Gene expression of cryptic-containing genes")
  
  
# Random thoughts:
# - Are 'enriched' genes over-represented in highly expressed genes? See where they pop up in curve, possibly histograms of expression
# - could make similar cumulative distribution plot (maybe just individual dots) coloured by event type as well

