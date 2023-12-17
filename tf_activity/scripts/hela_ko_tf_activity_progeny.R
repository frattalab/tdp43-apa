library(tidyverse)
library(decoupleR)
library(tidytext)
source("scripts/helpers.R")

# What is the overlap between ELK1/ELK4 ChIP-seq targets and known signalling pathways?
# Is inferred activity a consequence of signalling pathway induction upon TDP-43 KD, or due specifically to ELK1/ELK4?


progeny_top100 <- read_tsv("data/2023-12-16_progeny_top100_homosapiens.tsv")
progeny_top500 <- read_tsv("data/2023-12-16_progeny_top500_homosapiens.tsv")
progeny_all <- read_tsv("data/2023-12-16_progeny_full_homosapiens.tsv.gz")
ferguson_deseq <- read_csv("data/tdp43_kd_collection/deseq2_outputs/ferguson.hela.DESEQ2_results.csv")
# remove genes with NA padj value, keep one row per gene name
ferguson_deseq <- clean_deseq_df(ferguson_deseq)


# Select the strongest associated pathway for each gene

# select by min pvalue
progeny_all_p <- progeny_all %>%
  group_by(gene) %>%
  slice_min(p.value) %>%
  ungroup() %>%
  rename(target = gene, source = pathway)

# select by highest absolute weight
progeny_all_w <- progeny_all %>%
  group_by(gene) %>%
  slice_max(abs(weight)) %>%
  ungroup() %>%
  rename(target = gene, source = pathway)

#' quick and dirty function to load target lists from Rdata with a bunch of objects
load_target_list <- function(path) {
  
  load(path)
  return(chipseq_hela_target_lists)
  
}

# ELK1 + ELK4 ChIP-seq target lists
chipseq_hela_target_lists <- load_target_list("processed/ferguson_hela/2023-11-29_hela_ko_tf_activity_gsea.Rdata")

# combine into a df of gene lsit + target genes
chipseq_hela_target_df <- map(chipseq_hela_target_lists, ~
      enframe(.x, name = NULL, value = "target")) %>%
  bind_rows(.id = "gene_set")

# Assign targets to known signalling pathways

progeny_list <- list("progeny_top100" = progeny_top100,
     "progeny_top500" = progeny_top500,
     "progeny_all_p" = progeny_all_p,
     "progeny_all_w" = progeny_all_w) 

chipseq_hela_progeny_list <- progeny_list %>%
  map(~ left_join(chipseq_hela_target_df, .x, by = "target")
  )

# count number of overlaps for each signalling pathway and target gene set
chipseq_hela_progeny_counts <- map(chipseq_hela_progeny_list,
    ~ count(.x, gene_set, source, sort = T, .drop = FALSE))
  
# count the number of genes in each signalling pathway
progeny_counts <- map(progeny_list,
                      ~ count(.x, source, sort = T))

# combine the two and calculate fraction of genes in different signalling pathways
chipseq_hela_progeny_counts <- map2(.x = chipseq_hela_progeny_counts,
      .y = progeny_counts,
      ~ left_join(.x, .y, by = "source", suffix = c("", ".all")) %>%
        mutate(frac.all = n / n.all))

chipseq_hela_progeny_counts_df <- chipseq_hela_progeny_counts %>%
  bind_rows(.id = "progeny_set")

chipseq_hela_progeny_counts_df %>%
  mutate(source = reorder_within(source, frac.all, progeny_set),
         plot_gene_set = str_remove_all(gene_set, "^chipseq_hela_")) %>%
  ggplot(aes(x = source, y = frac.all, fill = plot_gene_set)) +
  facet_wrap("~ progeny_set", scales = "free") +
  scale_x_reordered() + geom_col(position = "dodge") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "top") +
  guides(fill = guide_legend(nrow = 2))


## main ones to watch out for are MAPK and EGFR

## Run progeny to infer activities

#
ferguson_splicing <- read_csv("data/tdp43_kd_collection/splicing/ferguson_collected_splicing.csv")

# any genes with 'significant' splicing change
ferguson_diff_spliced <- ferguson_splicing %>%
  filter(probability_changing > 0.95) %>%
  distinct(gene_name) %>%
  pull()


# get matrix of deseq results
ferguson_deseq_mtx <- ferguson_deseq %>%
  select(gene_name, log2FoldChange, stat) %>%
  drop_na(gene_name) %>%
  column_to_rownames(var = "gene_name") %>%
  as.matrix()

# repeat, but remove differentially spliced genes + TDP itself
ferguson_nospl_deseq_mtx <- ferguson_deseq %>%
  filter(!gene_name %in% c("TARDBP", ferguson_diff_spliced)) %>%
  select(gene_name, log2FoldChange, stat) %>%
  drop_na(gene_name) %>%
  column_to_rownames(var = "gene_name") %>%
  as.matrix()

# run progeny with diff gene lists
ferguson_deseq_progeny <- map(progeny_list,
    ~ run_mlm(mat=ferguson_deseq_mtx[, 'stat', drop=FALSE],
        net=.x,
        .source='source',
        .target='target',
        .mor='weight',
        minsize = 5)
    )

ferguson_deseq_progeny_nospl <-  map(progeny_list,
                                     ~ run_mlm(mat=ferguson_nospl_deseq_mtx[, 'stat', drop=FALSE],
                                         net=.x,
                                         .source='source',
                                         .target='target',
                                         .mor='weight',
                                         minsize = 5)
                                     )

# plot the activity scores for different progeny subsets
map2(.x = ferguson_deseq_progeny_nospl,
     .y = names(ferguson_deseq_progeny_nospl),
     ~ .x %>%
       mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
       rstatix::add_significance(p.col = "p_adj", output.col = "p_signif") %>%
      ggplot(aes(x = reorder(source, score), y = score, label = p_signif)) + 
  geom_bar(aes(fill = score), stat = "identity") +
    geom_text(nudge_y = -1) +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(title = .y,
       x = "Pathways")
  )


# TODO: for top 100 & top 500, remove target genes that are linked to significant pathways
# Remove/mask those genes from target lists & rerun GSEA - are associations still significant