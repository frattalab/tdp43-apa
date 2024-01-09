library(tidyverse)

## Script to output selected target gene lists to file (supplementary data)

load("processed/ferguson_hela/2023-11-29_hela_ko_tf_activity_gsea.Rdata")

names(elk1_target_list)

# All ELK1, unique ELk1, All ELK4, unique ELK4
lists_to_keep <- c("chipseq_hela_both_elk1", "chipseq_hela_unique_elk1", "chipseq_hela_unique_elk4", "chipseq_hela_all_elk4")

hela_chipseq_target_lists <- elk1_target_list[names(elk1_target_list) %in% lists_to_keep]

# append differentially spliced
hela_chipseq_target_lists[["majiq_differentially_spliced"]] <- ferguson_diff_spliced


# Want a df where column name = target list name, row values = gene_name

names(hela_chipseq_target_lists)
# [1] "chipseq_hela_both_elk1"       "chipseq_hela_unique_elk1"     "chipseq_hela_unique_elk4"     "chipseq_hela_all_elk4"        "majiq_differentially_spliced"

names(hela_chipseq_target_lists) <- c("ELK1_all", "ELK1_unique", "ELK4_unique", "ELK4_all", "majiq_differentially_spliced")

# create 2 col df of row_number + name of gene set
hela_chipseq_target_lists_dfs <- map2(.x = hela_chipseq_target_lists,
                                      .y = names(hela_chipseq_target_lists),
    ~ enframe(.x, name = "row_number", value = .y)) 


# find length of longest gene set
max_len <- hela_chipseq_target_lists %>%
  map(length) %>%
  reduce(max)

# create a combined df of gene sets, putting each as a column
hela_chipseq_targets_df <- c(list(gene_number = enframe(1:max_len, name = NULL, value = "row_number")),
  hela_chipseq_target_lists_dfs) %>%
  reduce(left_join, by = "row_number") %>%
  select(-row_number)
  
write_tsv(hela_chipseq_targets_df, "processed/ferguson_hela/2024-01-09_ferguson_hela_chipseq_target_gene_lists.tsv",col_names = T, na = "")

