library(tidyverse)

#' Get a named vector of genes ranked by score from smallest to largest value
get_ranked_gene_list <- function(df, score_col = "signed_padj", name_col = "gene_name") {
  
  # Make sure smallest first (how fgsea wants it ordered)
  df <- df %>% arrange(!!sym(score_col))
  
  set_names(df[[score_col]], df[[name_col]])
  
}