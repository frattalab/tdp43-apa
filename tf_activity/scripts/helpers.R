

#' Tidy up DESeq2 results table for analyses of enrichment
#' Removes low count genes (those filtered out by DESeq's independent filtering) & ensures 1 row per gene 
clean_deseq_df <- function(df, padj_col = "padj", id_col = "gene_name") {
  
  # Remove genes with NA padj - i.e. too low counts for stat testing &/ filtered out due to low counts (to maximize n sig)
  df %>%
    drop_na(!!sym(padj_col)) %>%
    distinct(!!sym(id_col), .keep_all = T)
  
}

#' Calculate a signed pvalue to rank genes for 'extent of up/downregulation' 
#' Sign of log2FoldChange * (-log10(adjusted pvalue))
#' This gives large values to most sig p-values & separates up/downreg genes
#' If gene has a p-value of 0, log10 transformation will give infinite value
#' In case of ties, genes with p-vals of 0 will be ranked by foldchange (smallest first) and will have rank added/subtracted to max signed pvalue (if up/downreg)
add_signed_pval <- function(df, out_col = "signed_pvalue", fc_col = "log2FoldChange", pval_col = "pvalue") {
  
  # Want to be able to rank from high +ve to low -ve - large pvals in the middle
  # -log10(pval) exaggerates small pvalues, so most sig have large values
  
  df_sign <- df %>%
    mutate("{out_col}" := sign(!!sym(fc_col)) * -log10(!!sym(pval_col))) %>%
    arrange(!!sym(out_col))
  
  # Get infinite signed_padj if pvalue = 0
  # for these, rank by FC (smallest to highest), assigning a value greater than max signed pvalue
  
  # Get max signed p-value   
  max_up <- filter(df_sign,
                   !is.infinite(!!sym(out_col))
  ) %>%
    pull(!!sym(out_col)) %>%
    max()
  
  # Min signed p-value
  max_down <- filter(df_sign,
                     !is.infinite(!!sym(out_col))
  ) %>%
    pull(!!sym(out_col)) %>%
    min()
  
  # Assign value for infinite signed pval above max/min value, with largest FC assigned largest value
  df_inf_adj <- df_sign %>%
    filter(is.infinite(!!sym(out_col))) %>%
    mutate(dirn = if_else(!!sym(fc_col) > 0, "up", "down")) %>%
    group_by(dirn) %>%
    # smallest FC placed first in group (assigned 1 to add/subract from fc rank)
    # largest FC will have largest value added/subtracted
    arrange(abs(!!sym(fc_col)), .by_group = T) %>%
    mutate(fc_rank = row_number()) %>%
    ungroup() %>%
    mutate("{out_col}":= if_else(dirn == "up",
                                 max_up + fc_rank,
                                 max_down - fc_rank)) %>%
    select(- fc_rank, -dirn)
  
  # Add infinite adjusted values back to df
  df_sign <- df_sign %>%
    filter(!is.infinite(!!sym(out_col))) %>%
    bind_rows(., df_inf_adj) %>%
    arrange(!!sym(out_col))
  
  return(df_sign)
}


#' Calculate a value of log2FoldChange scaled by the -log10(transformed p-value) to rank for 'extent of up/downregulation' 
#' Sign of log2FoldChange * (-log10(adjusted pvalue))
#' This gives large values to most sig p-values & largest fold changes
#' If gene has a p-value of 0, log10 transformation will give infinite value
#' In case of ties, genes with p-vals of 0 will be ranked by foldchange (smallest first) and will have rank added/subtracted to max signed pvalue (if up/downreg)
#' since purpose here is to rank, the actual values do not matter so much at the top end
add_log2fold_pval <- function(df, out_col = "pvalScaledLog2FoldChange", fc_col = "log2FoldChange", pval_col = "pvalue") {
  
  # Want to be able to rank from high +ve to low -ve - large pvals in the middle
  # -log10(pval) exaggerates small pvalues, so most sig have large values
  
  df_sign <- df %>%
    mutate("{out_col}" := !!sym(fc_col) * -log10(!!sym(pval_col))) %>%
    arrange(!!sym(out_col))
  
  # Get infinite signed_padj if pvalue = 0
  # for these, rank by FC (smallest to highest), assigning a value greater than max signed pvalue
  
  # Get max signed p-value   
  max_up <- filter(df_sign,
                   !is.infinite(!!sym(out_col))
  ) %>%
    pull(!!sym(out_col)) %>%
    max()
  
  # Min signed p-value
  max_down <- filter(df_sign,
                     !is.infinite(!!sym(out_col))
  ) %>%
    pull(!!sym(out_col)) %>%
    min()
  
  # Assign value for infinite signed pval above max/min value, with largest FC assigned largest value
  df_inf_adj <- df_sign %>%
    filter(is.infinite(!!sym(out_col))) %>%
    mutate(dirn = if_else(!!sym(fc_col) > 0, "up", "down")) %>%
    group_by(dirn) %>%
    # smallest FC placed first in group (assigned 1 to add/subract from fc rank)
    # largest FC will have largest value added/subtracted
    arrange(abs(!!sym(fc_col)), .by_group = T) %>%
    mutate(fc_rank = row_number()) %>%
    ungroup() %>%
    mutate("{out_col}":= if_else(dirn == "up",
                                 max_up + fc_rank,
                                 max_down - fc_rank)) %>%
    select(- fc_rank, -dirn)
  
  # Add infinite adjusted values back to df
  df_sign <- df_sign %>%
    filter(!is.infinite(!!sym(out_col))) %>%
    bind_rows(., df_inf_adj) %>%
    arrange(!!sym(out_col))
  
  return(df_sign)
}


#' Get a named vector of genes ranked by score from smallest to largest value
get_ranked_gene_list <- function(df, score_col = "signed_pvalue", name_col = "gene_name") {
  
  # Make sure smallest first (how fgsea wants it ordered)
  df <- df %>% arrange(!!sym(score_col))
  
  set_names(df[[score_col]], df[[name_col]])
  
}

#' Convert dorothea df to 
dorothea_to_gsea <- function(df, split_by_mor = T) {
  
  #
  if (split_by_mor) {
    
    group_cols <- c("tf", "mor")
    
  } else {
    group_cols <- c("tf")
  }
  
  # first group by tf & direction of regulation for each conf level
  df_grpd <- group_by(df,across(all_of(group_cols)))
  
  # generate a 'group name' - combine group column values with '_' separator

  if (split_by_mor) {
    # set to tf + direction of regulation
    grp_names <- mutate(group_keys(df_grpd), # df of tf & mor (in same order as split list of groups)
                      tf_mor = paste(tf, mor, sep = "_")) %>%
    pull(tf_mor)
    
  } else { 
    # just set to tf
    grp_names <- pull(group_keys(df_grpd)) # returns single col df
  }

  # generate list of tfs and target genes
  df_grpd %>%
    # split into list of dfs (by group)
    group_split() %>%
    # set the names
    set_names(grp_names) %>%
    map(~ unique(pull(.x, target)))

}
