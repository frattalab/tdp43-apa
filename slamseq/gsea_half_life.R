library(tidyverse)
library(fgsea)
# load in function to get ranked gene list
source("../riboseq/helpers.R")
set.seed(123)


slamseq_hl <- read_tsv("processed/2023-08-22_i3cortical_slamseq_grandr_kinetics.tsv")

# i3 cortical cryptic ale containing genes
cryptic_ale_gene_lists <- read_rds("../riboseq/processed/gsea_riboseq.ale_gene_lists.rds")


# rank experiment-wide absolute FCs and by up/downregulated
slamseq_hl <- slamseq_hl %>%
  mutate(hl_rank_desc_all = rank(desc(abs(log2FoldHalfLife))),
         hl_num_all = n(),
         frac_rank_all = hl_rank_desc_all / hl_num_all,
         dirn = sign(log2FoldHalfLife)) %>%
  group_by(dirn) %>%
  mutate(hl_rank_desc_dirn = rank(desc(abs(log2FoldHalfLife))),
         hl_num_dirn = n(),
         frac_rank_dirn = hl_rank_desc_dirn / hl_num_dirn)


# for each cryptic event, annotate with sign/direction of fold change
# convert lists to a df
cryptic_ale_gene_df <- map2(cryptic_ale_gene_lists,
     names(cryptic_ale_gene_lists),
     ~ enframe(.x, name = "rown", value = "gene_id") %>% select(-rown)) %>% 
  bind_rows(.id = "event_type")

cryptic_ale_slamseq_hl <- cryptic_ale_gene_df %>%
  left_join(slamseq_hl, by = c("gene_id" = "Gene"))

write_tsv(cryptic_ale_slamseq_hl, "processed/2023-10-19_cryptic_ale_grandr_kinetics_hl_ranks.tsv", col_names = T)

# sorted vector of log2fold half lives
hl_fc_ranks <- get_ranked_gene_list(slamseq_hl, score_col = "log2FoldHalfLife", name_col = "Gene")




# run standard GSEA
gsea_slamseq_cryp <- fgsea(pathways = cryptic_ale_gene_lists,
                           stats = hl_fc_ranks)

# collapse leading edge genes to comma separated string, then join to results df

# first generate id2name so can have insert readable gene symbols
id2name <- distinct(slamseq_hl, gene_id = Gene, gene_name)

pway_le <- gsea_slamseq_cryp %>%
  unnest_longer(leadingEdge) %>%
  left_join(id2name, by = c("leadingEdge" = "gene_id")) %>%
  group_by(pathway) %>%
  summarise(leadingEdgeName = paste(gene_name, collapse = ";"),
            leadingEdgeID = paste(leadingEdge, collapse = ";"),
  )

# join collapsed leading edge to df
gsea_slamseq_cryp_df <- gsea_slamseq_cryp %>%
  select(-leadingEdge) %>%
  left_join(pway_le, by = "pathway")

write_tsv(gsea_slamseq_cryp_df, "processed/2023-10-19_slamseq_gsea_ale_types_results.tsv", col_names = T)
saveRDS(hl_fc_ranks, "processed/gsea_slamseq.fc_ranks.rds")