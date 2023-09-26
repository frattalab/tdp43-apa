library(tidyverse)
library(fgsea)
source("helpers.R")
set.seed(123)




### GSEA on FC of different event types
# i.e. are the subset of ALE event type genes that show consistent directionality significant with respect to random gene sets of the same size?

plot_df_cryp_volc <- read_tsv("processed/2023-09-26_riboseq_volcano_plot_df.tsv")
deseq_res_df <- read_tsv("processed/2023-05-08_i3_cortical_riboseq_deseq2_results.tsv")


# GSEA on FC ranks
# make gene list of different event types
# NB: use gene names to avoid erroneous duplicate gene name ("HSPA14" annotated to ENSG00000284024.2 & ENSG00000187522.16)
tmp_grp <- plot_df_cryp_volc %>%
  filter(!simple_event_type %in% c("other", "complex")) %>%
  group_by(simple_event_type)

cryptic_ale_gene_lists <- tmp_grp %>%
  group_split() %>%
  set_names(group_keys(tmp_grp) %>% pull()) %>% # group keys returns 1 col df in case of single group
  map(~ unique(.x$gene_id))


# get a named vector of genes ranked by ribo-seq fold change 
riboseq_fc_ranks <- get_ranked_gene_list(deseq_res_df, score_col = "log2FoldChangeShrink", name_col = "gene_id")

# check for duplicate names
dup_genes <- names(riboseq_fc_ranks)[duplicated(names(riboseq_fc_ranks))]
length(dup_genes)

# check for duplicated values (fold changes)
dup_fcs <- duplicated(riboseq_fc_ranks)
which(dup_fcs) / length(riboseq_fc_ranks)
# [1] 0.3479498 0.3763135 0.4342261 0.4490796 0.4491586 0.4651181 0.4722288 0.4747571 0.4802876 0.4828158 0.4907166
# [12] 0.5097574 0.5154460 0.5200284 0.5226357 0.5372521 0.5634037 0.5733586 0.5893182 0.5997472 0.6183140 0.6194991
# [23] 0.6210002 0.6612151

# most are near middle of gene list... unlikely to affect enrichments too greatly


# run standard GSEA
gsea_riboseq_cryp <- fgsea(pathways = cryptic_ale_gene_lists,
                           stats = riboseq_fc_ranks)

# collapse leading edge genes to comma separated string, then join to results df

# first generate id2name so can have insert readable gene symbols
id2name <- distinct(deseq_res_df, gene_id, gene_name)

pway_le <- gsea_riboseq_cryp %>%
  unnest_longer(leadingEdge) %>%
  left_join(id2name, by = c("leadingEdge" = "gene_id")) %>%
  group_by(pathway) %>%
  summarise(leadingEdgeName = paste(gene_name, collapse = ";"),
            leadingEdgeID = paste(leadingEdge, collapse = ";"),
            )

# join collapsed leading edge to df
gsea_riboseq_cryp_df <- gsea_riboseq_cryp %>%
  select(-leadingEdge) %>%
  left_join(pway_le, by = "pathway")

write_tsv(gsea_riboseq_cryp_df, "processed/2023-09-26_riboseq_gsea_ale_types_results.tsv", col_names = T)
saveRDS(cryptic_ale_gene_lists, "processed/gsea_riboseq.ale_gene_lists.rds")
saveRDS(riboseq_fc_ranks, "processed/gsea_riboseq.fc_ranks.rds")


