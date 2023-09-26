library(tidyverse)
library(fgsea)
set.seed(123)

#' Get a named vector of genes ranked by score from smallest to largest value
get_ranked_gene_list <- function(df, score_col = "signed_padj", name_col = "gene_name") {
  
  # Make sure smallest first (how fgsea wants it ordered)
  df <- df %>% arrange(!!sym(score_col))
  
  set_names(df[[score_col]], df[[name_col]])
  
}


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



# # generate dot plot of GSEA results - gene sets, NES as x-y, size of dot indicates significance
# plot_df_gsea_riboseq <- gsea_riboseq_cryp_df %>%
#   mutate(plot_pathway = case_when(pathway == "spliced" ~ "AS-ALE",
#                                   pathway == "distal_3utr_extension" ~ "3'UTR-ALE",
#                                   pathway == "bleedthrough" ~ "Bleedthrough-ALE")) %>%
#   mutate(sig = padj < 0.05,
#          plot_pathway = fct_reorder(plot_pathway, -log10(padj)),
#          plot_size = -log10(padj)
#   ) 
# 
# # gsea dot plot
# plot_df_gsea_riboseq %>%
#   ggplot(aes(x = NES, y = plot_pathway, colour = sig, size = plot_size)) +
#   geom_point() +
#   geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
#   scale_x_continuous(limits = c(-3, 3)) +
#   theme_bw(base_size = 20) +
#   scale_colour_manual(values = c("#d95f02", "#1b9e77")) +
#   labs(x = "GSEA normalised enrichment score",
#        y = "Gene set",
#        colour = "padj < 0.05",
#        size = "-log10(padj)")
# 
# ggsave("processed/2023-06-28_gsea_cryptics_dotplot.png",
#        height = 8,
#        width = 12,
#        units = "in")
# 
# 
# # classic enrichment plots
# plotEnrichment(cryptic_ale_gene_lists[["distal_3utr_extension"]],
#                riboseq_fc_ranks) +
#   labs(title = "3'UTR-ALE",
#        subtitle = glue::glue("NES = {round(plot_df_gsea_riboseq[2,6], 4)}, padj = {round(plot_df_gsea_riboseq[2,3], 4)}"),
#        x = "Gene Ranks",
#        y = "Enrichment Score") +
#   theme(title = element_text(size = rel(1.5)),
#         axis.text = element_text(size = rel(1.5)))
# 
# ggsave("processed/2023-06-28_gsea_3utr_ext_cryptics_enrichplot.png",
#        device = "png",
#        height = 8,
#        width = 8,
#        units = "in",
#        dpi = "retina")