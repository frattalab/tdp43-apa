library(tidyverse)
library(fgsea)
library(ggrepel)
set.seed(123)

#' Get a named vector of genes ranked by score from smallest to largest value
get_ranked_gene_list <- function(df, score_col = "signed_padj", name_col = "gene_name") {
  
  # Make sure smallest first (how fgsea wants it ordered)
  df <- df %>% arrange(!!sym(score_col))
  
  set_names(df[[score_col]], df[[name_col]])
  
}



normed_count_mtx <- read_tsv("processed/2023-05-08_i3_cortical_riboseq_sf_normed_count_matrix.tsv")
deseq_res_df <- read_tsv("processed/2023-05-08_i3_cortical_riboseq_deseq2_results.tsv")
papa_cryp_et <- read_tsv("../preprocessing/processed/cryptics_summary_all_events_complex.tsv")

# ELK1 riboseq counts plot
plot_df_elk1_nc <- normed_count_mtx %>%
filter(gene_name == "ELK1") %>%
pivot_longer(cols = contains("sh"), values_to = "normed_counts", names_to = "sample_id") %>%
mutate(condition = if_else(str_detect(sample_id, "tdp"), "KD", "CTL"))

plot_df_elk1_nc %>%
  ggplot(aes(x = condition, y = normed_counts, colour = condition)) +
  geom_jitter(width = 0.25, size = rel(2.5)) +
  scale_colour_manual(name = "Condition",
                      values = c("#1b9e77", "#d95f02")) +
  scale_y_continuous(breaks = seq(20,120,10)) +
  labs(title = "ELK1 Ribo-seq CDS counts",
       subtitle = "Counts adjusted for variable library depth/composition with DESeq2's size factor method",
       x = "Condition",
       y = "CDS counts") +
  theme_bw() +
  theme(axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        title = element_text(size = rel(1.5)))


ggsave(last_plot(),
       filename = "2023-05-08_elk1_riboseq_counts_plot_all.png",
       path = "processed/",
       device = "png",
       width = 10,
       height = 10,
       units = "in",
       dpi = "retina")

# ELK1 
plot_df_elk1_volc <- deseq_res_df %>%
  drop_na(padj) %>%
  # mutate(padj = replace_na(padj, 1)) %>%
  mutate(plot_label = if_else(gene_name %in% "ELK1",
                              gene_name,
                              ""),
         plot_padj = if_else(-log10(padj) > 10, 10, -log10(padj)),
         plot_alpha = case_when(plot_label == "ELK1" ~ 1,
                                plot_padj > -log10(0.05) ~ 0.8,
                                TRUE ~ 0.1
         ),
         plot_colour = if_else(plot_label == "ELK1", "ELK1", "other")
  )

plot_df_elk1_volc %>%
  ggplot(aes(x = log2FoldChangeShrink,
             y = plot_padj,
             colour = plot_colour,
             label=plot_label,
             alpha=plot_alpha)) + 
  geom_point() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", "alpha" = 0.5) +
  geom_text_repel(max.overlaps = 50,
                  force = 55,size = rel(8)
  ) +
  scale_colour_manual(values = c("#d95f02", "#bdbdbd")) +
  theme_bw(base_size = 16) +
  scale_x_continuous(limits = c(-4,4),
                     breaks = seq(-4,4,1)) +
  guides(alpha = "none", colour = "none") +
  labs(title = "i3 Cortical Neuron Riboseq Differential expression",
       x = "Log2FoldChange (KD / WT)",
       y = "-log10(padj)")


# Same plot but with distal 3'UTR extensions highlighted
plot_df_d3utr_volc <- deseq_res_df %>%
  drop_na(padj) %>%
  # mutate(padj = replace_na(padj, 1)) %>%
  mutate(plot_label = if_else(gene_name %in% c("ELK1","SIX3", "TLX1"),
                              gene_name,
                              ""),
         plot_padj = if_else(-log10(padj) > 10, 10, -log10(padj)),
         plot_alpha = case_when(plot_label %in% c("ELK1","SIX3", "TLX1") ~ 1,
                                plot_padj > -log10(0.05) ~ 0.8,
                                TRUE ~ 0.1
         ),
         plot_colour = if_else(plot_label %in% c("ELK1","SIX3", "TLX1"), "Cryptic_3'UTR-ALE", "other")
  )


plot_df_d3utr_volc %>%
  ggplot(aes(x = log2FoldChangeShrink,
             y = plot_padj,
             colour = plot_colour,
             label=plot_label,
             alpha=plot_alpha)) + 
  geom_point() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", "alpha" = 0.5) +
  geom_text_repel(max.overlaps = 50,
                  force = 55,size = rel(8),seed = 123
  ) +
  scale_colour_manual(values = c("#d95f02", "#bdbdbd")) +
  theme_bw(base_size = 16) +
  scale_x_continuous(limits = c(-4,4),
                     breaks = seq(-4,4,1)) +
  guides(alpha = "none", colour = "none") +
  labs(title = "i3 Cortical Neuron Ribo-seq Differential expression",
       x = "Log2FoldChange (KD / WT)",
       y = "-log10(padj)")

ggsave("processed/2023-06-29_riboseq_de_volcano_elk1_six3_tlx1.png", device = "png", height = 8, width = 12, units = "in",dpi = "retina")



### volcnaos with cryptic event types

papa_cryp_et_i3 <- papa_cryp_et %>%
  filter(str_detect(experiment_name, "i3_cortical"))


# load in manual curation (of bleedthroughs) & get a df of failing IDs
bleedthrough_mv <- read_tsv("../preprocessing/processed/bleedthrough_manual_validation.tsv")
bleedthrough_f <- bleedthrough_mv %>%
  filter(event_manual_validation != "yes")

# Generate volcano plot with points highlighted if cryptic ALE containing (split by event type)
plot_df_cryp_volc <- deseq_res_df %>%
  # add event types
  left_join(papa_cryp_et_i3, by = "gene_name") %>%
  drop_na(padj) %>%
  # remove genes containing bleedthrough that failed manual curation
  filter(!gene_name %in% bleedthrough_f$gene_name) %>%
  mutate(simple_event_type = if_else(!is.na(simple_event_type) & str_detect(simple_event_type, ","),
                                     "complex",
                                     simple_event_type)) %>%
  # make columns for plotting
  mutate(simple_event_type = if_else(is.na(simple_event_type), "other", simple_event_type),
         plot_alpha = case_when(simple_event_type != "other" ~ 1,
                                padj < 0.05 ~ 0.5,
                                TRUE ~ 0.00001
         ),
         plot_label = if_else(simple_event_type != "other" & padj < 0.05,
                              gene_name,
                              ""),
         plot_padj = if_else(-log10(padj) > 10, 10, -log10(padj)),
  ) 

# volcabno plot with cryptic ALE genes annotated
# add background and cryptic ALE points in layers so different event types are plotted on top of background
ggplot(filter(plot_df_cryp_volc, simple_event_type == "other"),
       aes(x = log2FoldChangeShrink,
           y = plot_padj,
           colour = simple_event_type,
           label=plot_label,
           alpha=plot_alpha)) + 
  geom_point() +
  geom_point(data = filter(plot_df_cryp_volc, simple_event_type != "other")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", "alpha" = 0.5) +
  scale_x_continuous(limits = c(-4,4),
                     breaks = seq(-4,4,1)) +
  scale_y_continuous(limits = c(0,10.5),
                     breaks = seq(0,10,1)) +
  scale_colour_manual(values = c("#d95f02", "#a6cee3", "#33a02c", "#bdbdbd", "#1f78b4")) +
  geom_text_repel(data = filter(plot_df_cryp_volc, simple_event_type != "other"),
                  max.overlaps = 10000,
                  force = 30,
                  size = rel(3),
  ) +
  theme_bw(base_size = 16) +
  guides(alpha = "none") +
  labs(x = "Log2FoldChange (KD / WT)",
       y = "-log10(padj)",
       colour = "Event type")

ggsave("processed/2023-06-28_riboseq_ale_events_volcano_clean.png",
       height = 8,
       width = 12,
       units = "in",
       dpi = "retina"
)


### GSEA on FC of different event types
# i.e. are the subset of ALE event type genes that show consistent directionality significant with respect to random gene sets of the same size?


# GSEA on FC ranks
# make gene list of different event types
tmp_grp <- plot_df_cryp_volc %>%
  filter(!simple_event_type %in% c("other", "complex")) %>%
  group_by(simple_event_type)

cryptic_ale_gene_lists <- tmp_grp %>%
  group_split() %>%
  set_names(group_keys(tmp_grp) %>% pull()) %>% # group keys returns 1 col df in case of single group
  map(~ unique(.x$gene_name))


# get a named vector of genes ranked by ribo-seq fold change 
riboseq_fc_ranks <- get_ranked_gene_list(deseq_res_df, score_col = "log2FoldChangeShrink", name_col = "gene_name")

# run standard GSEA
gsea_riboseq_cryp <- fgsea(pathways = cryptic_ale_gene_lists,
                           stats = riboseq_fc_ranks)

# collapse leading edge genes to comma separated string, then join to results df
pway_le <- gsea_riboseq_cryp %>%
  unnest_longer(leadingEdge) %>%
  group_by(pathway) %>%
  summarise(leadingEdge = paste(leadingEdge, collapse = ","))

gsea_riboseq_cryp_df <- gsea_riboseq_cryp %>%
  select(-leadingEdge) %>%
  left_join(pway_le, by = "pathway")

# generate dot plot of GSEA results - gene sets, NES as x-y, size of dot indicates significance
plot_df_gsea_riboseq <- gsea_riboseq_cryp_df %>%
  mutate(plot_pathway = case_when(pathway == "spliced" ~ "AS-ALE",
                                  pathway == "distal_3utr_extension" ~ "3'UTR-ALE",
                                  pathway == "bleedthrough" ~ "Bleedthrough-ALE")) %>%
  mutate(sig = padj < 0.05,
         plot_pathway = fct_reorder(plot_pathway, -log10(padj)),
         plot_size = -log10(padj)
  ) 

# gsea dot plot
plot_df_gsea_riboseq %>%
  ggplot(aes(x = NES, y = plot_pathway, colour = sig, size = plot_size)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(limits = c(-3, 3)) +
  theme_bw(base_size = 20) +
  scale_colour_manual(values = c("#d95f02", "#1b9e77")) +
  labs(x = "GSEA normalised enrichment score",
       y = "Gene set",
       colour = "padj < 0.05",
       size = "-log10(padj)")

ggsave("processed/2023-06-28_gsea_cryptics_dotplot.png",
       height = 8,
       width = 12,
       units = "in")


# classic enrichment plots
plotEnrichment(cryptic_ale_gene_lists[["distal_3utr_extension"]],
               riboseq_fc_ranks) +
  labs(title = "3'UTR-ALE",
       subtitle = glue::glue("NES = {round(plot_df_gsea_riboseq[2,6], 4)}, padj = {round(plot_df_gsea_riboseq[2,3], 4)}"),
       x = "Gene Ranks",
       y = "Enrichment Score") +
  theme(title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5)))

ggsave("processed/2023-06-28_gsea_3utr_ext_cryptics_enrichplot.png",
       device = "png",
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina")
