library(tidyverse)
library(ggrepel)
set.seed(123)

#' Tidy up DESeq2 results table for analyses of enrichment
#' Removes low count genes (those filtered out by DESeq's independent filtering) & ensures 1 row per gene 
clean_deseq_df <- function(df, padj_col = "padj", id_col = "gene_name") {
  
  # Remove genes with NA padj - i.e. too low counts for stat testing &/ filtered out due to low counts (to maximize n sig)
  df %>%
    drop_na(!!sym(padj_col)) %>%
    distinct(!!sym(id_col), .keep_all = T)
  
}


seddighi_df <- read_csv("data/seddighi.ipscCortical_neuron.DESEQ2_results.csv")

# manual curation
bleedthrough_mv <- read_tsv("../preprocessing/processed/bleedthrough_manual_validation.tsv")
# event_type_complex_mc <- read_tsv("../postmortem/data/cryptics_summary_complex_manual_curation.tsv")

# use final definition of cryptic events
cryp_df <- read_tsv("../preprocessing/processed/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv")




# genes with le_ids that fail manual validation
bleedthrough_mv_f_nm <- bleedthrough_mv %>% filter(event_manual_validation != "yes") %>% pull(gene_name)
# event_type_complex_mc_f_nm <- event_type_complex_mc %>% filter(event_manual_validation != "yes") %>% pull(gene_name)
# mv_fail_gn <- c(bleedthrough_mv_f_nm, event_type_complex_mc_f_nm)


seddighi_df <- clean_deseq_df(seddighi_df)


# add in information about whether a cryptic containing gene
seddighi_df_all <- cryp_df %>%
  # subset for i3 cortical cryptics only
  filter(str_detect(experiment_name, "i3_cortical")) %>%
  distinct(gene_name, simple_event_type) %>%
  left_join(seddighi_df, ., by = "gene_name")

seddighi_df <- cryp_df %>%
  # subset for seddighi i3 cortical cryptics only
  filter(str_detect(experiment_name, "seddighi_i3_cortical")) %>%
  distinct(gene_name, simple_event_type) %>%
  left_join(seddighi_df, ., by = "gene_name")

# annotate as cryptic containing, sig diff exprn + direction
seddighi_df <- seddighi_df %>%
  mutate(cryptic = !is.na(simple_event_type),
         sig = padj < 0.05,
         dirn = sign(log2FoldChange)) 

seddighi_df_all <- seddighi_df_all %>%
  mutate(cryptic = !is.na(simple_event_type),
         sig = padj < 0.05,
         dirn = sign(log2FoldChange)) 

# get number of ALE containing genes that are up/downregulated
de_cryptics_gene_counts <- seddighi_df %>%
  filter(cryptic) %>%
  count(sig, dirn)

de_cryptics_all_gene_counts <- seddighi_df_all %>%
  filter(cryptic) %>%
  count(sig, dirn)

# get differentially expressed gene counts split by cryptic type
de_cryptics_gene_counts_event <- seddighi_df %>%
  filter(cryptic & sig) %>%
  count(simple_event_type, sig, dirn)


de_cryptics_all_gene_counts_event <- seddighi_df_all %>%
  filter(cryptic & sig) %>%
  count(simple_event_type, sig, dirn)

write_tsv(de_cryptics_gene_counts, "processed/2023-12-22_seddighi_cryptic_seddighi_diff_expressed_gene_counts.tsv")
write_tsv(de_cryptics_all_gene_counts, "processed/2023-12-22_i3_cryptic_seddighi_diff_expressed_gene_counts.tsv")
write_tsv(de_cryptics_gene_counts_event, "processed/2023-12-22_seddighi_cryptic_seddighi_diff_expressed_gene_counts_event_type.tsv",col_names = T)
write_tsv(de_cryptics_all_gene_counts_event, "processed/2023-12-22_i3_cryptic_seddighi_diff_expressed_gene_counts_event_type.tsv",col_names = T)




# RNA volcano with all APA genes highlighted
plot_seddighi_df_all <- seddighi_df_all %>%
  mutate(plot_padj = if_else(-log10(padj) > 50, 50, -log10(padj)),
         plot_alpha = case_when(cryptic & sig ~ 1,
                                cryptic ~ 0.3,
                                plot_padj > -log10(0.05) ~ 0.1,
                                TRUE ~ 0.01
                                ),
         plot_event_type = case_when(simple_event_type == "distal_3utr_extension" ~ "3'Ext",
                                     simple_event_type == "bleedthrough" ~ "IPA",
                                     simple_event_type == "spliced" ~ "ALE",
                                     T ~ "Other"),
         plot_event_type = factor(plot_event_type, levels = c("ALE", "IPA", "3'Ext", "Other"))
         ) 


seddighi_volcano_cryptic_all <- ggplot(filter(plot_seddighi_df_all, plot_event_type == "Other"),
       aes(x = log2FoldChange,
           y = plot_padj,
           colour = plot_event_type,
           alpha=plot_alpha)) +
  geom_point() +
  geom_point(data = filter(plot_seddighi_df_all, plot_event_type != "Other"), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", "alpha" = 0.5) +
  scale_colour_manual(values = c("#d95f02","#1f78b4", "#33a02c", "#bdbdbd")) +
  theme_bw(base_size = 20) +
  scale_x_continuous(limits = c(-5,5),
                     breaks = seq(-10,10,1)) +
  guides(alpha = "none") +
  labs(colour = "Event Type",
       x = "Log2FoldChange (KD / WT)",
       y = "-log10(padj)") +
  theme(legend.position = "top")

seddighi_volcano_cryptic_all

ggsave("processed/2023-12-22_seddighi_rna_de_volcano_all_apa_colour_nolab.svg",
       device = svg,
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")


# volcano with specific genes highlighted
# here - all 3'UTR extension cryptics with increased overall translation
highlight_genes <- c("ELK1", "SIX3", "TLX1", "BRINP2")

plot_seddighi_cryp_3utr <- seddighi_df %>%
  mutate(plot_label = if_else(gene_name %in% highlight_genes,
                              gene_name,
                              ""),
         plot_padj = if_else(-log10(padj) > 50, 50, -log10(padj)),
         plot_alpha = case_when(plot_label %in% highlight_genes ~ 5,
                                plot_padj > -log10(0.05) ~ 0.1,
                                TRUE ~ 0.01
         ),
         plot_colour = if_else(plot_label %in% highlight_genes, "3'UTR-ALE cryptics", "other")
  ) 

min(plot_seddighi_cryp_3utr$log2FoldChange)
max(plot_seddighi_cryp_3utr$log2FoldChange)


ggplot(filter(plot_seddighi_cryp_3utr, plot_colour == "other"),
         aes(x = log2FoldChange,
             y = plot_padj,
             colour = plot_colour,
             label=plot_label,
             alpha=plot_alpha)) +
  geom_point() +
  geom_point(data = filter(plot_seddighi_cryp_3utr, plot_colour != "other")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", "alpha" = 0.5) +
  geom_text_repel(data = filter(plot_seddighi_cryp_3utr, plot_colour != "other"),
                  max.overlaps = 1000,
                  force = 10,
                  size = rel(8),
                  seed = 123
  ) +
  scale_colour_manual(values = c("#d95f02", "#bdbdbd")) +
  theme_bw(base_size = 20) +
  scale_x_continuous(limits = c(-5,5),
                     breaks = seq(-10,10,1)) +
  guides(alpha = "none", colour = "none") +
  labs(
       x = "Log2FoldChange (KD / WT)",
       y = "-log10(padj)")

ggsave("processed/2023-10-10_seddighi_rna_de_volcano_cryptic_3utrs_lab.png",
       device = "png",
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")

ggsave("processed/2023-10-10_seddighi_rna_de_volcano_cryptic_3utrs_lab.svg",
       device = svg,
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")



## Volcano plot labelling all diff translated genes (from ribo-seq) on RNA levels
riboseq <- read_tsv("data/2023-10-18_i3_cryptic_genes_riboseq_deseq_summary.tsv")
riboseq_sig <- filter(riboseq, diff_translated & !(gene_name %in% mv_fail_gn))


# df specifying plot criteria - labels for diff translated genes containing cryptics
plot_seddighi_cryp_all <- seddighi_df %>%
  left_join(select(riboseq_sig, gene_name, simple_event_type, diff_translated), by = "gene_name") %>%
  mutate(plot_label = if_else(diff_translated & padj < 0.05 & !(gene_name %in% mv_fail_gn),
                              gene_name,
                              ""),
         plot_padj = if_else(-log10(padj) > 50, 50, -log10(padj)),
         plot_alpha = case_when(diff_translated & !(gene_name %in% mv_fail_gn) ~ 5,
                                plot_padj > -log10(0.05) ~ 0.1,
                                TRUE ~ 0.01
         ),
         plot_event_type = case_when(simple_event_type == "bleedthrough" ~ "IPA",
                                     simple_event_type == "distal_3utr_extension" ~ "3'Ext",
                                     simple_event_type == "spliced" ~ "ALE"),
         plot_event_type = factor(plot_event_type, levels = c("ALE","IPA","3'Ext"))
         # plot_colour = if_else(plot_label %in% highlight_genes, "3'UTR-ALE cryptics", "other")
  ) %>%
  arrange(desc(diff_translated))

# base plot with just colours by event type
base_ev_type_volc <- ggplot(filter(plot_seddighi_cryp_all, is.na(diff_translated)),
       aes(x = log2FoldChange,
           y = plot_padj,
           colour = plot_event_type,
           label=plot_label,
           alpha=plot_alpha)) +
  geom_point(colour = "#bdbdbd") +
  geom_point(data = filter(plot_seddighi_cryp_all, diff_translated), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", "alpha" = 0.5) +
  scale_colour_manual(name = "",
                      values = c( "#1f78b4","#33a02c","#d95f02")) +
  scale_x_continuous(limits = c(-5,5),
                     breaks = seq(-10,10,1)) +
  theme_bw(base_size = 20) +
  guides(alpha = "none", label = "none") +
  labs(
    x = "Log2FoldChange (KD / WT)",
    y = "-log10(padj)") +
  theme(legend.position = "top")

base_ev_type_volc

# colour all event types, label all genes
all_gn_ev_type_volc <- base_ev_type_volc +
  geom_text_repel(data = filter(plot_seddighi_cryp_all, diff_translated),
                max.overlaps = 1000,
                force = 10,
                size = rel(6),
                seed = 123
) 

# COLOUR ALL event types but just label the 3'UTR-ALEs
utr3_gn_ev_type_volc <- base_ev_type_volc +
  geom_text_repel(data = filter(plot_seddighi_cryp_all, diff_translated & gene_name %in% highlight_genes),
                  max.overlaps = 1000,
                  force = 10,
                  size = rel(8),
                  seed = 123
  )

base_ev_type_volc
all_gn_ev_type_volc
utr3_gn_ev_type_volc

# base plot, no labels
ggsave("processed/2023-12-13_seddighi_rna_de_volcano_cryptic_no_lab.png",
       plot = base_ev_type_volc,
       device = "png",
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")

ggsave("processed/2023-12-13_seddighi_rna_de_volcano_cryptic_no_lab.svg",
       plot = base_ev_type_volc,
       device = svg,
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")

# all labelled
ggsave("processed/2023-12-13_seddighi_rna_de_volcano_cryptic_all_lab.png",
       plot = all_gn_ev_type_volc,
       device = "png",
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")

ggsave("processed/2023-12-13_seddighi_rna_de_volcano_cryptic_all_lab.svg",
       plot = all_gn_ev_type_volc,
       device = svg,
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")

# only increased 3'UTRs up
ggsave("processed/2023-12-13_seddighi_rna_de_volcano_cryptic_3utr_lab.png",
       plot = utr3_gn_ev_type_volc,
       device = "png",
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")

ggsave("processed/2023-12-13_seddighi_rna_de_volcano_cryptic_3utr_lab.svg",
       plot = utr3_gn_ev_type_volc,
       device = svg,
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")
