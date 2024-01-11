library(tidyverse)
library(fgsea)
library(ggrepel)
library(ggrastr)
source("helpers.R")
set.seed(123)





normed_count_mtx <- read_tsv("processed/2023-05-08_i3_cortical_riboseq_sf_normed_count_matrix.tsv")
deseq_res_df <- read_tsv("processed/2023-05-08_i3_cortical_riboseq_deseq2_results.tsv")
papa_cryp_et <- read_tsv("../preprocessing/processed/cryptics_summary_all_events_complex.tsv")
# manual curation of complex event types
papa_pm_et_mc <- read_tsv("../postmortem/processed/2023-09-11_liu_facs_decoys_per_sample_delta_ppau.cryptics.tsv")
bleedthrough_mv <- read_tsv("../preprocessing/processed/bleedthrough_manual_validation.tsv")
event_type_complex_mc <- read_tsv("../postmortem/data/cryptics_summary_complex_manual_curation.tsv")

# riboseq gsea results 
gsea_riboseq_cryp_df <- read_tsv("processed/2023-09-26_riboseq_gsea_ale_types_results.tsv")


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

ggsave(last_plot(),
       filename = "2023-05-08_elk1_riboseq_counts_plot_all.svg",
       path = "processed/",
       device = svg,
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
       x = "Ribo-seq Log2FoldChange (TDP43KD / CTRL)",
       y = "-log10(padj)")


# Same plot but with (select) distal 3'UTR extensions highlighted
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
                  force = 55,size = rel(8),
                  seed = 123
  ) +
  scale_colour_manual(values = c("#d95f02", "#bdbdbd")) +
  theme_bw(base_size = 16) +
  scale_x_continuous(limits = c(-4,4),
                     breaks = seq(-4,4,1)) +
  guides(alpha = "none", colour = "none") +
  labs(title = "i3 Cortical Neuron Ribo-seq Differential expression",
       x = "Ribo-seq Log2FoldChange (TDP43KD / CTRL)",
       y = "-log10(padj)")

ggsave("processed/2023-06-29_riboseq_de_volcano_elk1_six3_tlx1.png",
       device = "png", height = 8, width = 12, units = "in",dpi = "retina")

ggsave("processed/2023-06-29_riboseq_de_volcano_elk1_six3_tlx1.svg",
       device = svg, height = 8, width = 12, units = "in",dpi = "retina")


###------
# Volcanos with all cryptic event types
###------

papa_cryp_et_i3 <- papa_cryp_et %>%
  filter(str_detect(experiment_name, "i3_cortical"))

# # load in manual curation (of bleedthroughs) & get a df of failing IDs
# 
bleedthrough_f <- bleedthrough_mv %>%
  filter(event_manual_validation != "yes")

# have manual curation of event types
# combine the two tables, filtering for passing events
event_type_complex_mc <- event_type_complex_mc %>%
  mutate(plot_le_id = paste(gene_name, str_split_i(le_id, "_", 2), sep = "_"))

event_type_complex_mc_p <- event_type_complex_mc %>%
  select(le_id, plot_le_id, gene_name, event_manual_validation, simple_event_type = manual_simple_event_type) %>%
  # some events did not pass manual curation - filter out
  filter(event_manual_validation == "yes")


# Where have re-annotated complex events, update the event type annotation
papa_cryp_et_i3 <- left_join(papa_cryp_et_i3, select(event_type_complex_mc_p, gene_name, simple_event_type),
            by = c("gene_name")) %>%
   mutate(et_upd = if_else(is.na(simple_event_type.y),
                           simple_event_type.x,
                           simple_event_type.y)) %>% 
   select(-starts_with("simple_event_type")) %>%
   rename(simple_event_type = et_upd) 

# # get a vector of gene names that have failed manual validation
# union(bleedthrough_f$gene_name,
#       filter(event_type_complex_mc, event_manual_validation != "yes") %>% pull(gene_name))

# annotate input cryptics with Riboseq DESeq information
papa_cryp_et_i3_deseq <- papa_cryp_et_i3 %>%
  left_join(select(deseq_res_df, gene_name, baseMean, log2FoldChange, log2FoldChangeShrink, pvalue, padj), by = "gene_name") %>%
  mutate(diff_translated = padj < 0.05)

# counts of signif events in each event type
papa_cryp_et_i3_sig_counts <- count(papa_cryp_et_i3_deseq, simple_event_type, diff_translated) %>%
  group_by(simple_event_type) %>%
  mutate(frac_signif = n / sum(n)) %>%
  ungroup()

# write to file
write_tsv(papa_cryp_et_i3_deseq, "processed/2023-10-18_i3_cryptic_genes_riboseq_deseq_summary.tsv", col_names = T)
write_tsv(papa_cryp_et_i3_sig_counts, "processed/2023-10-18_i3_cryptic_genes_riboseq_eventtype_sig_summary.tsv", col_names = T)


# Generate volcano plot with points highlighted if cryptic ALE containing (split by event type)
plot_df_cryp_volc <- deseq_res_df %>%
  # add event types for cryptics only
  left_join(papa_cryp_et_i3, by = "gene_name") %>%
  drop_na(padj) %>%
  # remove genes containing bleedthrough that failed manual curation
  filter(!gene_name %in% bleedthrough_f$gene_name) %>%
  mutate(simple_event_type = if_else(!is.na(simple_event_type) & str_detect(simple_event_type, ","),
                                     "complex",
                                     simple_event_type)) %>%
  # make columns for plotting
  #1. does gene contain a cryptic ALE
  #2. if not make more transparent
  #3. label gene name if Ribo-seq is significant and a cryptic ALE gene
  #4. shrink extremely small padj values for plotting
  mutate(simple_event_type = if_else(is.na(simple_event_type), "other", simple_event_type),
         plot_alpha = case_when(simple_event_type != "other" & padj < 0.05 ~ 1,
                                simple_event_type != "other" ~ 0.25,
                                padj < 0.05 ~ 0.25,
                                TRUE ~ 0.00001
         ),
         plot_label = if_else(simple_event_type != "other" & padj < 0.05,
                              gene_name,
                              ""),
         plot_padj = if_else(-log10(padj) > 10, 10, -log10(padj)),
         plot_event_type = case_when(simple_event_type == "bleedthrough" ~ "IPA",
                                     simple_event_type == "complex" ~ "Complex",
                                     simple_event_type == "distal_3utr_extension" ~ "3'Ext",
                                     # simple_event_type == "proximal_3utr_extension" ~ "Prox 3'Ext",
                                     simple_event_type == "other" ~ "None",
                                     simple_event_type == "spliced" ~ "ALE"),
         plot_event_type = factor(plot_event_type, levels = c("3'Ext", "ALE", "IPA", "Complex", "None"))
  ) 

# volcabno plot with cryptic ALE genes annotated
# add background and cryptic ALE points in layers so different event types are plotted on top of background
base_volcano <- ggplot(filter(plot_df_cryp_volc, simple_event_type == "other"),
       aes(x = log2FoldChangeShrink,
           y = plot_padj,
           colour = plot_event_type,
           label=plot_label,
           alpha=plot_alpha)) + 
  geom_point() +
  # overlay cryptic event types
  geom_point(data = filter(plot_df_cryp_volc, simple_event_type != "other"),  size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", "alpha" = 0.5) +
  scale_x_continuous(limits = c(-4,4),
                     breaks = seq(-4,4,1)) +
  scale_y_continuous(limits = c(0,10.5),
                     breaks = seq(0,10,1)) +
  scale_colour_manual(values = c("#d95f02","#1f78b4","#a6cee3", "#33a02c", "#bdbdbd")) +
  theme_bw(base_size = 20) +
  guides(alpha = "none") +
  labs(x = "Ribo-seq Log2FoldChange (TDP43KD / CTRL)",
       y = "-log10(padj)",
       colour = "Event Type") +
  theme(legend.position = "top")

base_volcano
base_volcano_rast <- rasterise(base_volcano, layers = "Point", dpi = 300)

ggsave("processed/2024-01-11_riboseq_ale_events_volcano_clean_no_gn.png",
       plot = base_volcano,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina"
)

ggsave("processed/2024-01-11_riboseq_ale_events_volcano_clean_no_gn.svg",
       plot = base_volcano,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina",
       device = svg
)

ggsave("processed/2024-01-11_riboseq_ale_events_volcano_clean_no_gn_rast.svg",
       plot = base_volcano_rast,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina",
       device = svg
)

# volcano with no colour legend
base_volcano_noguide <- base_volcano + guides(colour = "none")
base_volcano_noguide_rast <- rasterise(base_volcano_noguide, layers = "Point",dpi = 300)

ggsave("processed/2024-01-11_riboseq_ale_events_volcano_clean_no_gn_no_leg.png",
       plot = base_volcano_noguide,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina"
)

ggsave("processed/2024-01-11_riboseq_ale_events_volcano_clean_no_gn_no_leg.svg",
       plot = base_volcano_noguide,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina",
       device = svg
)

ggsave("processed/2024-01-11_riboseq_ale_events_volcano_clean_no_gn_no_leg_rast.svg",
       plot = base_volcano_noguide_rast,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina",
       device = svg
)

# volcano with text labels
base_volcano_gn <- base_volcano +
  geom_text_repel(data = filter(plot_df_cryp_volc, simple_event_type != "other"),
                  max.overlaps = 10000,
                  force = 30,
                  size = rel(4),
                  min.segment.length = 0,
                  seed = 123
  )

base_volcano_gn_rast <- rasterise(base_volcano_gn, layers = "Point", dpi = 300)

ggsave("processed/2024-01-11_riboseq_ale_events_volcano_clean_gn.png",
       plot = base_volcano_gn,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina"
)

ggsave("processed/2024-01-11_riboseq_ale_events_volcano_clean_gn.svg",
       plot = base_volcano_gn,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina",
       device = svg
)

ggsave("processed/2024-01-11_riboseq_ale_events_volcano_clean_gn_rast.svg",
       plot = base_volcano_gn_rast,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina",
       device = svg
)


# save plotting df to file
write_tsv(plot_df_cryp_volc,file = "processed/2023-09-26_riboseq_volcano_plot_df.tsv", col_names = T)

# volcano with no text labels, but larger points 
# base_volcano +
#   geom_point(data = filter(plot_df_cryp_volc, simple_event_type != "other"), size = 2.5) +
#   theme_bw(base_size = 16) +
#   guides(alpha = "none") +
#   labs(x = "Log2FoldChange (TDP43KD / CTRL)",
#        y = "-log10(padj)",
#        colour = "Event type")

# ggplot(filter(plot_df_cryp_volc, simple_event_type == "other"),
#        aes(x = log2FoldChange,
#            y = plot_padj,
#            colour = simple_event_type,
#            label=plot_label,
#            alpha=plot_alpha)) + 
#   geom_point() +
#   geom_point(data = filter(plot_df_cryp_volc, simple_event_type != "other"), size = 3) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", "alpha" = 0.5) +
#   geom_vline(xintercept = 0, linetype = "dashed", "alpha" = 0.5) +
#   scale_x_continuous(limits = c(-4,4),
#                      breaks = seq(-4,4,1)) +
#   scale_y_continuous(limits = c(0,10.5),
#                      breaks = seq(0,10,1)) +
#   scale_colour_manual(values = c("#33a02c", "#a6cee3", "#d95f02", "#bdbdbd", "#1f78b4")) +
#   theme_bw(base_size = 16) +
#   guides(alpha = "none") +
#   labs(x = "Log2FoldChange (TDP43KD / CTRL)",
#        y = "-log10(padj)",
#        colour = "Event type")
# 
# 
# ggsave("processed/2023-07-10_riboseq_no_gene_name_labels_larger_points.png",
#        height = 8,
#        width = 12,
#        units = "in",
#        dpi = "retina"
# )


volc_ribo_gn_base <- ggplot(filter(plot_df_cryp_volc, simple_event_type == "other"),
       aes(x = log2FoldChange,
           y = plot_padj,
           colour = simple_event_type,
           label=plot_label,
           alpha=plot_alpha)) + 
  geom_point() +
  geom_point(data = filter(plot_df_cryp_volc, simple_event_type != "other"), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", "alpha" = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", "alpha" = 0.5) +
  scale_x_continuous(limits = c(-4,4),
                     breaks = seq(-4,4,1)) +
  scale_y_continuous(limits = c(0,10.5),
                     breaks = seq(0,10,1)) +
  scale_colour_manual(values = c("#33a02c", "#a6cee3", "#d95f02", "#bdbdbd", "#1f78b4")) +
  theme_bw(base_size = 20) +
  guides(alpha = "none") +
  labs(x = "Ribo-seq Log2FoldChange (TDP43KD / CTRL)",
       y = "-log10(padj)",
       colour = "Event type") +
  theme(legend.position = "top")

volc_ribo_gn_base

# just elk1, six3 & tlx1
volc_ribo_gn_3utr_tfs <- volc_ribo_gn_base +
  geom_text_repel(data = filter(plot_df_cryp_volc, gene_name %in% c("ELK1", "SIX3", "TLX1")),
                  max.overlaps = 10000,
                  force = 30,
                  size = rel(8),
                  min.segment.length = 0,
                  seed = 123
  )

volc_ribo_gn_3utr_up <- volc_ribo_gn_base +
  geom_text_repel(data = filter(plot_df_cryp_volc, gene_name %in% c("ELK1", "SIX3", "TLX1", "BRINP2")),
                max.overlaps = 10000,
                force = 30,
                size = rel(8),
                min.segment.length = 0,
                seed = 123
                )

volc_ribo_gn_3utr_up

d3utr_sig_tr_gn <- plot_df_cryp_volc %>%
  filter(padj < 0.05 & simple_event_type == "distal_3utr_extension") %>%
  pull(gene_name)

# label all d3utr reg genes
volc_ribo_gn_3utr_all <- volc_ribo_gn_base +
  geom_text_repel(data = filter(plot_df_cryp_volc, gene_name %in% d3utr_sig_tr_gn),
                  max.overlaps = 10000,
                  force = 30,
                  size = rel(8),
                  min.segment.length = 0,
                  seed = 123
  )

volc_ribo_gn_3utr_all

ggsave("processed/2023-10-31_riboseq_no_labels_big_larger_points.png",
       plot = volc_ribo_gn_base,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina"
)

ggsave("processed/2023-10-31_riboseq_d3utr_elk1_six3_tlx1_labels_big_larger_points.png",
       plot = volc_ribo_gn_3utr_tfs,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina"
)

ggsave("processed/2023-10-31_riboseq_d3utr_up_labels_big_larger_points.png",
       plot = volc_ribo_gn_3utr_up,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina"
)

ggsave("processed/2023-10-31_riboseq_d3utr_all_labels_big_larger_points.png",
       plot = volc_ribo_gn_3utr_all,
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina"
)

####-----
# GSEA plots
####-----

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
  ggplot(aes(x = NES, y = plot_pathway, colour = plot_pathway, size = plot_size)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme_bw(base_size = 20) +
  scale_colour_manual(values = c("#33a02c", "#d95f02", "#1f78b4")) +
  guides(colour = "none") +
  labs(x = "GSEA normalised enrichment score",
       y = "",
       size = "-log10(padj)") +
  theme(legend.position = "top")

ggsave("processed/2023-10-20_gsea_cryptics_dotplot.png",
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina")

ggsave("processed/2023-10-20_gsea_cryptics_dotplot.svg",
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina",
       device = svg)

ggsave("processed/2023-10-20_gsea_cryptics_dotplot_wide.svg",
       height = 5.33,
       width = 8,
       units = "in",
       dpi = "retina",
       device = svg)

# classic GSEA enrichment plots

# read in GSEA input gene lists for different ALE types (as )
cryptic_ale_gene_lists <- readRDS("processed/gsea_riboseq.ale_gene_lists.rds")
riboseq_fc_ranks <- readRDS("processed/gsea_riboseq.fc_ranks.rds")


plot_gsea_line <- function(gene_list, ranks, nes, padj, plot_title = "",round_digits = 5) {
  
  # print(gsea_df)
  # # extract NES & padj from gsea_df
  # pway_df <- gsea_df %>%
  #   dplyr::filter(pathway %in% plot_pathway)
  # 
  # padj <- round(pway_df$padj, digits = round_digits)
  # nes <- round(pway_df$NES, digits = round_digits)
  # print(padj)
  # print(nes)
  # print(pway_df)
  
  plotEnrichment(gene_list,
                 ranks) +
  labs(title = plot_title,
         subtitle = glue::glue("NES = {round(nes, round_digits)}, padj = {round(padj, round_digits)}"),
         x = "Gene Ranks",
         y = "Enrichment Score") +
    theme(title = element_text(size = rel(1.5)),
          axis.text = element_text(size = rel(1.5))
          )
    
}


plot_gsea_line(gene_list = cryptic_ale_gene_lists[["distal_3utr_extension"]],
              ranks = riboseq_fc_ranks,nes = 1.5, padj = 0.03,plot_title = "3'UTR-ALE")

tmp_grpd <- plot_df_gsea_riboseq %>%
  group_by(plot_pathway)

gsea_pway_grpd <- group_split(tmp_grpd) %>%
  set_names(group_keys(tmp_grpd) %>% pull())  # group keys returns 1 col df in case of single group

gsea_enrichplots <- pmap(list(x = cryptic_ale_gene_lists,
          y = gsea_pway_grpd,
          z = names(gsea_pway_grpd)),
     function(x,y, z) plot_gsea_line(x, riboseq_fc_ranks, y$NES, y$padj, plot_title = z)
     )

gsea_enrichplots <- gsea_enrichplots %>%
  set_names(names(gsea_pway_grpd))

walk2(.x = gsea_enrichplots,
      .y = names(gsea_enrichplots),
      ~ ggsave(paste("processed/2023-09-26_gsea_cryptics_enrichplot.",
               str_replace_all(.y, "'|-", "_"),
               ".png",
               sep = ""),
               plot = .x,
         device = "png",
       height = 8,
       width = 12,
       units = "in",
       dpi = "retina")
      )

walk2(.x = gsea_enrichplots,
      .y = names(gsea_enrichplots),
      ~ ggsave(paste("processed/2023-09-26_gsea_cryptics_enrichplot.",
                     str_replace_all(.y, "'|-", "_"),
                     ".svg",
                     sep = ""),
               plot = .x,
               device = svg,
               height = 8,
               width = 12,
               units = "in",
               dpi = "retina")
)


## TSVs


####-------
# Are ribo-seq significant cryptics also changed in corresponding direction on RNA level?
####-------
seddighi_rna <- read_csv("data/rnaseq/seddighi.ipscCortical_neuron.DESEQ2_results.csv")

plot_df_cryp_volc_rna <- plot_df_cryp_volc %>%
  # filter for cryptic-containing genes that are significant only
  left_join(select(seddighi_rna, gene_name, log2FoldChange, padj), suffix = c(".ribo", ".rna"), by = "gene_name") %>%
  mutate(has_cryptic = if_else(simple_event_type == "other", "no_cryptic", "cryptic"),
         riboseq_sig = if_else(padj.ribo < 0.05, T, F),
         rnaseq_sig = if_else(padj.rna < 0.05, T, F)) 

plot_df_cryp_volc_rna %>%
  # drop genes where rna levels filtered by DESEQ2 independent filtering
  # drop_na(padj.rna) %>%
  count(has_cryptic, riboseq_sig, rnaseq_sig) 

# A tibble: 11 × 4
# has_cryptic riboseq_sig rnaseq_sig     n
# <chr>       <lgl>       <lgl>      <int>
#   1 cryptic     FALSE       FALSE         35
# 2 cryptic     FALSE       TRUE          52
# 3 cryptic     FALSE       NA             1
# 4 cryptic     TRUE        FALSE          1
# 5 cryptic     TRUE        TRUE          20
# 6 no_cryptic  FALSE       FALSE       4712
# 7 no_cryptic  FALSE       TRUE        5239
# 8 no_cryptic  FALSE       NA           287
# 9 no_cryptic  TRUE        FALSE        149
# 10 no_cryptic  TRUE        TRUE         626
# 11 no_cryptic  TRUE        NA            36

# All but one cryptic genes sig on ribo-seq are also sig on differential expression
# Enrichment also true for genes not containing cryptics, but to slightly lesser extent 

# what about directionality of change? Are RNA and translation changes consistent or opposing?
# Focus on riboseq DE only
plot_df_cryp_volc_rna %>%
  # filter for cryptic-containing genes that are significant only
  filter(plot_label != "") %>%
  left_join(select(seddighi_rna, gene_name, log2FoldChange, padj),suffix = c(".ribo", ".rna"), by = "gene_name") %>%
  mutate(rna_padj_lt_05 = if_else(padj.rna < 0.05, T, F)) %>%
  ggplot(aes(y = log2FoldChange.ribo, x = log2FoldChange.rna, colour = simple_event_type, shape = rna_padj_lt_05)) +
  geom_point(size = rel(3)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(limits = c(-3, 3)) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_colour_manual(values = c("#d95f02","#33a02c", "#1f78b4")) +
  theme_bw(base_size = 14) +
  labs(
       x = "Ribo-seq Log2FoldChange(KD/WT)",
       y = "RNA-seq Log2FoldChange(KD/WT)",
       shape = "RNA-seq DE padj < 0.05",
       colour = "ALE Event Type")

ggsave("processed/2023-09-05_riboseq_cryptic_sig_rna_fc_direction_scatter.png",
       device = "png",
       height = 8,
       width = 8,
       units = "in",
       dpi = "retina")

# Should repeat with cryptic containing but ns ribo-seq - are they differentially coinciding with significant rna changes?
# 52 cryptics sig on total RNA levels but not riboseq, 35 ns on both. So by eye clearly a shift...
tibble(riboseq_sig = c("riboseq_sig", "riboseq_ns"), rnaseq_sig = c(20, 52), rnaseq_ns = c(1, 35)) %>%
column_to_rownames("riboseq_sig") %>%
 fisher.test()
# FET suggests that cryptics with riboseq sig changes are enriched for differential expression on RNA level
# Consistent directionality?

