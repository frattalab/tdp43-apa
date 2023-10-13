library(tidyverse)
library(ggrepel)
library(ggpubr)
library(glue)

calc_means <- function(df, sample_cols, out_col) {
  
  df[, out_col] <- rowMeans(df[, sample_cols])
  df

}

prep_means_df_plot <- function(df, mean_col) {
  
  df %>%
    mutate(gene_name = fct_reorder(gene_name, !!sym(mean_col)),
           plot_gn = if_else(enr_cryp, gene_name, ""),
           plot_x = log2(!!sym(mean_col) + 1)) 
}

#' Generate a ggplot of log2 transformed counts for each gene in ascending order
#' must contain gene_name, enr_cryp column (T/F - whether to highlight) & a column of your choice with mean value
base_dotplot <- function(df, mean_col, plot_title = "Cryptic-containing gene expression") {
  
  df %>%
    prep_means_df_plot(., mean_col) %>%
    ggplot(aes(x = plot_x, y = gene_name, colour = enr_cryp, label = plot_gn)) +
    geom_point(size = rel(0.3)) +
    geom_vline(xintercept = log2(11), linetype = "dashed") +
    geom_vline(xintercept = log2(26), linetype = "dashed") +
    geom_vline(xintercept = log2(51), linetype = "dashed") +
    scale_colour_manual(values = c("#bdbdbd", "#d95f02")) +
    theme_classic(base_size = 14) +
    labs(title = plot_title,
         subtitle = "Vertical lines correpond to a mean of 10, 25 & 50 reads respectively",
         x = "log2(mean count + 1)",
         y = "Gene",
         colour = "Passes enrichment criteria") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "top")
  
}

#'
#' must contain gene_name, enr_cryp column (T/F - whether to highlight) & a column of your choice with mean value
ge_density_plot <- function(df, mean_col, plot_title = "Cryptic-containing gene expression") {
  
  df %>%
    prep_means_df_plot(., mean_col) %>%
    ggplot(aes(x = plot_x, fill = enr_cryp, colour = enr_cryp)) +
    geom_density(alpha = 0.5) +
    scale_colour_manual(values = c("#bdbdbd", "#d95f02")) +
    scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
    theme_classic(base_size = 14) +
    labs(title = plot_title,
         x = "log2(mean count + 1)",
         fill = "Passes enrichment criteria",
         colour = "Passes enrichment criteria"
    ) +
    theme(legend.position = "top")
  
}

#' must contain gene_name, enr_cryp column (T/F - whether to highlight) & a column of your choice with mean value
ge_boxplot <- function(df, mean_col, plot_title = "Cryptic-containing gene expression") {

  df %>%
    prep_means_df_plot(., mean_col) %>%
    ggplot(aes(y = plot_x, x = enr_cryp, colour = enr_cryp)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() +
    ggpubr::stat_compare_means() +
    scale_colour_manual(values = c("#bdbdbd", "#d95f02")) +
    scale_fill_manual(values = c("#bdbdbd", "#d95f02")) +
    theme_classic(base_size = 14) +
    labs(title = plot_title,
         x = "Passes enrichment criteria",
         y = "log2(mean count + 1)"
    )  +
    guides(fill = "none",
           colour = "none")
  
}

# decoy tx quantification
ppau_delta_paired_median_all_cryp <- read_tsv("processed/2023-09-11_liu_facs_decoys_delta_ppau.all_samples.cryptics.tsv")

# gene expression matrices
liu_abundance <- read_tsv("processed/gene_exprn/2023-09-07_liu_facs_salmon_summarised_abundance.tsv")
liu_counts <- read_tsv("processed/gene_exprn/2023-09-07_liu_facs_salmon_summarised_counts.tsv")
liu_counts_norm <- read_tsv("processed/gene_exprn/2023-09-07_liu_facs_salmon_summarised_counts_deseq_normalised.tsv")

# cryptic-containing genes
cryptics_summary <- read_tsv("data/cryptics_summary_all_events.tsv")
cryptic_gn <- distinct(cryptics_summary, gene_name) %>% pull()


# all cols minus gene_id and gene_name
sample_cols <- colnames(liu_abundance)[3:ncol(liu_abundance)]
# all TDP-negative samples
neg_sample_cols <- colnames(liu_abundance)[str_detect(colnames(liu_abundance), "_negative$")]
# all TDP-positive samples
pos_sample_cols <- colnames(liu_abundance)[str_detect(colnames(liu_abundance), "_positive$")]

# calculate condition-wise & overall means for all three categories
liu_abundance <- calc_means(liu_abundance, sample_cols, "mean_all") 
liu_abundance <- calc_means(liu_abundance, neg_sample_cols, "mean_neg") 
liu_abundance <- calc_means(liu_abundance, pos_sample_cols, "mean_pos") 

liu_counts <- calc_means(liu_counts, sample_cols, "mean_all") 
liu_counts <- calc_means(liu_counts, neg_sample_cols, "mean_neg") 
liu_counts <- calc_means(liu_counts, pos_sample_cols, "mean_pos") 

liu_counts_norm <- calc_means(liu_counts_norm, sample_cols, "mean_all") 
liu_counts_norm <- calc_means(liu_counts_norm, neg_sample_cols, "mean_neg") 
liu_counts_norm <- calc_means(liu_counts_norm, pos_sample_cols, "mean_pos") 

# filter for cryptics
liu_abundance_cryp <- filter(liu_abundance, gene_name %in% cryptic_gn)
liu_counts_cryp <- filter(liu_counts, gene_name %in% cryptic_gn)
liu_counts_norm_cryp <- filter(liu_counts_norm, gene_name %in% cryptic_gn)

# combine into single df of means
liu_all_means <-  bind_rows(counts_norm = select(liu_counts_norm_cryp, gene_id, gene_name, starts_with("mean_")),
          counts = select(liu_counts_cryp, gene_id, gene_name, starts_with("mean_")),
          tpm = select(liu_abundance_cryp, gene_id, gene_name, starts_with("mean_")),
          .id = "abundance_type")


# standard counts filter - try a few different mean cutoffs (how many excluded?)

low_counts_mean_all_ngenes <- seq(0,100, 10) %>%
  set_names() %>%
  map(~ liu_all_means %>%
        filter(abundance_type == "counts_norm",
               mean_all <= .x) %>%
        summarise(n_genes = n_distinct(gene_name))
      )  %>%
  bind_rows(.id = "max_counts")
# max_counts n_genes
# <chr>        <int>
#   1 0                1
# 2 10              11
# 3 20              14
# 4 30              15
# 5 40              20
# 6 50              22
# 7 60              28
# 8 70              28
# 9 80              31
# 10 90              36
# 11 100             39


low_counts_mean_neg_ngenes <- seq(0,100, 10) %>%
  set_names() %>%
  map(~ liu_all_means %>%
        filter(abundance_type == "counts_norm",
               mean_neg <= .x) %>%
        summarise(n_genes = n_distinct(gene_name))
  )  %>%
  bind_rows(.id = "max_counts")

low_counts_mean_pos_ngenes <- seq(0,100, 10) %>%
  set_names() %>%
  map(~ liu_all_means %>%
        filter(abundance_type == "counts_norm",
               mean_pos <= .x) %>%
        summarise(n_genes = n_distinct(gene_name))
  )  %>%
  bind_rows(.id = "max_counts")


write_tsv(low_counts_mean_all_ngenes, "processed/gene_exprn/2023-10-13_gene_exprn_low_count_filter_mean_all.tsv", col_names = T)
write_tsv(low_counts_mean_neg_ngenes, "processed/gene_exprn/2023-10-13_gene_exprn_low_count_filter_mean_neg.tsv", col_names = T)
write_tsv(low_counts_mean_pos_ngenes, "processed/gene_exprn/2023-10-13_gene_exprn_low_count_filter_mean_pos.tsv", col_names = T)

# ~ 40 / 283 have GENE-level counts < 100. I.e. unique fragments across the whole gene, so could be very few reads over individual isoforms/exons. Hard to say what's enough coverage to say difficult to detect an individual exon
# abvoe is conventional filter for gene-level stats analysis, not necessarily event-level
# 11 genes have GE below conventional filter

# Try some plots to visualise
liu_all_means %>%
  filter(abundance_type == "counts_norm") %>%
  mutate(gene_name = fct_reorder(gene_name, mean_all)) %>%
  ggplot(aes(x = log2(mean_all + 1), y = gene_name)) +
  geom_point(size = rel(0.1)) +
  geom_vline(xintercept = log2(11), linetype = "dashed") +
  geom_vline(xintercept = log2(26), linetype = "dashed") +
  geom_vline(xintercept = log2(50), linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(title = "Cryptic-containing gene expression") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
  
liu_all_means %>%
    filter(abundance_type == "counts_norm") %>%
    mutate(gene_name = fct_reorder(gene_name, mean_all)) %>%
    ggplot(aes(x = log2(mean_all + 1))) +
  stat_ecdf() +
    geom_vline(xintercept = log2(11), linetype = "dashed") +
    geom_vline(xintercept = log2(26), linetype = "dashed") +
    geom_vline(xintercept = log2(51), linetype = "dashed") +
    theme_bw(base_size = 14) +
    labs(title = "Gene expression of cryptic-containing genes")
  

### Plots of enriched vs not enriched cryptic-containing genes

cryp_enriched_gn <- ppau_delta_paired_median_all_cryp %>%
  filter(median_paired_delta_ppau > 0.05) %>%
  pull(gene_name)

# label genes as passing enrichment criteria or not
liu_all_means <- liu_all_means %>%
  mutate(enr_cryp = if_else(gene_name %in% cryp_enriched_gn, T, F))

# for each mean calculation (all samples, -ve only, +ve only)
# make base dot plot highlighting expression of enriched cryptic-containing genes
base_dot_plots <- colnames(liu_all_means)[str_starts(colnames(liu_all_means), "mean_")] %>%
  set_names() %>%
  map(~ base_dotplot(filter(liu_all_means, abundance_type == "counts_norm"), .x, glue("{.x} - Cryptic gene expression")))

# add gene names to all plots
gn_dot_plots <- base_dot_plots %>%
  map(~ .x +
        geom_text_repel(size = rel(2.5),
                        direction = "x",
                        min.segment.length = 0,
                        force = 25, 
                        max.overlaps = 100, 
                        segment.size = 0.1, 
                        segment.linetype = 2,
                        seed = 123)
      )

# write to file
if (!dir.exists("processed/gene_exprn")) {dir.create("processed/gene_exprn", recursive = T)}

walk2(base_dot_plots,
      names(base_dot_plots),
      ~ ggsave(filename = glue("processed/gene_exprn/2023-09-13_cryptic_gene_exprn_dot_plot.{.y}.no_label.png"),
               plot = .x,
               height = 8,
               width = 8,
               units = "in",
               dpi = "retina")
      )

walk2(gn_dot_plots,
      names(gn_dot_plots),
      ~ ggsave(filename = glue("processed/gene_exprn/2023-09-13_cryptic_gene_exprn_dot_plot.{.y}.label.png"),
               plot = .x,
               height = 8,
               width = 8,
               units = "in",
               dpi = "retina")
)


### Summary plots of gene expression for enriched vs not enriched cryptics

cryp_enr_vs_ge_densitys <- colnames(liu_all_means)[str_starts(colnames(liu_all_means), "mean_")] %>%
  set_names() %>%
  map(~ ge_density_plot(filter(liu_all_means, abundance_type == "counts_norm"), .x, glue("{.x} - Cryptic gene expression")))

cryp_enr_vs_ge_boxplots <- colnames(liu_all_means)[str_starts(colnames(liu_all_means), "mean_")] %>%
  set_names() %>%
  map(~ ge_boxplot(filter(liu_all_means, abundance_type == "counts_norm"), .x, glue("{.x} - Cryptic gene expression")))

# write to file
walk2(cryp_enr_vs_ge_densitys,
      names(cryp_enr_vs_ge_densitys),
      ~ ggsave(filename = glue("processed/gene_exprn/2023-09-13_cryptic_gene_exprn_enr_vs_not_density_plot.{.y}.png"),
               plot = .x,
               height = 8,
               width = 8,
               units = "in",
               dpi = "retina")
)

walk2(cryp_enr_vs_ge_boxplots,
      names(cryp_enr_vs_ge_boxplots),
      ~ ggsave(filename = glue("processed/gene_exprn/2023-09-13_cryptic_gene_exprn_enr_vs_not_box_plot.{.y}.png"),
               plot = .x,
               height = 8,
               width = 8,
               units = "in",
               dpi = "retina")
)




#### PLAYGROUND
# Random thoughts:
# - Are 'enriched' genes over-represented in highly expressed genes? See where they pop up in curve, possibly histograms of expression
# - could make similar cumulative distribution plot (maybe just individual dots) coloured by event type as well


# liu_all_means %>%
#   filter(abundance_type == "counts_norm") %>%
#   mutate(gene_name = fct_reorder(gene_name, mean_all),
#          enr_cryp = if_else(gene_name %in% cryp_enriched_gn, T, F),
#          plot_gn = if_else(enr_cryp, gene_name, "")) %>%
#   ggplot(aes(x = log2(mean_all + 1), y = gene_name, colour = enr_cryp, label = plot_gn)) +
#   geom_point(size = rel(0.1)) +
#   geom_vline(xintercept = log2(11), linetype = "dashed", alpha = 0.5) +
#   geom_vline(xintercept = log2(26), linetype = "dashed", alpha = 0.5) +
#   geom_vline(xintercept = log2(51), linetype = "dashed", alpha = 0.5) +
#   geom_text_repel(size = rel(2),
#                   direction = "x",
#                   min.segment.length = 0,
#                   force = 25, max.overlaps = 100, segment.size = 0.1, segment.linetype = 2,
#                   seed = 123) +
#   scale_colour_manual(values = c("#bdbdbd", "#d95f02")) +
#   theme_classic(base_size = 14) +
#   labs(title = "Cryptic-containing gene expression",
#        subtitle = "lines left-right are mean of 10, 25 & 50 reads in gene",
#        colour = "Passes enrichment criteria") +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         legend.position = "top")
# 
# 
# liu_all_means %>%
#   filter(abundance_type == "counts_norm") %>%
#   mutate(gene_name = fct_reorder(gene_name, mean_all),
#          enr_cryp = if_else(gene_name %in% cryp_enriched_gn, T, F)) %>%
#   ggplot(aes(x = log2(mean_all + 1), y = gene_name, colour = enr_cryp)) +
#   geom_point(size = rel(0.1)) +
#   geom_vline(xintercept = log2(11), linetype = "dashed") +
#   geom_vline(xintercept = log2(26), linetype = "dashed") +
#   geom_vline(xintercept = log2(51), linetype = "dashed") +
#   theme_bw(base_size = 14) +
#   labs(title = "Cryptic-containing gene expression") +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())
# 
# 
# liu_all_means %>%
#   filter(abundance_type == "counts_norm") %>%
#   mutate(gene_name = fct_reorder(gene_name, mean_all),
#          enr_cryp = if_else(gene_name %in% cryp_enriched_gn, T, F)) %>%
#   ggplot(aes(x = log2(mean_all + 1), fill = enr_cryp, colour = enr_cryp)) +
#   geom_density(alpha = 0.5)
# 
# liu_all_means %>%
#   filter(abundance_type == "counts_norm") %>%
#   mutate(gene_name = fct_reorder(gene_name, mean_all),
#          enr_cryp = if_else(gene_name %in% cryp_enriched_gn, T, F)) %>%
#   ggplot(aes(y = log2(mean_all + 1), x = enr_cryp, colour = enr_cryp)) +
#   geom_boxplot(outlier.shape = NA) + 
#   geom_jitter(alpha = 0.5) +
#   ggpubr::stat_compare_means()
# 
# # atte,pt to plot gene name alongside cryptics
# liu_all_means %>%
#   filter(abundance_type == "counts_norm") %>%
#   mutate(gene_name = fct_reorder(gene_name, mean_all),
#          enr_cryp = if_else(gene_name %in% cryp_enriched_gn, T, F),
#          plot_gn = if_else(enr_cryp, gene_name, "")) %>%
#   ggplot(aes(x = log2(mean_all + 1), y = gene_name, colour = enr_cryp, label = plot_gn)) +
#   geom_point(size = rel(0.1)) +
#   geom_text(size = rel(2), nudge_x = 5) +
#   geom_vline(xintercept = log2(11), linetype = "dashed") +
#   geom_vline(xintercept = log2(26), linetype = "dashed") +
#   geom_vline(xintercept = log2(51), linetype = "dashed") +
#   scale_colour_manual(values = c("#bdbdbd", "#d95f02")) +
#   theme_classic(base_size = 14) +
#   labs(title = "Cryptic-containing gene expression") +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())
#   
# # plot gene name with ggrepel
#   liu_all_means %>%
#     filter(abundance_type == "counts_norm") %>%
#     mutate(gene_name = fct_reorder(gene_name, mean_all),
#            enr_cryp = if_else(gene_name %in% cryp_enriched_gn, T, F),
#            plot_gn = if_else(enr_cryp, gene_name, "")) %>%
#     ggplot(aes(x = log2(mean_all + 1), y = gene_name, colour = enr_cryp, label = plot_gn)) +
#     geom_point(size = rel(0.1)) +
#     geom_vline(xintercept = log2(11), linetype = "dashed", alpha = 0.5) +
#     geom_vline(xintercept = log2(26), linetype = "dashed", alpha = 0.5) +
#     geom_vline(xintercept = log2(51), linetype = "dashed", alpha = 0.5) +
#     geom_text_repel(size = rel(2),
#                     direction = "x",
#                     min.segment.length = 0,
#                     force = 25, max.overlaps = 100, segment.size = 0.1, segment.linetype = 2,
#                     seed = 123) +
#     scale_colour_manual(values = c("#bdbdbd", "#d95f02")) +
#     theme_classic(base_size = 14) +
#     labs(title = "Cryptic-containing gene expression",
#          subtitle = "lines left-right are mean of 10, 25 & 50 reads in gene",
#          colour = "Passes enrichment criteria") +
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "top")
