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
seddighi_df <- clean_deseq_df(seddighi_df)

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
  theme_bw(base_size = 16) +
  scale_x_continuous(limits = c(-5,5),
                     breaks = seq(-10,10,1)) +
  guides(alpha = "none", colour = "none") +
  labs(
       x = "Log2FoldChange (KD / WT)",
       y = "-log10(padj)")

ggsave("processed/2023-09-26_seddighi_rna_de_volcano_cryptic_3utrs_lab.png",
       device = "png",
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")

ggsave("processed/2023-09-26_seddighi_rna_de_volcano_cryptic_3utrs_lab.svg",
       device = svg,
       height = 8,
       width = 8,
       dpi = "retina",
       units = "in")
