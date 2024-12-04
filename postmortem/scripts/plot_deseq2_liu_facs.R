library(tidyverse)
library(ggrepel)

df <- read_tsv("processed/gene_exprn/deseq2/2023-20-11_deseq2_liu_facs_results.tsv")
outdir <- "processed/gene_exprn/deseq2"

# add ranks
df <- df %>%
  mutate(abs_fc_shrink = abs(log2FoldChangeShrink),
         rank_abs_fc_shrink = min_rank(desc(abs_fc_shrink)),
         rank_padj = min_rank(padj))


# prepare for volcano plot
genes_to_label <- c("ELK1")

# log10 padj
plot_df <- df %>%
  mutate(plot_padj = -log10(padj),
         plot_padj = if_else(plot_padj > 40, 40, plot_padj), # shrink extreme pvalues to reduce plotting window
         plot_label = if_else(gene_name %in% genes_to_label, gene_name, ""),
         plot_alpha = case_when(plot_label != "" ~ 1,
                                padj < 0.05 ~ 0.5,
                                T ~ 0.01
                                ),
         plot_size = if_else(plot_label != "", 1, 0),
         plot_colour = if_else(plot_label != "", T, F)
         )

# volcano with ELK1 labelled
volcano_elk1 <- plot_df %>%
  ggplot(aes(x = log2FoldChangeShrink,
             y = plot_padj,
             label = plot_label,
             alpha = plot_alpha,
             colour = plot_colour)) +
  geom_point(size = 2) +
  geom_text_repel(seed = 123, size = 6, min.segment.length = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(breaks = seq(-4, 4 , 1)) +
  scale_colour_manual(values = c("#bdbdbd", "#d95f02")) +
  theme_bw(base_size = 12) +
  labs(x = "RNA-seq Log2FoldChange\n(TDPnegative / TDPpositive)",
       y = "-log10(padj)") +  
  guides(alpha = "none", colour = "none", size = "none")

volcano_elk1

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

ggsave(filename = "2024-11-27_liu_facs_deseq2_volcano.elk1.png",
       plot = volcano_elk1,
       path = outdir,
       width = 100,
       height = 100,
       units = "mm",
       dpi = "retina")

