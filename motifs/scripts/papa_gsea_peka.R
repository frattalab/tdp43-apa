library(tidyverse)
library(fgsea)
source("../riboseq/helpers.R") # load in get_ranked_gene_list
set.seed(123)

dbrn_tbl_ale <- read_tsv("processed/peka/papa/2023-11-27_papa_cryptics_kmer6_window_250_distal_window_500_relpos_0.cleaned_6mer_distribution_genome_simple.tsv")
dbrn_tbl_3utr <- read_tsv("processed/peka/papa/2023-11-16_papa_cryptics_kmer6_window_250_distal_window_500_relpos_0.cleaned_6mer_distribution_genome_simple.tsv")

# map comparison names to cleaned event types & region types
plot_clean_names <- tibble(comparison_name = c("bleedthrough_exonstart", "bleedthrough_pas", "spliced_exonstart", "spliced_pas", "d3utr_pas_proximal", "d3utr_pas_distal"),
                           plot_event_type = c("Bleedthrough-ALE", "Bleedthrough-ALE", "AS-ALE", "AS-ALE", "3'UTR-ALE", "3'UTR-ALE"),
                           plot_region_type = factor(c("Intron Start", "PAS", "Exon Start", "PAS", "Proximal", "Distal"),
                                                     levels = c("Intron Start", "Exon Start", "PAS", "Proximal", "Distal")
                           )
                           )


# subset to just SH-SY5Y background
dbrn_tbl_ale <- dbrn_tbl_ale %>%
  filter(str_ends(comparison_name, "background_shsy5y$")) %>%#
  # remove background suffix
  mutate(comparison_name = str_remove(comparison_name, "\\.background_shsy5y$")) %>%
  # add in cleaned category names
  left_join(plot_clean_names, by = "comparison_name")

# subset to just 3'UTR-ALE, add in cleaned names
dbrn_tbl_3utr <- dbrn_tbl_3utr %>%
  filter(str_starts(comparison_name, "^d3utr")) %>%
  left_join(plot_clean_names, by = "comparison_name")


dbrn_tbl_comb <- bind_rows(dbrn_tbl_3utr, dbrn_tbl_ale)

# calculate ranks of kmers for each comparison
dbrn_tbl_ranks <- dbrn_tbl_comb %>%
  group_by(comparison_name) %>%
  arrange(desc(abs(PEKA_score)), .by_group = T) %>%
  mutate(peka_rank = row_number(),
         peka_rank_frac = peka_rank / n()) %>%
  ungroup()

# 6mer groups enriched in Halleger et al. TDP-43 + CTD MUT iCLIP
# Create list ready for use with fgsea
tdp_motif_groups <- list("YG-containing-motifs" = c("UGUGUG", "GUGUGU","UGUGCG", "UGCGUG","CGUGUG","GUGUGC"),
                         "YA-containing-motifs" = c("AUGUGU", "GUAUGU", "GUGUAU", "UGUGUA", "UGUAUG", "UGCAUG"),
                         "AA-containing-motifs" = c("GUGUGA", "AAUGAA", "GAAUGA", "UGAAUG", "AUGAAU", "GUGAAU", "GAAUGU", "UUGAAU")
                         )

tdp_motif_groups[["combined-motifs"]] <- unlist(tdp_motif_groups,use.names = F)
tdp_motif_groups

# for each comparison, get a named vector sorted from highest to lowest absolute PEKA score (most to least enriched kmer)

tmp_grpd <- dbrn_tbl_comb %>%
  mutate(PEKA_score_abs = abs(PEKA_score)) %>%
  group_by(comparison_name)

comparison_kmer_ranks <- tmp_grpd %>%
  group_split() %>%
  set_names(group_keys(tmp_grpd) %>% pull()) %>%
  map(~ get_ranked_gene_list(drop_na(.x, PEKA_score_abs), "PEKA_score_abs", "kmer"))

# also make ranks just using std PEKA scroe
# comparison_kmer_ranks_std_peka <- tmp_grpd %>%
#   group_split() %>%
#   set_names(group_keys(tmp_grpd) %>% pull()) %>%
#   map(~ get_ranked_gene_list(drop_na(.x, PEKA_score), "PEKA_score", "kmer"))

# run fgsea on kmer ranks of each comparison, using tdp motif groups as target list
# Q - are TDP-43 kmers over-represented at the top of the kmer list with respect to random chance?
# will use a 1 sided test here, as only interested if appears in most enriched kmers
gsea_comparisons <- map(comparison_kmer_ranks,
                        ~ fgsea(pathways = tdp_motif_groups,
                                stats = .x,
                                scoreType = "pos")
                        ) %>%
  bind_rows(.id = "comparison_name")

# gsea_comparisons_std_peka <- map(comparison_kmer_ranks_std_peka,
#                         ~ fgsea(pathways = tdp_motif_groups,
#                                 stats = .x,
#                                 scoreType = "pos")
# ) %>%
#   bind_rows(.id = "comparison_name")

# padj across all comparisons
gsea_comparisons <- mutate(gsea_comparisons, padj.all = p.adjust(pval, method = "BH"))


# ties in abs PEKA score - check where occur
ties_abs_peka <- count(tmp_grpd, PEKA_score_abs) %>% ungroup() %>% filter(n > 1)


# make a dot plot of the NES scores for each motif group
dotplot_gsea <- gsea_comparisons %>% 
  left_join(plot_clean_names, by = "comparison_name") %>%
  mutate(plot_pathway = factor(pathway, levels = c("combined-motifs", "YG-containing-motifs", "YA-containing-motifs", "AA-containing-motifs")),
         plot_size = -log10(padj.all),
         sig = padj < 0.05) %>%
  ggplot(aes(x = NES, y = plot_pathway, colour = sig, size = plot_size * 2)) +
  facet_wrap("plot_event_type ~ plot_region_type", ncol = 2) +
  geom_point() +
  theme_bw(base_size = 20) +
  labs(title = "Halleger et al TDP-43 motif groups + cryptic regions",
       x = "GSEA normalised enrichment score",
       y = "",
       colour = "padj < 0.05",
       size = "-log10(padj)") +
  theme(legend.position = "top") +
  scale_size_continuous(labels = ~ .x / 2)

dotplot_gsea

ggsave("2023-11-28_peka_cryptics_halleger_motif_groups_gsea_dotplot.png",
       plot = dotplot_gsea,
       path = "processed/peka/papa/plots/",
       device = "png",
       width = 12,
       height = 10,
       units = "in", dpi = "retina"
       )
