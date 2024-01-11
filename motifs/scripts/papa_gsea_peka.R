library(tidyverse)
library(fgsea)
source("../riboseq/helpers.R") # load in get_ranked_gene_list
set.seed(123)

dbrn_tbl <- read_tsv("processed/peka/papa/2023-11-27_papa_cryptics_kmer6_window_250_distal_window_500_relpos_0.cleaned_6mer_distribution_genome_simple.tsv")
# remove bleedtrhough events, too low n to be reliable
dbrn_tbl <- filter(dbrn_tbl, str_starts(comparison_name, "bleedthrough", negate = T))



# map comparison names to cleaned event types & region types
plot_clean_names <- tibble(comparison_name = c("bleedthrough_exonstart", "bleedthrough_pas", "spliced_exonstart", "spliced_pas", "d3utr_pas_proximal", "d3utr_pas_distal"),
                           plot_event_type = c("IPA", "IPA", "ALE", "ALE", "3'Ext", "3'Ext"),
                           plot_region_type = factor(c("Intron Start", "PAS", "Exon Start", "PAS", "Proximal", "Distal"),
                                                     levels = c("Intron Start", "Exon Start", "PAS", "Proximal", "Distal")
                           )
                           )


# subset to just SH-SY5Y background & add cleaned names
dbrn_tbl <- dbrn_tbl %>%
  filter(str_ends(comparison_name, "background_shsy5y$")) %>%#
  # remove background suffix
  mutate(comparison_name = str_remove(comparison_name, "\\.background_shsy5y$")) %>%
  # add in cleaned category names
  left_join(plot_clean_names, by = "comparison_name")

# calculate ranks of kmers for each comparison
dbrn_tbl_ranks <- dbrn_tbl %>%
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

# add a set of all 20 kmers combined
tdp_motif_groups[["combined-motifs"]] <- unlist(tdp_motif_groups,use.names = F)
tdp_motif_groups

# for each comparison, get a named vector sorted from highest to lowest PEKA score (most to least enriched kmer)
# NB: use standard (but 1 sided so as not to confound depleted kmers with TDP-43 binding
tmp_grpd <- dbrn_tbl %>%
  mutate(PEKA_score_abs = abs(PEKA_score)) %>%
  group_by(comparison_name)

comparison_kmer_ranks <- tmp_grpd %>%
  group_split() %>%
  set_names(group_keys(tmp_grpd) %>% pull()) %>%
  map(~ get_ranked_gene_list(drop_na(.x, PEKA_score), "PEKA_score", "kmer"))

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

## collapse leading edge gene names to a comma separated string
gsea_comparisons <- gsea_comparisons %>%
  # expand vector to single row per kmer
  unnest_longer(leadingEdge) %>%
  group_by(comparison_name, pathway) %>%
  summarise(leadingEdgeName = paste(leadingEdge, collapse = ",")) %>%
  ungroup() %>%
  # join back to original df
  left_join(gsea_comparisons, ., by = c("comparison_name", "pathway")) %>% 
  select(-leadingEdge) %>%
  rename(leadingEdge = leadingEdgeName)


# ties in PEKA score - at what values do they occur?
# todo: FIX this definition, see if can rank by different metric + rerun gsea with different approach
ties_abs_peka <- count(drop_na(tmp_grpd, PEKA_score), PEKA_score) %>% ungroup() %>% filter(n > 1)

# Are they over-represented in certain groups?
count(ties_abs_peka, comparison_name)

# Where are the values in the PEKA_score dbrn - tend to be low?

duplicated_kmer_ranks <- map(comparison_kmer_ranks, ~ enframe(.x, "kmer", "PEKA_score")) %>%
  bind_rows(.id = "comparison_name") %>%
  # add counts for each value
  left_join(ties_abs_peka, by = c("comparison_name", "PEKA_score")) %>%
  # add column to label duplicated score
  mutate(score_duplicated = !is.na(n))

duplicated_peka_density <- duplicated_kmer_ranks %>%
  ggplot(aes(x = PEKA_score, fill = score_duplicated)) +
  facet_wrap("~ comparison_name", scales = "free_y") +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#7570b3", "#1b9e77")) + 
  theme_bw(base_size = 20) + 
  labs(fill = "Duplicated PEKA score",
       x = "PEKA score",
       y = "Density") +
  theme(legend.position = "top")

duplicated_peka_density

# generally have similar normal-like shape to distribution, perhaps slightly narrower around 0
# Exception is distal 3'UTRs where sharp peak around 0, which is group that has considerably most duplicated values

duplicated_kmer_ranks %>%
  filter(score_duplicated) %>%
  left_join(select(dbrn_tbl, comparison_name, kmer, PEKA_score, etxn, artxn, aroxn), by = c("comparison_name", "kmer")) %>%
  mutate(na_extn = is.na(etxn),
         zero_artxn = artxn == 0,
         zero_aroxn = aroxn == 0) %>%
  count(comparison_name, na_extn, zero_artxn, zero_aroxn)

# All have some occurrence in background regions (aroxn)
# Many not found once in cryptic regions (artxn), bias strongest for distal pas
# Remaining are 'true duplicates' - some occurrence in both cryptic and bg, but same PEKA-score
# for 'true duplicates', check whether could rank differently e.g. 
# A tibble: 8 × 5
# comparison_name    na_extn zero_artxn zero_aroxn     n
# <chr>              <lgl>   <lgl>      <lgl>      <int>
#   1 d3utr_pas_distal   FALSE   FALSE      FALSE        446
# 2 d3utr_pas_distal   TRUE    TRUE       FALSE       1127
# 3 d3utr_pas_proximal FALSE   FALSE      FALSE         52
# 4 d3utr_pas_proximal TRUE    TRUE       FALSE         31
# 5 spliced_exonstart  FALSE   FALSE      FALSE         55
# 6 spliced_exonstart  TRUE    TRUE       FALSE          9
# 7 spliced_pas        FALSE   FALSE      FALSE         55
# 8 spliced_pas        TRUE    TRUE       FALSE         30

# could non-zero duplicated be ranked secondly by a different metric (e.g. relative enrichment, average occurrence in cryptic etc.?)
duplicated_kmer_ranks %>%
  filter(score_duplicated) %>%
  left_join(select(dbrn_tbl, comparison_name, kmer, PEKA_score, etxn, artxn, aroxn), by = c("comparison_name", "kmer")) %>%
  mutate(na_extn = is.na(etxn),
         zero_artxn = artxn == 0,
         zero_aroxn = aroxn == 0) %>%
  filter(!na_extn & !zero_artxn & !zero_aroxn) %>%
  arrange(comparison_name, PEKA_score.x)


# make a dot plot of the NES scores for each motif group
dotplot_gsea <- gsea_comparisons %>% 
  left_join(plot_clean_names, by = "comparison_name") %>%
  mutate(plot_pathway = if_else(pathway == "combined-motifs", true = "Combined motifs", false = pathway),
         plot_pathway = str_replace_all(plot_pathway, "-motifs$", " motifs"),
         plot_pathway = factor(plot_pathway, levels = c("Combined motifs", "YG-containing motifs", "YA-containing motifs", "AA-containing motifs")),
         plot_size = -log10(padj.all),
         sig = padj < 0.05) %>%
  ggplot(aes(x = NES, y = plot_pathway, colour = sig, size = plot_size * 5)) +
  facet_wrap("plot_event_type ~ plot_region_type", ncol = 2) +
  geom_point() +
  theme_bw(base_size = 20) +
  labs(
       x = "Normalised Enrichment Score",
       y = "",
       colour = "padj < 0.05",
       size = "-log10(padj)") +
  theme(legend.position = "top") +
  scale_size_continuous(labels = ~ .x / 5) +
  scale_colour_manual(values = c("#7570b3", "#1b9e77"))

dotplot_gsea

ggsave("2024-01-11_peka_cryptics_halleger_motif_groups_gsea_dotplot.png",
       plot = dotplot_gsea,
       path = "processed/peka/papa/plots/",
       device = "png",
       width = 12.5,
       height = 7.5,
       units = "in",
       dpi = "retina"
       )

ggsave("2024-01-11_peka_cryptics_halleger_motif_groups_gsea_dotplot.svg",
       plot = dotplot_gsea,
       path = "processed/peka/papa/plots/",
       device = svg,
       width = 12.5,
       height = 7.5,
       units = "in",
       dpi = "retina"
)


