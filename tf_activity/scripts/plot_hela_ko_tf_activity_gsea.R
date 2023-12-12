library(tidyverse)
library(patchwork)
library(fgsea)
load("processed/ferguson_hela/2023-11-29_hela_ko_tf_activity_gsea.Rdata")


find_cross_0_idx <- function(x) {
  
  # Find the index of the last positive value
  last_positive_index <- tail(which(x > 0), 1)
  
  # Find the index of the first negative value
  first_negative_index <- head(which(x < 0), 1)
  
  # report value halfway between two indexes
  mean(c(last_positive_index, first_negative_index))
}

# find_cross_0_idx(c(10, 8, 6, 3, 2, 0, -1, -3, -5, -7)) # should return 6 (0)
# find_cross_0_idx(c(5, 3, 2, -1, -3, -5)) # should return 3.5 (2 = 3rd, -1 = 4th idx)

plot_gsea_line <- function(gene_list, ranks, nes, padj, plot_title = "", round_digits = 6, zero_line = FALSE) {

  # base plot   
  p <- plotEnrichment(gene_list,
                 ranks) +
    labs(title = plot_title,
         subtitle = glue::glue("NES = {round(nes, round_digits)}, padj = {round(padj, round_digits)}"),
         x = "Gene Ranks",
         y = "Enrichment Score") +
    theme(title = element_text(size = rel(1.5)),
          axis.text = element_text(size = rel(1.5))
    )
  
  if (zero_line) {
    # report point that is mean of ranks of last positive and first negative values
    cross_0_idx <- find_cross_0_idx(ranks)
    p <- p + geom_vline(xintercept = cross_0_idx, linetype = "dashed")
  }
  
  p
}

# make plots for just chip-seq

# results df of just chipseq
gsea_ferguson_chipseq <- gsea_ferguson_all %>%
  filter(str_starts(pathway, "chipseq"))
# target list of just chipseq
elk1_target_list_chipseq <- elk1_target_list["chipseq_hela_both_elk1"]

# want to make separate plots per stat
# so split into group
tmp_grpd <- gsea_ferguson_chipseq %>%
  mutate(plot_pathway = "ELK1 HeLA ChIP-seq") %>%
  group_by(plot_pathway, score_type)

gsea_pway_grpd <- group_split(tmp_grpd) %>%
  set_names(group_keys(tmp_grpd) %>% mutate(id = paste(plot_pathway, score_type, sep = " - ")) %>% pull())  # group keys returns 1 col df in case of single group

# create classic gsea enrichment plots for each stat on HELA KO
gsea_enrichplots_ferguson_chipseq <- pmap(list(x = rev(ferguson_deseq_ranks_nospl), # here stat comes first, in other two pval is first
                              y = gsea_pway_grpd,
                              z = names(gsea_pway_grpd)),
                         function(x, y, z) plot_gsea_line(elk1_target_list_chipseq$chipseq_hela_both_elk1,
                                                          x, y$NES, y$padj, plot_title = z, zero_line = T)
)

# repeat for no diff spliced
# results df of just chipseq
gsea_ferguson_nospl_chipseq <- gsea_ferguson_nospl_all %>%
  filter(str_starts(pathway, "chipseq_hela_both"))

# want to make separate plots per stat
# so split into group
tmp_grpd <- gsea_ferguson_nospl_chipseq %>%
  mutate(plot_pathway = "ELK1 HeLA ChIP-seq") %>%
  group_by(plot_pathway, score_type)

gsea_pway_grpd <- group_split(tmp_grpd) %>%
  set_names(group_keys(tmp_grpd) %>% mutate(id = paste(plot_pathway, score_type, sep = " - ")) %>% pull())


gsea_enrichplots_ferguson_chipseq_nopsl <- pmap(list(x = rev(ferguson_deseq_ranks), # here stat comes first, in other two pval is first
                                                    y = gsea_pway_grpd,
                                                    z = names(gsea_pway_grpd)),
                                               function(x, y, z) plot_gsea_line(elk1_target_list_chipseq$chipseq_hela_both_elk1,
                                                                                x, y$NES, y$padj, plot_title = paste(z, "nospl", sep = " - "),
                                                                                zero_line = T)
                                               )

# want to make separate plots per stat
# so split into group
# finally for spliced removed, should make a function...
tmp_grpd <- gsea_ferguson_all %>%
  mutate(plot_pathway = pathway) %>%
  arrange(plot_pathway) %>%
  group_by(plot_pathway, score_type)

gsea_pway_grpd <- group_split(tmp_grpd) %>%
  set_names(group_keys(tmp_grpd) %>% mutate(id = paste(plot_pathway, score_type, sep = " - ")) %>% pull())

# need 1 set of target lists for each statistic
# first match up target lists with those present in GSEA results
gsea_targets_bool <- !names(elk1_target_list) %in% setdiff(names(elk1_target_list), unique(gsea_ferguson_all$pathway))


tmp_trgt_list <- rep(elk1_target_list[gsea_targets_bool],
                     2)
# make sure sort oreder matches the gsea pathway order
tmp_trgt_list <- tmp_trgt_list[group_keys(tmp_grpd) %>% pull(plot_pathway)]

# names(gsea_pway_grpd)
# names(tmp_trgt_list)
 
# for (i in seq(1,24)) {
#   
#   print(paste(names(gsea_pway_grpd)[i], names(tmp_trgt_list)[i], sep = " | "))
# 
# }

gsea_enrichplots_ferguson_all <- pmap(list(x = rep(ferguson_deseq_ranks, 12),
                                                 y = gsea_pway_grpd,
                                                 z = names(gsea_pway_grpd),
                                                 t = tmp_trgt_list),
                                            function(x, y, z, t) plot_gsea_line(t,
                                                                                x, y$NES, y$padj, 
                                                                                plot_title = z,
                                                                                zero_line = T)
) %>%
  set_names(names(gsea_pway_grpd))

# gsea_enrichplots_ferguson_all

# finally for spliced remvoed, should make a function...
tmp_grpd <- gsea_ferguson_nospl_all %>%
  mutate(plot_pathway = pathway) %>%
  arrange(plot_pathway) %>%
  group_by(plot_pathway, score_type)

gsea_pway_grpd <- group_split(tmp_grpd) %>%
  set_names(group_keys(tmp_grpd) %>% mutate(id = paste(plot_pathway, score_type, sep = " - ")) %>% pull())

names(gsea_pway_grpd)

# need 1 set of target lists for each statistic
gsea_targets_bool <- !names(elk1_target_list) %in% setdiff(names(elk1_target_list), unique(gsea_ferguson_nospl_all$pathway))

tmp_trgt_list <- rep(elk1_target_list[gsea_targets_bool],
                     2)
# make sure sort oreder matches the gsea pathway order
tmp_trgt_list <- tmp_trgt_list[group_keys(tmp_grpd) %>% pull(plot_pathway)]


gsea_enrichplots_ferguson_all_nospl <- pmap(list(x = rep(ferguson_deseq_ranks, 12),
                                           y = gsea_pway_grpd,
                                           z = names(gsea_pway_grpd),
                                           t = tmp_trgt_list),
                                      function(x, y, z, t) plot_gsea_line(t,
                                                                       x, y$NES, y$padj, plot_title = paste(z, "nospl", sep = " - "),
                                                                       zero_line = T)
                                      ) %>%
  set_names(names(gsea_pway_grpd))

# gsea_enrichplots_ferguson_chipseq
# gsea_enrichplots_ferguson_chipseq_nopsl
# gsea_enrichplots_ferguson_all
# gsea_enrichplots_ferguson_all_nospl

gsea_enrichplots_ferguson_comb_stat <- gsea_enrichplots_ferguson_chipseq$stat / gsea_enrichplots_ferguson_chipseq_nopsl$stat
gsea_enrichplots_ferguson_comb_signedp <- gsea_enrichplots_ferguson_chipseq$signed_pvalue / gsea_enrichplots_ferguson_chipseq_nopsl$signed_pvalue

gsea_enrichplots_ferguson_comb_all <- map2(.x = gsea_enrichplots_ferguson_all,
                                           .y = gsea_enrichplots_ferguson_all_nospl,
                                           ~ .x / .y)

# gsea_enrichplots_ferguson_comb_stat
# gsea_enrichplots_ferguson_comb_signedp
# gsea_enrichplots_ferguson_comb_all

if (!dir.exists("processed/ferguson_hela/svgs")) {dir.create("processed/ferguson_hela/svgs", recursive = T)}

ggsave(filename = "2023-12-12_ferguson_hela_chipseq_gsea_enrichplot_spl_nospl_stat.png",
       plot = gsea_enrichplots_ferguson_comb_stat,
       device = "png",
       path = "processed/ferguson_hela/",
       height = 8,
       width = 12,
       units = "in", 
       dpi = "retina")


ggsave(filename = "2023-12-12_ferguson_hela_chipseq_gsea_enrichplot_spl_nospl_signedp.png",
       plot = gsea_enrichplots_ferguson_comb_signedp,
       device = "png",
       path = "processed/ferguson_hela/",
       height = 8,
       width = 12,
       units = "in", 
       dpi = "retina")

# for others, plot in multipage pdfs


# Set the A4 paper size dimensions in inches
a4_width <- 8.27
a4_height <- 11.69

pdf("processed/ferguson_hela/2023-12-12_ferguson_hela_all_targets_gsea_enrichplot_spl_nospl.pdf",
    width = a4_height, height = a4_width )

# Loop through each ggplot object and print it to the PDF file
for (i in seq_along(gsea_enrichplots_ferguson_comb_all)) {
  print(gsea_enrichplots_ferguson_comb_all[[i]])
}

# Close the PDF file
dev.off()

# construct SVGs for each plot
# first for spl included only
walk2(.x = gsea_enrichplots_ferguson_all,
      .y = str_replace_all(names(gsea_enrichplots_ferguson_all), " - ", "."),
      ~ ggsave(filename = paste("2023-12-12_gsea_enrichplot_spl", .y, "svg",sep = "."),
               plot = .x,
               device = svg,
               path = "processed/ferguson_hela/svgs",
               height = 5,
               width = 15,
               dpi = "retina"
               ),
      .progress = T
      )
     

# next for no spl removed
walk2(.x = gsea_enrichplots_ferguson_all_nospl,
      .y = str_replace_all(names(gsea_enrichplots_ferguson_all_nospl), " - ", "."),
      ~ ggsave(filename = paste("2023-12-12_gsea_enrichplot_nospl", .y, "svg",sep = "."),
               plot = .x,
               device = svg,
               path = "processed/ferguson_hela/svgs",
               height = 5,
               width = 15,
               dpi = "retina"
               ),
      .progress = T
      )

# next for both combined
walk2(.x = gsea_enrichplots_ferguson_comb_all,
      .y = str_replace_all(names(gsea_enrichplots_ferguson_comb_all), " - ", "."),
      ~ ggsave(filename = paste("2023-12-12_gsea_enrichplot_both", .y, "svg",sep = "."),
               plot = .x,
               device = svg,
               path = "processed/ferguson_hela/svgs",
               height = 7.5,
               width = 15,
               dpi = "retina"
      ),
      .progress = T
)
