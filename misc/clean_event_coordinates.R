library(tidyverse)
library(plyranges)

# Quick script - read in quantification GTF, collapse to unique intervals for every le_id and output to file

papa_quant_gtf <- read_gff2("data/novel_ref_combined.quant.last_exons.gtf")
outdir <- "processed"

# summarise coordinates to single string
le_id_coords <-  papa_quant_gtf %>%
  as_tibble() %>%
  select(le_id, start, end) %>%
  distinct(.keep_all = T) %>%
  arrange(le_id, start, end) %>%
  group_by(le_id) %>%
  summarise(start = paste(start, collapse = ";"),
            end = paste(end, collapse = ";")
            )

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}
write_tsv(le_id_coords, file.path(outdir, "le_id_collapsed_coords.quant.last_exons.tsv"))
