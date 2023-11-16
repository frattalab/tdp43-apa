library(decoupleR)
library(readr)

collectri_hs <- get_collectri(organism='human', split_complexes=FALSE)

write_tsv(collectri_hs, "data/2023-11-15_collectri_homosapiens.tsv")
