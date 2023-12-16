library(decoupleR)
library(progeny)
library(readr)


progeny_full <- progeny::model_human_full
# decoupleR vignette uses top 500
progeny_top500 <- get_progeny(organism = 'human', top = 500)
# progeny vignette recommends top 100 for each pathway
progeny_top100 <- get_progeny(organism = 'human', top = 100)

if (!dir.exists("data")) {dir.create("data")}

write_tsv(progeny_full, "data/2023-12-16_progeny_full_homosapiens.tsv.gz",col_names = T)
write_tsv(progeny_top100, "data/2023-12-16_progeny_top100_homosapiens.tsv", col_names = T )
write_tsv(progeny_top500, "data/2023-12-16_progeny_top500_homosapiens.tsv", col_names = T)