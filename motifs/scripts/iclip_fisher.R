library(tidyverse)
library(broom)

#' Perform fisher's exact test on 'summary count' df output by scripts/iclip_summarise_coverage.R
counts_to_fisher <- function(data) {
  
  # Split data by region_type into a list of tibbles
  data_grouped <- data %>%
    group_by(region_type) 
  
  region_groups <- data_grouped %>%
    group_split(.keep = F)
  
  # name the list elements according to the group key
  region_names <- data_grouped %>%
    group_keys() %>%
    pull(region_type)
  
  region_groups <- set_names(region_groups, region_names)
  
  # For each group/tibble, transform to 2x2 contigency matrix and perform fisher's exact test
  # Combine across groups into final df of Fisher's exact test results
  results <- purrr::map(region_groups, ~ {
      .x %>%
        pivot_wider(names_from = binding, values_from = n) %>%
        arrange(desc(cryptic_status)) %>%
        relocate(`1`, .before = `0`) %>%
        # select(-region_type) %>%
        column_to_rownames("cryptic_status") %>%
        fisher.test() %>%
        broom::tidy()
    }) %>%
    bind_rows(.id = "region_type")
  
  return(results)
}


outdir <- "processed/iclip_maps/coverage/background_shsy5y/summary"
spliced_counts <- read_tsv(file.path(outdir, "2024-11-13_counts.any_binding.window_500.spliced.tsv"))
d3utr_counts <- read_tsv(file.path(outdir, "2024-11-13_counts.any_binding.window_500.d3utr.tsv"))

fisher_d3utr <- counts_to_fisher(d3utr_counts)
fisher_spliced <- counts_to_fisher(spliced_counts)

write_tsv(fisher_d3utr, file.path(outdir, "2024-11-13_fisher_exact.any_binding.window_500.d3utr.tsv"))
write_tsv(fisher_spliced, file.path(outdir, "2024-11-13_fisher_exact.any_binding.window_500.spliced.tsv"))