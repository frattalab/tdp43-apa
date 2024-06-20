library(tidyverse)

#' Factory function to process a chunked 'coverage ratio' BED-like file produced by pipeline (to be used with readr::read_tsv_chunked())
process_chunk <- function(le_id_filter = NULL) {
  
  
  # https://stackoverflow.com/questions/49238163/how-to-pass-arguments-to-a-callback-function-for-readrread-csv-chunked
  # if want to parameterise chunked function, need this to be a function factory
  
  if (is.null(le_id_filter)) {
    
    return(function(df, pos) {
      mutate(df,
             le_id = str_split_i(name, "\\|", 1),
             up_down_ratio = if_else(up_down_ratio == "inf", Inf, as.numeric(up_down_ratio)),
             up_down_ratio = as.numeric(up_down_ratio)
             )
      })
    
  } else {
    
    return(function(df, pos) {
      mutate(df, le_id = str_split_i(name, "\\|", 1)) %>% 
        filter(le_id %in% le_id_filter) %>%
        mutate(up_down_ratio = if_else(up_down_ratio == "inf", Inf, as.numeric(up_down_ratio)),
               up_down_ratio = as.numeric(up_down_ratio)
               )
      })
    
  }
  

  
}


cryptics_df <- read_tsv("processed/2023-12-10_cryptics_summary_all_events.tsv")
cryptics_df_mv <- read_tsv("processed/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv")

outdir <- "processed/curation/coverage_ratios"
if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

# extract i3 cortical cryptics
cryptic_le_ids_i3 <- cryptics_df_mv %>%
  filter(str_detect(experiment_name, "i3_cortical")) %>%
  distinct(le_id) %>%
  pull()

cryptic_le_ids_sh <- cryptics_df_mv %>%
  filter(str_detect(experiment_name, "shsy5y")) %>%
  distinct(le_id) %>%
  pull()  

# read in coverage paths and extract sample names
i3_paths <- list.files(path = "data/pas_curation/window_coverages/seddighi_i3_cortical/coverage",
                       pattern = ".pas_windows.summarised_coverage.ratio.bed$",
                       full.names = T) %>%
  set_names(str_remove_all(basename(.), ".pas_windows.summarised_coverage.ratio.bed$"))

shsy5y_paths <- list.files(path = "data/pas_curation/window_coverages/zanovello_shsy5y_chx/coverage",
                           pattern = ".pas_windows.summarised_coverage.ratio.bed$",
                           full.names = T) %>%
  set_names(str_remove_all(basename(.), ".pas_windows.summarised_coverage.ratio.bed$"))


# read in all coverage dfs, filtering for cryptics only and combining into single df
cov_cryptics_i3 <- map(i3_paths,
    ~ read_tsv_chunked(.x,
                 callback = DataFrameCallback$new(process_chunk(le_id_filter = cryptic_le_ids_i3)),
                 col_names = c("chrom", "start", "end", "name", "score", "strand", "up_cov", "down_cov", "up_down_ratio"),
                 show_col_types = F,
                 progress = T
                 ),
    .progress = T
    ) %>%
  bind_rows(.id = "sample_name")

cov_cryptics_sh <- map(shsy5y_paths,
                       ~ read_tsv_chunked(.x,
                                          callback = DataFrameCallback$new(process_chunk(le_id_filter = cryptic_le_ids_sh)),
                                          col_names = c("chrom", "start", "end", "name", "score", "strand", "up_cov", "down_cov", "up_down_ratio"),
                                          show_col_types = F,
                                          progress = T
                       ),
                       .progress = T
) %>%
  bind_rows(.id = "sample_name")

id_cols <- str_split("chrom, start, end, name, score, strand, le_id", pattern = ", ")[[1]]

# Calculate summary statistics (mean, sd, median) across samples for every interval and coverage metric
cov_cryptics_i3_summ <- cov_cryptics_i3 %>%
  group_by(across(all_of(id_cols))) %>%
  summarise(across(.cols = all_of(c("up_cov", "down_cov", "up_down_ratio")),
                   .fns = list(mean = ~ mean(.x),
                               sd = ~ sd(.x),
                               median = ~ median(.x)
                               ),
                   .names = "{.fn}__{.col}"
                   )
            ) %>%
  ungroup() %>%
  # Label an event as originating from PAPA / putatively updated
  mutate(end_origin = if_else(str_ends(name, "original$"), "PAPA", "PATR"))


cov_cryptics_sh_summ <- cov_cryptics_sh %>%
  group_by(across(all_of(id_cols))) %>%
  summarise(across(.cols = all_of(c("up_cov", "down_cov", "up_down_ratio")),
                   .fns = list(mean = ~ mean(.x),
                               sd = ~ sd(.x),
                               median = ~ median(.x)
                   ),
                   .names = "{.fn}__{.col}"
  )
  ) %>%
  ungroup() %>%
  # Label an event as originating from PAPA / putatively updated
  mutate(end_origin = if_else(str_ends(name, "original$"), "PAPA", "PATR"))


write_tsv(cov_cryptics_i3_summ, file.path(outdir, "2024-06-20_summarise_coverage_ratios.seddighi_i3_cortical.cryptics.tsv"))
write_tsv(cov_cryptics_sh_summ, file.path(outdir, "2024-06-20_summarise_coverage_ratios.zanovello_shsy5y_chx.cryptics.tsv"))


# repeat for all events
cov_cryptics_i3_all_summ <- map(i3_paths,
                                ~ read_tsv_chunked(.x,
                                                   callback = DataFrameCallback$new(process_chunk(le_id_filter = NULL)),
                                                   col_names = c("chrom", "start", "end", "name", "score", "strand", "up_cov", "down_cov", "up_down_ratio"),
                                                   show_col_types = F,
                                                   progress = T
                                ),
                                .progress = T
) %>%
  bind_rows(.id = "sample_name") %>%
  group_by(across(all_of(id_cols))) %>%
  summarise(across(.cols = all_of(c("up_cov", "down_cov", "up_down_ratio")),
                   .fns = list(mean = ~ mean(.x),
                               sd = ~ sd(.x),
                               median = ~ median(.x)
                   ),
                   .names = "{.fn}__{.col}"
  )
  ) %>%
  ungroup() %>%
  # Label an event as originating from PAPA / putatively updated
  mutate(end_origin = if_else(str_ends(name, "original$"), "PAPA", "PATR"))

write_tsv(cov_cryptics_i3_all_summ, file.path(outdir, "2024-06-20_summarise_coverage_ratios.seddighi_i3_cortical.all.tsv"))

cov_cryptics_sh_all_summ <- map(shsy5y_paths,
                       ~ read_tsv_chunked(.x,
                                          callback = DataFrameCallback$new(process_chunk(le_id_filter = NULL)),
                                          col_names = c("chrom", "start", "end", "name", "score", "strand", "up_cov", "down_cov", "up_down_ratio"),
                                          show_col_types = F,
                                          progress = T
                       ),
                       .progress = T
) %>%
  bind_rows(.id = "sample_name") %>%
  group_by(across(all_of(id_cols))) %>%
  summarise(across(.cols = all_of(c("up_cov", "down_cov", "up_down_ratio")),
                   .fns = list(mean = ~ mean(.x),
                               sd = ~ sd(.x),
                               median = ~ median(.x)
                   ),
                   .names = "{.fn}__{.col}"
  )
  ) %>%
  ungroup() %>%
  # Label an event as originating from PAPA / putatively updated
  mutate(end_origin = if_else(str_ends(name, "original$"), "PAPA", "PATR"))

write_tsv(cov_cryptics_sh_all_summ, file.path(outdir, "2024-06-20_summarise_coverage_ratios.zanovello_shsy5y_chx.all.tsv"))
