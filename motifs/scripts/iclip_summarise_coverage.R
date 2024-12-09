library(tidyverse)
source("scripts/fncs_plot_iclip.R")

#####

#' Summarise iCLIP coverage information to binary column across each interval
#'
#' This function processes a single input path to compute binding coverage statistics. 
#' It parses the coverage data, groups it by `name`, and determines if any positions
#' have a coverage of 1 to consider as 'bound'
#'
#' @param input_path A character vector (length 1) representing the input path or data to be processed.
#' @param flank_interval Single length integer denoting the width of extension either side of the point features
#' @return A tibble summarizing the binding coverage for each group (name field in BED-like file), with columns:
#' \itemize{
#'   \item \code{name}: The grouping variable.
#'   \item \code{binding}: A logical value indicating whether any position within the group has a coverage of 1.
#' }
#'
#' @import dplyr
#' @export
summarise_binding <- function(input_path, flank_interval) {
  input_path %>%
    parse_coverage(., flank_interval, average = FALSE) %>%
    group_by(name) %>%
    summarise(binding = any(coverage == 1)) %>%
    ungroup()
}

#' Summarise iCLIP coverage information across series of files and combine the results
#' 
#' This function processes a single input path to compute binding coverage statistics. 
#' It parses the coverage data, groups it by `name`, and determines if any positions
#' have a coverage of 1.
#'
#' @param file_paths A named character vector representing the coverage files to summarise and merge
#' @return A tibble summarizing the binding coverage for each group, with columns:
#' \itemize{
#'   \item \code{name}: The grouping variable.
#'   \item \code{binding}: A logical value indicating whether any position within the group has a coverage of 1.
#'   \item \code{origin}: Origin file-path/tibble - corresponds to 'name' field of input vector
#' }
#'
#' @import dplyr
#' @export
map_process_binding <- function(file_paths, flank_interval) {

  if (is.null(names(file_paths)) | any(names(file_paths) == "")) {
    stop("Error: The input `file_paths` must be a named vector.")
  }
  
  file_paths %>%
    map(~ summarise_binding(.x, flank_interval)) %>% # Apply summarise_binding to each file
    bind_rows(.id = "origin")     # Combine results into a single dataframe
}


#' Separate 'coverage origin' column into individual units
#'
#' This function takes a dataframe and separates the `origin` column into four
#' new columns: `background_type`, `type`, `region_type`, and `cryptic`, 
#' based on a period (`.`) delimiter.
#'
#' @param df A dataframe containing an `origin` column to be separated.
#' @return A dataframe with the `origin` column split into four columns: 
#' `background_type`, `type`, `region_type`, and `cryptic`.
#' @examples
#' input_df <- data.frame(origin = c("bg1.type1.region1.crypt1", "bg2.type2.region2.crypt2"))
#' separate_origin(input_df)
#' 
#' @import dplyr
#' @import tidyr
#' @export
separate_origin <- function(df) {
  if (!"origin" %in% colnames(df)) {
    stop("Error: The input dataframe must contain a column named 'origin'.")
  }
  
  df %>%
    separate(
      origin,
      into = c("background_type", "event_type", "region_type", "cryptic_status"),
      sep = "\\.",
      remove = TRUE
    )
}


#' Standardise 'region_type' column across event types
#'
#' This function standardises the `region_type` column in a dataframe. 
#' Values of `region_type` that are "proximal" or "le_start" are replaced with "start", 
#' while all other values are replaced with "end".
#'
#' @param df A dataframe containing a `region_type` column.
#' @return A dataframe with the `region_type` column standardised.
#' @examples
#' input_df <- data.frame(region_type = c("proximal", "distal", "le_start", "other"))
#' standardise_region_type(input_df)
#' 
#' @import dplyr
#' @export
standardise_region_type <- function(df) {
  if (!"region_type" %in% colnames(df)) {
    stop("Error: The input dataframe must contain a column named 'region_type'.")
  }
  
  df %>%
    mutate(
      region_type = if_else(
        region_type %in% c("proximal", "le_start"),
        "start",
        "end"
      )
    )
}


#####

# Find coverage files for each event type and region
d3utr_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/d3utr", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 
spliced_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/spliced", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 
bleedthrough_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/bleedthrough_uniq", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 
d3utrprox_paths <- list.files(path = "processed/iclip_maps/coverage/background_shsy5y/d3utr_proximal", pattern = "\\.txt\\.gz$", recursive = T, full.names = T) 

# dir for output TSVs
outdir <- "processed/iclip_maps/coverage/background_shsy5y/summary"


# Extract info from file path to construct IDs
# in format background_type.event_type.region_type.cryptic/bg

# processed/iclip_maps/coverage/background_shsy5y/d3utr/distal/regions.flank_500.coverage.bg.txt.gz -> background_shsy5y.d3utr.distal.bg 
d3utr_nm <- paste(str_remove(str_replace_all(dirname(d3utr_paths), "\\/", "."), "processed.iclip_maps.coverage."),
                  str_extract(basename(d3utr_paths), pattern = "bg|cr"),
                  sep = "."
)

# processed/iclip_maps/coverage/background_shsy5y/d3utr/distal/regions.flank_500.coverage.bg.txt.gz -> background_shsy5y.d3utr.distal.bg 
d3utrprox_nm <- paste(str_remove(str_replace_all(dirname(d3utrprox_paths), "\\/", "."),
                                 "processed.iclip_maps.coverage."),
                      str_extract(basename(d3utrprox_paths), pattern = "bg|cr"),
                      sep = "."
)

# background_shsy5y/spliced/pas/regions.flank_500.coverage.bg.txt.gz" -> background_shsy5y.spliced.pas.bg
spliced_nm <- paste(str_remove(str_replace_all(dirname(spliced_paths), "\\/", "."), "processed.iclip_maps.coverage."),
                    str_extract(basename(spliced_paths), pattern = "bg|cr"),
                    sep = "."
)

bleedthrough_nm <- paste(str_remove(str_replace_all(dirname(bleedthrough_paths), "\\/", "."), "processed.iclip_maps.coverage."),
                         str_extract(basename(bleedthrough_paths), pattern = "bg|cr"),
                         sep = "."
)


d3utr_paths <- set_names(d3utr_paths, d3utr_nm)
d3utrprox_paths <- set_names(d3utrprox_paths, d3utrprox_nm)
spliced_paths <- set_names(spliced_paths, spliced_nm)
bleedthrough_paths <- set_names(bleedthrough_paths, bleedthrough_nm)

# Read in all dfs per event type, summarising binding info to yes/no in complete search window
d3utr_cov_summ <- map_process_binding(d3utr_paths, 500) %>%
  separate_origin()

d3utrprox_cov_summ <- map_process_binding(d3utrprox_paths, 500) %>%
  separate_origin()

spliced_cov_summ <- map_process_binding(spliced_paths, 500) %>%
  separate_origin()

bleedthrough_cov_summ <- map_process_binding(bleedthrough_paths, 500) %>%
  separate_origin()


# summarise to wide format with standardised names across events

# standardise region type columns across event types
# all event types have 2 positions - 'start' (e.g. le_start, prox pas) & 'end' (pas, distal pas)
d3utr_cov_summ <- standardise_region_type(d3utr_cov_summ) %>%
  arrange(desc(cryptic_status), name)
d3utrprox_cov_summ <- standardise_region_type(d3utrprox_cov_summ) %>%
  arrange(desc(cryptic_status), name)
spliced_cov_summ <- standardise_region_type(spliced_cov_summ) %>%
  arrange(desc(cryptic_status), name)
bleedthrough_cov_summ <- standardise_region_type(bleedthrough_cov_summ) %>%
  arrange(desc(cryptic_status), name)

# Generate summary counts over each region type and cryptic status
d3utr_cov_summ_counts <- d3utr_cov_summ %>%
  count(region_type, cryptic_status, binding)

d3utrprox_cov_summ_counts <- d3utrprox_cov_summ %>%
  count(region_type, cryptic_status, binding)

spliced_cov_summ_counts <- spliced_cov_summ %>%
  count(region_type, cryptic_status, binding)

bleedthrough_cov_summ_counts <- bleedthrough_cov_summ %>%
  count(region_type, cryptic_status, binding)

# output TSVs to file

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}

# Output per-event summaries

write_tsv(mutate(d3utr_cov_summ, binding = as.integer(binding)),
          file = file.path(outdir, "2024-11-13_event.any_binding.window_500.d3utr.tsv"),
          col_names = T)

write_tsv(mutate(d3utrprox_cov_summ, binding = as.integer(binding)),
          file = file.path(outdir, "2024-11-13_event.any_binding.window_500.d3utr_proximal.tsv"),
          col_names = T)

write_tsv(mutate(spliced_cov_summ, binding = as.integer(binding)),
          file = file.path(outdir, "2024-11-13_event.any_binding.window_500.spliced.tsv"),
          col_names = T)

write_tsv(mutate(bleedthrough_cov_summ, binding = as.integer(binding)),
          file = file.path(outdir, "2024-11-13_event.any_binding.window_500.bleedthrough_uniq.tsv"),
          col_names = T)


# summary counts
write_tsv(mutate(d3utr_cov_summ_counts, binding = as.integer(binding)),
          file = file.path(outdir, "2024-11-13_counts.any_binding.window_500.d3utr.tsv"),
          col_names = T)

write_tsv(mutate(d3utrprox_cov_summ_counts, binding = as.integer(binding)),
          file = file.path(outdir, "2024-11-13_counts.any_binding.window_500.d3utr_proximal.tsv"),
          col_names = T)

write_tsv(mutate(spliced_cov_summ_counts, binding = as.integer(binding)),
          file = file.path(outdir, "2024-11-13_counts.any_binding.window_500.spliced.tsv"),
          col_names = T)

write_tsv(mutate(bleedthrough_cov_summ_counts, binding = as.integer(binding)),
          file = file.path(outdir, "2024-11-13_counts.any_binding.window_500.bleedthrough_uniq.tsv"),
          col_names = T)
