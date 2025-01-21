library(tidyverse)

#' Calculate nucleotide proportions across positions from a wide-format position-frequency matrix (output of get_position_nucleotide_frequency.py)
#' 
#' @description
#' Takes a position frequency matrix in tibble format and calculates nucleotide 
#' proportions across all positions. Optionally filters out 'N' nucleotides and 
#' returns data in long format ready for plotting.
#'
#' @param pfm_tibble A tibble containing nucleotide frequencies. Should have one column
#'   for nucleotides and remaining columns representing positions.
#' @param remove_Ns Logical indicating whether to remove 'N' nucleotides before
#'   calculating proportions. Default is TRUE.
#' @param return_long Logical indicating whether to return data in long format.
#'   Default is FALSE.
#' @param nucleotide_col Character string specifying the name of the nucleotide column.
#'   Default is "nucleotide".
#'
#' @return A tibble containing nucleotide proportions. If return_long is TRUE,
#'   returns data in long format with columns for nucleotide, position, and fraction.
#'   Otherwise returns data in wide format.
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' pfm <- tibble(
#'   nucleotide = c("A", "C", "G", "T", "N"),
#'   pos1 = c(10, 20, 15, 5, 2),
#'   pos2 = c(12, 18, 16, 4, 1)
#' )
#' calc_nucleotide_proportions(pfm)
#' calc_nucleotide_proportions(pfm, return_long = TRUE)
#'
#' @export
calc_nucleotide_proportions <- function(pfm, 
                                        remove_Ns = TRUE,
                                        return_long = TRUE,
                                        nucleotide_col = "nucleotide") {
  
  # Initial data processing
  pfm
  
  # Filter out Ns if requested
  if (remove_Ns) {
    pfm <- pfm %>%
      filter(!!sym(nucleotide_col) != "N")
  }
  
  # Calculate proportions
  result <- pfm %>%
    mutate(across(-!!sym(nucleotide_col), ~ .x / sum(.x)))
  
  # Convert to long format if requested
  if (return_long) {
    result <- result %>%
      pivot_longer(-!!sym(nucleotide_col), 
                   names_to = "position",
                   values_to = "fraction") %>%
      mutate(position = as.numeric(position))
  }
  
  return(result)
}


##

kd_patrs_pfm <- read_tsv("processed/curation/patr_internal_priming/nuc_freq.polya_clusters.extend_50_both.align_center.tsv")
kd_patrs_min3_pfm <- read_tsv("processed/curation/patr_internal_priming/nuc_freq.polya_clusters.extend_50_both.min_3_reads.align_center.tsv")
kd_patrs_min5_pfm <- read_tsv("processed/curation/patr_internal_priming/nuc_freq.polya_clusters.extend_50_both.min_5_reads.align_center.tsv")
polyadb_pfm <- read_tsv("processed/curation/patr_internal_priming/nuc_freq.polyadb_v3.extend_50_both.align_center.tsv")

# Across all positions, calculate nucleotide proportions
# do separately with and without Ns
pg_kd_pfm_prop_withn <- kd_patrs_pfm %>%
  mutate(across(-nucleotide, ~ .x / sum(.x)))

pg_kd_pfm_prop_withoutn <- kd_patrs_pfm %>%
  filter(nucleotide != "N") %>%
  mutate(across(-nucleotide, ~ .x / sum(.x)))

kd_pfm_prop_withn
kd_pfm_prop_withoutn

# as difference so small, continue without Ns
long_kd_pfm_prop_withoutn <- calc_nucleotide_proportions(kd_pfm_prop_withoutn)
long_kd_min3_pfm_prop_withoutn <- calc_nucleotide_proportions(kd_patrs_min3_pfm)
long_kd_min5_pfm_prop_withoutn <- calc_nucleotide_proportions(kd_patrs_min5_pfm)
long_polyadb_pfm_prop_withoutn <- calc_nucleotide_proportions(polyadb_pfm)
long_kd_pfm_prop_withoutn
long_polyadb_pfm_prop_withoutn

long_kd_pfm_prop_withoutn %>%
  ggplot(aes(x = position, y = fraction, colour = nucleotide, group = nucleotide)) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-50,50,5)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

long_polyadb_pfm_prop_withoutn %>%
  ggplot(aes(x = position, y = fraction, colour = nucleotide, group = nucleotide)) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-50,50,5)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

# try some smoothing
long_kd_pfm_prop_withoutn %>%
  ggplot(aes(x = position, y = fraction, colour = nucleotide, group = nucleotide)) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-50,50,5)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")


# combined plot
bind_rows(patrs = long_kd_pfm_prop_withoutn, polyadb = long_polyadb_pfm_prop_withoutn, .id = "origin") %>%
  ggplot(aes(x = position, y = fraction, colour = nucleotide, group = nucleotide)) +
  facet_wrap(~ origin, ncol = 2) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-50,50,5)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

# repeat filtering out those with only a single read
bind_rows(patrs_min3 = long_kd_min3_pfm_prop_withoutn, polyadb = long_polyadb_pfm_prop_withoutn, .id = "origin") %>%
  ggplot(aes(x = position, y = fraction, colour = nucleotide, group = nucleotide)) +
  facet_wrap(~ origin, ncol = 2) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-50,50,5)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

# slightly more stringent
bind_rows(patrs_min5 = long_kd_min5_pfm_prop_withoutn, polyadb = long_polyadb_pfm_prop_withoutn, .id = "origin") %>%
  ggplot(aes(x = position, y = fraction, colour = nucleotide, group = nucleotide)) +
  facet_wrap(~ origin, ncol = 2) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-50,50,5)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

