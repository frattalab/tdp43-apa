library(tidyverse)

parse_coverage <- function(file, flank_interval) {
  f <- read_tsv(file, col_names = c("chr", "start", "end", "ID", ".", "strand", "position", "coverage"))
  f$position[f$strand == "-"] <- (2 + flank_interval*2) - f$position[f$strand == "-"]
  
  average_coverage <- f %>%
    group_by(position) %>%
    summarize(avg_coverage = mean(coverage),
              se = sd(coverage) / sqrt(n()),
              # 95 % confidence intervals
              lwr = avg_coverage - 1.96*se, 
              upr = avg_coverage + 1.96*se,
              n_overlaps = sum(coverage),
              n_events = n(),
              frac_overlaps = sum(coverage) / n()
    )
  average_coverage
  
}

#