library(tidyverse)

parse_coverage <- function(file, flank_interval) {
  f <- read_tsv(file, col_names = c("chr", "start", "end", "ID", ".", "strand", "position", "coverage"))
  # make sure positions are strand-aware
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

plot_coverage <- function(df, ci_se_mult = 1.96, event_col = "plot_type", group_col = "cryptic", facet_ncol = 2, fill_colours = c("#000000", "#d95f02"), line_colours = c("#000000", "#d95f02")) {
  
  group_cols <- c(event_col, group_col)
  
  # generate confidence interval values
  plot_df <- df %>%
    mutate(plot_ymin = avg_coverage - (ci_se_mult*se),
           plot_ymax = avg_coverage + (ci_se_mult*se)) %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(ymin_smooth = stats::predict(loess(plot_ymin~position,span=0.1)),
           ymax_smooth = stats::predict(loess(plot_ymax~position,span=0.1))) %>%
    ungroup()
  
  plot_df %>%
    ggplot(aes(x = position, y = avg_coverage, color=!!sym(group_col), fill=!!sym(group_col))) +
    geom_smooth(method = "loess", span=0.2, se = F) +
    geom_ribbon(aes(ymin = ymin_smooth,  ymax = ymax_smooth, fill = !!sym(group_col)),
                alpha = 0.31) +
    xlab("Position") +
    ylab("Average Coverage") +
    geom_vline(xintercept = 500, linetype = "dashed", alpha = 0.5) +
    facet_wrap(event_col, ncol = facet_ncol, scales = "fixed") +
    scale_x_continuous(
      limits = c(0, 1001),
      breaks = seq(0,1000,100),
      labels = as.character(seq(-500,500,100))
    ) +
    scale_fill_manual(values = fill_colours) +
    scale_color_manual(values = line_colours) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top")
  
  
}

d3utr_paths <- list.files(path = "processed/iclip_maps/d3utr", full.names = T) 

# split by '.', extract the two elements before '.txt.gz'
d3utr_nm <- str_split(basename(d3utr_paths), "\\.",simplify = T) %>% 
  apply(MARGIN = 1,
        function(x) {paste(c(x[length(x)-3], x[length(x)-2]),
                           collapse = ".")
          })

# read in coverages
d3utr_average_coverage <- d3utr_paths %>%
  set_names(d3utr_nm) %>%
  map(~ parse_coverage(.x, 500)) %>%
  bind_rows(.id = "origin") %>%
  # pull out site type & status
  separate(origin, into = c("type", "cryptic"), sep = "\\.", remove = T)

# make df plot ready
d3utr_average_coverage <- d3utr_average_coverage %>%
  mutate(plot_type = if_else(type == "prox", "Proximal",
                             "Distal"),
         plot_type = factor(plot_type, levels = c("Proximal", "Distal"))
         )



plot_coverage(d3utr_average_coverage)
plot_coverage(d3utr_average_coverage, ci_se_mult = 1, facet_ncol = 1)

# average_coverage_prox_bg <- parse_coverage(prox_bg_file)
# average_coverage_dist_bg <- parse_coverage(dist_bg_file)
# average_coverage_prox_cr <- parse_coverage(prox_cr_file)
# average_coverage_dist_cr <- parse_coverage(dist_cr_file)
# 
# average_coverage_prox_cr$cryptic <- "cryptic"
# average_coverage_dist_cr$cryptic <- "cryptic"
# average_coverage_prox_bg$cryptic <- "background"
# average_coverage_dist_bg$cryptic <- "background"
# 
# average_coverage_prox_cr$type <- "proximal"
# average_coverage_dist_cr$type <- "distal"
# average_coverage_prox_bg$type <- "proximal"
# average_coverage_dist_bg$type <- "distal"


### Playground

# 
# d3utr_average_coverage %>%
#   ggplot(aes(x = position, y = avg_coverage, color=cryptic, fill=cryptic)) +
#   geom_smooth(method = "loess", span=0.2, se = F) +
#   geom_ribbon(data = d3utr_average_coverage  %>% group_by(cryptic) %>%
#                 mutate(ymin_smooth = stats::predict(loess(lwr~position,span=0.1)),
#                        ymax_smooth = stats::predict(loess(upr~position,span=0.1))),
#               aes(ymin = ymin_smooth,  ymax = ymax_smooth, fill = cryptic),
#               alpha = 0.31) +
#   xlab("Position") +
#   ylab("Average Coverage") +
#   facet_wrap(~plot_type, ncol = 2, scales = "fixed")+
#   scale_x_continuous(
#     limits = c(0, 1001),
#     breaks = seq(0,1000,100),
#     labels = as.character(seq(-500,500,100))
#   )
# 
# d3utr_average_coverage %>%
#   ggplot(aes(x = position, y = avg_coverage, color=cryptic, fill=cryptic)) +
#   geom_smooth(method = "loess", span=0.2, se = F) +
#   geom_ribbon(data = d3utr_average_coverage  %>% group_by(cryptic) %>%
#                 mutate(ymin_smooth = stats::predict(loess(lwr~position,span=0.1)),
#                        ymax_smooth = stats::predict(loess(upr~position,span=0.1))),
#               aes(ymin = ymin_smooth,  ymax = ymax_smooth, fill = cryptic),
#               alpha = 0.31) +
#   xlab("Position") +
#   ylab("Average Coverage") +
#   facet_wrap(~plot_type, ncol = 2, scales = "fixed")+
#   scale_x_continuous(
#     limits = c(0, 1001),
#     breaks = seq(0,1000,100),
#     labels = as.character(seq(-500,500,100))
#   )
# 
# d3utr_average_coverage %>%
#   ggplot(aes(x = position, y = avg_coverage, color=cryptic, fill=cryptic)) +
#   geom_smooth(method = "loess", span=0.2, se = F) +
#   geom_ribbon(data = d3utr_average_coverage  %>% group_by(across(all_of(c("plot_type", "cryptic")))) %>%
#                 mutate(ymin_smooth = stats::predict(loess(lwr~position,span=0.1)),
#                        ymax_smooth = stats::predict(loess(upr~position,span=0.1))),
#               aes(ymin = ymin_smooth,  ymax = ymax_smooth, fill = cryptic),
#               alpha = 0.31) +
#   xlab("Position") +
#   ylab("Average Coverage") +
#   geom_vline(xintercept = 500, linetype = "dashed", alpha = 0.5) +
#   facet_wrap(~plot_type, ncol = 1, scales = "fixed") +
#   scale_x_continuous(
#     limits = c(0, 1001),
#     breaks = seq(0,1000,100),
#     labels = as.character(seq(-500,500,100))
#   ) +
#   scale_fill_manual(values = c("#000000", "#d95f02")) +
#   scale_color_manual(values = c("#000000", "#d95f02")) +
#   theme_bw(base_size = 14) +
#   theme(legend.position = "top")