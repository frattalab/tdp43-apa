library(tidyverse)

plot_coverage_df <- function(df, ci_se_mult = 1.96, event_col = "plot_type", group_col = "plot_cryptic") {
  
  group_cols <- c(event_col, group_col)
  
  # generate confidence interval values
  df %>%
    mutate(plot_ymin = avg_coverage - (ci_se_mult*se),
           plot_ymax = avg_coverage + (ci_se_mult*se)) %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(ymin_smooth = stats::predict(loess(plot_ymin~position, span=0.1)),
           ymax_smooth = stats::predict(loess(plot_ymax~position, span=0.1))) %>%
    ungroup()
  
}

plot_coverage <- function(df, ci_se_mult = 1.96, event_col = "plot_type", group_col = "plot_cryptic", facet_ncol = 2, fill_colours = c("#000000", "#d95f02"), line_colours = c("#000000", "#d95f02"), fill_lab = "", colour_lab = "", title_lab = "") {
  
  # generate confidence interval values
  plot_df <- plot_coverage_df(df, ci_se_mult, event_col, group_col)
  
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
    scale_y_continuous(limits = c(NA, 0.1),
                       breaks = seq(0, 0.1, 0.02)) +
    scale_fill_manual(values = fill_colours) +
    scale_color_manual(values = line_colours) +
    theme_bw(base_size = 20) +
    theme(legend.position = "top",
          axis.text = element_text(size = rel(1.25)),
          axis.title = element_text(size = rel(1.25)),
          strip.text = element_text(size = rel(1.25))
    ) +
    labs(fill = fill_lab, colour = colour_lab, title = title_lab)
  
  
}
