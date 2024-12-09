library(tidyverse)


bed_dir <- "data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/rel_quant_beds"
bed_regex <- "_04\\.bed$"
ctrl_regex <- "^NT_"
kd_regex <- "^TDP43_"
outdir <- "processed/dapars_comparison/"

# list full paths to files in a directory matching a regex (and name vector with filename as name)
get_file_paths <- function(path, regex) {
  # List files in the directory that match both prefix and suffix patterns
  files <- list.files(path = path, pattern = regex, full.names = TRUE)
  # Name the vector with the basename of each file for reference
  names(files) <- basename(files)
  return(files)
}

# Wrapper function to read in BED file and add sample_name as column
read_bed_file <- function(file_path) {
  # Read the TSV file, assign columns, and add filename as sample_name
  read_tsv(file_path, col_names = c("chr", "start", "end", "name", "score", "strand"), show_col_types = F) %>%
    mutate(sample_name = basename(file_path))  # Add sample_name column
}


# Get file paths for control and knockdown files
ctrl_files <- get_file_paths(bed_dir, paste0(ctrl_regex, ".*", bed_regex))
kd_files <- get_file_paths(bed_dir, paste0(kd_regex, ".*", bed_regex))

# Read and combine files for each condition
ctrl_data <- map(ctrl_files, read_bed_file) %>%
  bind_rows() %>%
  mutate(condition = "CTRL")

kd_data <- map(kd_files, read_bed_file) %>%
  bind_rows() %>%
  mutate(condition = "TDP-43 KD")

# Combine both conditions into a single tibble for analysis
comb_bed <- bind_rows(ctrl_data, kd_data) %>%
  mutate(sample_name = str_remove_all(sample_name, "_DaPars2_04.bed$")) %>% 
  separate(name, into = c("tx_prefix", "tx_id", "gene_name", "pas"), sep = "_", remove = F)

# Plot the distal PAS usage across conditions for each of the genes (when provided with )
distal_pdui_plot <- comb_bed %>%
  filter(pas == "distal") %>%
  ggplot(aes(x = gene_name, y = score*100, colour = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(seed = 123), size = 0.75) +
  scale_y_continuous(breaks = seq(0,80,10)) +
  scale_colour_manual(values = c("#999999", "#E69F00")) +
  labs(x = "",
       y = "Distal PAS usage (%, DaPars2)",
       colour = "Condition") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")


distal_pdui_plot

if (!dir.exists(outdir)) {dir.create(outdir, recursive = T)}


ggsave(filename = file.path(outdir, "2024-11-11_updated.elk1_six3_tlx1.pdui.seddighi_i3_cortical.png"),
       plot = distal_pdui_plot,
       width = 100, height = 100, units = "mm", dpi = "retina"
       )

ggsave(filename = file.path(outdir, "2024-11-11_updated.elk1_six3_tlx1.pdui.seddighi_i3_cortical.pdf"),
       plot = distal_pdui_plot,
       width = 100, height = 100, units = "mm", dpi = "retina"
)

