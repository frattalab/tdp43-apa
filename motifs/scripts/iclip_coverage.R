library(tidyverse)

#' run bedtools slop on a BED file to extend interval by user-specified distance
extend_bed <- function(bed, chromsizes, flank_interval, outfile, bedtools_path = "/home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools") {
  
  if (!dirname(outfile) == ".") {
    system(paste("mkdir -p", dirname(outfile)), intern = T)
  }  
  
  # run bedtools slop
  system(paste(bedtools_path, "slop -b", flank_interval, "-i", bed, "-g", chromsizes, ">", outfile), intern = T)
}




#' run bedtools coverage between extended BED file & a bed file of iCLIP peaks
#' coverage args - c('-d') is non-negotiable (reports cov for every position in interval)
get_coverage <- function(extended_bed,
                         iclip, 
                         outfile = "coverage_output.txt",
                         coverage_args = c("-d"),
                         bedtools_path = "/home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools") {

  if (!dirname(outfile) == ".") {
    system(paste("mkdir -p", dirname(outfile)), intern = T)
  }  
  
  # get a space separated string of additional bedtools coverage arguments
  cov_args <- paste(coverage_args, collapse = " ")
  
  # coverage_file <- paste(out_prefix, ".txt",sep = "")
  # # system(paste("mkdir -p", bed_dir), intern = T)
  # # system(paste(bedtools_path, "slop -b", flank_interval, "-i", bed_file, "-g", chromsizes, ">", bed_extend_file), intern = T)
  system(paste(bedtools_path, "coverage",
               cov_args,
               "-a", extended_bed,
               "-b", pooled_iclip,
               ">", outfile),
         intern = T)
  
}

#' Split output of bedtools coverage into separate files of foreground/background coverage
split_coverage <- function(coverage_file,
                           outfile_foreground_suffix,
                           outfile_background_suffix,
                           outfile_prefix = "coverage_output",
                           foreground_grep = "grep cryptic",
                           background_grep = "grep background") {
  
  # construct output file paths
  out_foreground <- paste(outfile_prefix, outfile_foreground_suffix, sep = "")
  out_background <- paste(outfile_prefix, outfile_background_suffix, sep = "")

  # # bg_file <- file.path(bed_dir,  "coverage_output.prox.bg.txt")
  # # cr_file <- file.path(bed_dir, "coverage_output.prox.cr.txt")
  # bg_file <- paste(bed_dir,)
  
  # dist_bg_file <- file.path(bed_dir, "coverage_output.dist.bg.txt")
  # dist_cr_file <- file.path(bed_dir, "coverage_output.dist.cr.txt")
  # 
  system(paste(foreground_grep, coverage_file, ">", out_foreground), intern = T)
  system(paste(background_grep, coverage_file, ">", out_background), intern = T)
  
  # system(paste("grep background", coverage_file, "| grep prox >", prox_bg_file), intern = T)
  # system(paste("grep background", coverage_file, "| grep dist >", dist_bg_file), intern = T)
  # system(paste("grep crypt", coverage_file, "| grep prox >", prox_cr_file), intern = T)
  # system(paste("grep crypt", coverage_file, "| grep dist >", dist_cr_file), intern = T)
  
}

iclip1 <- "data/iCLIP/tardbp-shsy5y-1-20210701-mh_mapped_to_genome_single_peaks.sort.bed"
iclip2 <- "data/iCLIP/tardbp-shsy5y-2-20210701-mh_mapped_to_genome_single_peaks.sort.bed"
pooled_iclip <- "data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed"
chromsizes <- "data/GRCh38.primary_assembly.genome.chromsizes.txt"
bedtools_path <- "/home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools"

# run bedtools slop for PAS BED files
d3utr_pas_bed <- "data/iCLIP/2023-07-04_papa_cryptic_d3utr.pas.bed"
spliced_pas_bed <- "data/iCLIP/2023-07-04_papa_cryptic_spliced.pas.bed"
bleedthrough_pas_bed <- "data/iCLIP/2023-07-04_papa_cryptic_bleedthrough.pas.bed"

extend_bed(d3utr_pas_bed, chromsizes, flank_interval = 500, outfile = "data/iCLIP/2023-07-04_papa_cryptic_d3utr.pas.flank_500.bed")
extend_bed(spliced_pas_bed, chromsizes, flank_interval = 500, outfile = "data/iCLIP/2023-07-04_papa_cryptic_spliced.pas.flank_500.bed")
extend_bed(bleedthrough_pas_bed, chromsizes, flank_interval = 500, outfile = "data/iCLIP/2023-07-04_papa_cryptic_bleedthrough.pas.flank_500.bed")

# Run coverage on PAS BED files
cov_outdir <- "processed/iclip_maps/"
if (!dir.exists(cov_outdir)) {dir.create(cov_outdir, recursive = T)}

# since coverage files can get pretty big (2*flank_interval lines for each interval in input file), will eventually want to delete
tmp_d3utr_cov_bed <- file.path(cov_outdir, "d3utr", "2023-07-04_papa_cryptic_spliced.pas.flank_500.coverage.txt")
tmp_spliced_cov_bed <- file.path(cov_outdir, "spliced", "2023-07-04_papa_cryptic_spliced.pas.flank_500.coverage.txt")
tmp_bleedthrough_cov_bed <- file.path(cov_outdir, "bleedthrough", "2023-07-04_papa_cryptic_bleedthrough.pas.flank_500.coverage.txt")

# 
get_coverage("data/iCLIP/2023-07-04_papa_cryptic_d3utr.pas.flank_500.bed",
             iclip = pooled_iclip,
             outfile = tmp_d3utr_cov_bed,
             )

# first split into proximal/distal
split_coverage(coverage_file = tmp_d3utr_cov_bed,
               outfile_foreground_suffix = ".prox.txt",
               outfile_background_suffix = ".dist.txt",
               outfile_prefix = str_remove(tmp_d3utr_cov_bed, ".txt$"),
               foreground_grep = "grep prox",
               background_grep = "grep dist")

# delete original/temp cov bed
file.remove(tmp_d3utr_cov_bed)

# now split both proximal/distal into cryptic/background
split_coverage(coverage_file = paste(str_remove(tmp_d3utr_cov_bed, ".txt$"), ".prox.txt", sep = ""),
               outfile_foreground_suffix = ".cr.txt",
               outfile_background_suffix = ".bg.txt",
               outfile_prefix = paste(str_remove(tmp_d3utr_cov_bed, ".txt$"), ".prox", sep = "")
               )

# delete proximal bed
file.remove(paste(str_remove(tmp_d3utr_cov_bed, ".txt$"), ".prox.txt", sep = ""))

split_coverage(coverage_file = paste(str_remove(tmp_d3utr_cov_bed, ".txt$"), ".dist.txt", sep = ""),
               outfile_foreground_suffix = ".cr.txt",
               outfile_background_suffix = ".bg.txt",
               outfile_prefix = paste(str_remove(tmp_d3utr_cov_bed, ".txt$"), ".dist", sep = "")
               )

# now delete distal all bed
file.remove(paste(str_remove(tmp_d3utr_cov_bed, ".txt$"), ".dist.txt", sep = ""))

# repeat for spliced
get_coverage("data/iCLIP/2023-07-04_papa_cryptic_spliced.pas.flank_500.bed",
             iclip = pooled_iclip,
             outfile = tmp_spliced_cov_bed,
)

# extract cryptic/bg
split_coverage(coverage_file = tmp_spliced_cov_bed,
               outfile_foreground_suffix = ".cr.txt",
               outfile_background_suffix = ".bg.txt",
               outfile_prefix = str_remove(tmp_spliced_cov_bed, ".txt$"))

file.remove(tmp_spliced_cov_bed)

# repeat for bleedthrough
get_coverage("data/iCLIP/2023-07-04_papa_cryptic_bleedthrough.pas.flank_500.bed",
             iclip = pooled_iclip,
             outfile = tmp_bleedthrough_cov_bed,
)

# extract cryptic/bg
split_coverage(coverage_file = tmp_bleedthrough_cov_bed,
               outfile_foreground_suffix = ".cr.txt",
               outfile_background_suffix = ".bg.txt",
               outfile_prefix = str_remove(tmp_bleedthrough_cov_bed, ".txt$"))

file.remove(tmp_bleedthrough_cov_bed)

# gzip all the output files
list.files(cov_outdir,
           "\\.((bg)|(cr))\\.txt$",
           full.names = T,
           recursive = T) %>%
  walk(~ system(paste("gzip", .x), intern = T))


## repeat for 5'ends of bleedthroughs and spliced events

# input beds
spliced_start_bed <- "data/iCLIP/2023-07-04_papa_cryptic_spliced.le_start.bed"
bleedthrough_start_bed <- "data/iCLIP/2023-07-04_papa_cryptic_bleedthrough.le_start.bed"

# extend single nuc coords in either direction
extend_bed(spliced_start_bed, chromsizes, flank_interval = 500, outfile = "data/iCLIP/2023-07-04_papa_cryptic_spliced.le_start.flank_500.bed")
extend_bed(bleedthrough_start_bed, chromsizes, flank_interval = 500, outfile = "data/iCLIP/2023-07-04_papa_cryptic_bleedthrough.le_start.flank_500.bed")

tmp_spliced_cov_bed <- file.path(cov_outdir, "spliced", "2023-07-04_papa_cryptic_spliced.le_start.flank_500.coverage.txt")
tmp_bleedthrough_cov_bed <- file.path(cov_outdir, "bleedthrough", "2023-07-04_papa_cryptic_bleedthrough.le_start.flank_500.coverage.txt")

# Get coverage for spliced events
get_coverage("data/iCLIP/2023-07-04_papa_cryptic_spliced.le_start.flank_500.bed",
             iclip = pooled_iclip,
             outfile = tmp_spliced_cov_bed,
)

# extract cryptic/bg
split_coverage(coverage_file = tmp_spliced_cov_bed,
               outfile_foreground_suffix = ".cr.txt",
               outfile_background_suffix = ".bg.txt",
               outfile_prefix = str_remove(tmp_spliced_cov_bed, ".txt$"))

file.remove(tmp_spliced_cov_bed)

# repeat for bleedthrough
get_coverage("data/iCLIP/2023-07-04_papa_cryptic_bleedthrough.le_start.flank_500.bed",
             iclip = pooled_iclip,
             outfile = tmp_bleedthrough_cov_bed,
)

# extract cryptic/bg
split_coverage(coverage_file = tmp_bleedthrough_cov_bed,
               outfile_foreground_suffix = ".cr.txt",
               outfile_background_suffix = ".bg.txt",
               outfile_prefix = str_remove(tmp_bleedthrough_cov_bed, ".txt$"))

file.remove(tmp_bleedthrough_cov_bed)

# gzip all the output files
list.files(cov_outdir,
           "\\.((bg)|(cr))\\.txt$",
           full.names = T,
           recursive = T) %>%
  walk(~ system(paste("gzip", .x), intern = T))
