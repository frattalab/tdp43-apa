library(optparse)

option_list <- list(make_option(c("-r", "--regions"),
                                type="character",
                                help="Path to BED file of single-nucleotide intervals for which to extend and get coverage from iCLIP peaks file"),
                    make_option(c("-i","--iclip"),
                                type="character",
                                help = "Path to BED file containing iCLIP peaks"),
                    make_option(c("-c","--chromsizes"),
                                type="character",
                                help = "Path to TSV file containing chromosome sizes (chromosome\\tsize)"),
                    make_option(c("-w","--window"),
                                type="integer",
                                default = 500,
                                help = "Window size in both directions around intervals in regions BED file to search for coverage ([default= %default])"),
                    make_option(c("-o", "--output-dir"),
                                dest = "output_dir",
                                default = "iclip_coverage",
                                help = "name of/path to output directory ([default= %default])"),                  
                    make_option(c("-p", "--bedtools-path"),
                                dest = "bedtools_path",
                                default = "bedtools",
                                type = "character",
                                help = "Path to bedtools executable ([default= %default])")
                    )

opt_parser <- OptionParser(option_list = option_list)

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  stop()
}

opt <- parse_args(opt_parser)

# Now running analysis code
# define minimal functions
suppressPackageStartupMessages(library(tidyverse))

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
               "-b", iclip,
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
  
  # split BEDs into two separate groups 
  system(paste(foreground_grep, coverage_file, ">", out_foreground), intern = T)
  system(paste(background_grep, coverage_file, ">", out_background), intern = T)
  

}


# create output directory if required
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = T)
}  





# extend the bed intervals by specified width in either direction

# construct the extended BED output file name
# bed_prefix <- str_remove(basename(opt$regions, "\\.bed$"))
out_bed_ext <- file.path(opt$output_dir, paste("regions.flank_", opt$window, ".bed",sep = ""))

# run bedtools slop on regions
message("Running bedtools slop to extend intervals in regions file...")
extend_bed(opt$regions, opt$chromsizes, flank_interval = 500, outfile = out_bed_ext, bedtools_path = opt$bedtools_path)

# Get single-nucleotide iCLIP coverage over extended intervals

# since coverage files can get pretty big (2*flank_interval lines for each interval in input file), will eventually want to delete
tmp_cov_bed <- file.path(opt$output_dir, paste("regions.flank_", opt$window, ".coverage.txt",sep = ""))

# run bedtools coverage on extended windows
message("Running bedtools coverage on extended windows...")
get_coverage(out_bed_ext,
             iclip = opt$iclip,
             outfile = tmp_cov_bed,
             bedtools_path = opt$bedtools_path
             )


# now split both proximal/distal into cryptic/background
message("splitting into cryptic and background regions...")
split_coverage(coverage_file = tmp_cov_bed,
               outfile_foreground_suffix = ".cr.txt",
               outfile_background_suffix = ".bg.txt",
               outfile_prefix = str_remove(tmp_cov_bed, "\\.txt$")
               )

# delete proximal bed
message("removing temporary coverage TXT file...")
file.remove(tmp_cov_bed)

# finally gzip the output TXT files to save on disk space
# gzip all the output files
message("gzipping the output coverage TXT files...")
list.files(opt$output_dir,
           "\\.((bg)|(cr))\\.txt$",
           full.names = T,
           recursive = T) %>%
  walk(~ system(paste("gzip", .x), intern = T))



