#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import sys
import argparse


def main(coverage_tsvs: list,
         region_bed: str,
         output_bed: str):
    '''_summary_

    Parameters
    ----------
    coverage_tsvs : list
        Paths to coverage TSV files output by megadepth ('<prefix>.annotation.tsv') for 'plus' and 'minus' strand intervals
    region_bed : str
        Original BED of intervals passed to megadepth (containing both plus and minus strand intervals)
    output_bed : str
        Name of output BED file containing region_bed with coverage value added to 'score' field of BED file (5th field)
    '''

    # read in coverage TSVs output by megadepth (1 per strand)
    cov_tsvs = [pd.read_csv(cov_path, sep="\t",names=["Chromosome", "Start", "End", "Score"]) for cov_path in coverage_tsvs]

    # combine to single df of all regions
    cov = pd.concat(cov_tsvs, ignore_index=True)

    # read in combined regions file used as input
    regions = pr.read_bed(region_bed, as_df=True)
    regions.drop(columns="Score", inplace=True)

    # Add coverage column as score for each interval
    regions = regions.merge(cov, on=["Chromosome", "Start", "End"], how="left", validate="one_to_one")
    
    # convert to BED format and output to file
    regions = pr.PyRanges(regions, int64=True).sort()

    # output to file
    regions.to_bed(output_bed)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Add summarised coverage values output per interval + strand by megadepth to the score field of the original interval BED file')

    parser.add_argument(
        'coverage_tsvs', 
        nargs='+', 
        help="Paths to coverage TSV files output by megadepth ('<prefix>.annotation.tsv') for 'plus' and 'minus' strand intervals"
    )
    parser.add_argument('-r', '--regions', 
        help="Original BED of intervals passed to megadepth (containing both plus and minus strand intervals)"
    )
    parser.add_argument('-o', '--output_bed', 
        help="Name of output BED file containing region_bed with coverage value added to 'score' field of BED file (5th field)"
    )

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.coverage_tsvs, args.regions, args.output_bed)