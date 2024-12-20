#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import sys
import argparse
import os

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

    strand2strand = {"plus": "+", "minus": "-"}

    # read in coverage TSVs output by megadepth (1 per strand)
    cov_tsvs = [pd.read_csv(cov_path, sep="\t",names=["Chromosome", "Start", "End", "Score"]) for cov_path in coverage_tsvs]
    
    # infer strand from the filename ('plus'/'minus', output by snakemake rule as <sample_name>.<megadepth>.<strand>.annotation.tsv)
    strands_tsvs = [os.path.basename(cov_path).split(".")[-3] for cov_path in coverage_tsvs]
    assert all(strand in strand2strand.keys() for strand in strands_tsvs), print(f"Invalid inferred strand from coverage TSV file names, should follow <prefix>.<plus/minus>.annotation.tsv")
    
    # convert plus/minus to +/-
    strands_tsvs = [strand2strand[strand] for strand in strands_tsvs]

    # Add strand as column
    cov_tsvs = [df.assign(Strand=pd.Series([strand]*len(df))) for df, strand in zip(cov_tsvs, strands_tsvs)]

    # combine to single df of all regions
    cov = pd.concat(cov_tsvs, ignore_index=True)

    # dups can occur if same PAS for different regions/identifiers
    cov.drop_duplicates(inplace=True)

    # read in combined regions file used as input
    regions = pr.read_bed(region_bed, as_df=True)
    regions.drop(columns="Score", inplace=True)

    # Add coverage column as score for each interval
    # Note: multiple matches for same PAS coordinate are permitted (e.g. if different Name identifiers)
    regions = regions.merge(cov, on=["Chromosome", "Start", "End", "Strand"], how="left")
    
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