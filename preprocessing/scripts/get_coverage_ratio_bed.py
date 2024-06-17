#!/usr/bin/env python3


import pyranges as pr
import pandas as pd
import numpy as np
import sys
import argparse

def main(pas_bed: str,
         upstream_bed: str,
         downstream_bed: str,
         output_bed: str,
         window_size: int = 200,
         exclude_pas: bool = False
         ) -> None:
    '''_summary_

    Parameters
    ----------
    pas_bed : str
        BED file of PAS (all PAS that were extended up/downstream)
    upstream_bed : str
        BED file of PAS-upstream extended regions with mean coverage added to Score field ( produced by add_coverage_to_bed.py)
    downstream_bed : str
        BED file of PAS-downstream extended regions with mean coverage added to Score field ( produced by add_coverage_to_bed.py)
    output_bed : str
        Name of output BED-like containing regions from pas_bed with upstream, downstream and upstream:downstream coverage appended as 7th, 8th and 9th fields
    window_size : int, optional
        Size of extension added to intervals from pas_bed to generate upstream_bed and downstream_bed, by default 200
    exclude_pas : bool, optional
        _description_, by default False
    '''

    # read in BED files
    pas = pr.read_bed(pas_bed, as_df=True)
    up_bed = pr.read_bed(upstream_bed, as_df=True)
    down_bed = pr.read_bed(downstream_bed, as_df=True)

    # strip upstream/downstream suffix from Name fields
    up_bed.loc[:, "Name"] = up_bed["Name"].str.removesuffix("|upstream")
    down_bed.loc[:,"Name"] = down_bed["Name"].str.removesuffix("|downstream")

    # check that both dfs have completely the same regions
    pas2up = set(pas.Name).difference(set(up_bed.Name))
    pas2down = set(pas.Name).difference(set(down_bed.Name))

    assert len(pas2up) == 0, "missing intervals between input PAS bed and upstream bed (according to 'Name' field)"
    assert len(pas2down) == 0, "missing intervals between input PAS bed and downstream bed (according to 'Name' field)"

    # Subset window intervals to the PAS coordinate (upstream = 3'most coordinate, downstream = 5'most coordinate)
    up_bed = pr.PyRanges(up_bed).three_end()
    down_bed = pr.PyRanges(down_bed).five_end()

    # Shift coordinate by 1 to match initial PAS coordinate if excluded from windows
    if exclude_pas:
        # PAS is at end of interval, shift downstream by 1 (strand-aware) 
        up_bed = up_bed.extend({"3": 1}).subsequence(start=-1, strand=True)
        # PAS is at start of interval, shift upstream by 1 (strand-aware)
        down_bed = down_bed.extend({"5": 1}).subsequence(end=1, strand=True)

    up_bed = up_bed.as_df()
    down_bed = down_bed.as_df()

    # combine region BEDs with original BED file, adding Score column (mean coverage for given window)
    pas = pas.merge(up_bed, how="left", on=["Chromosome", "Start", "End", "Strand", "Name"], validate="one_to_one", suffixes=[None, "_up"])
    pas = pas.merge(down_bed, how="left", on=["Chromosome", "Start", "End", "Strand", "Name"], validate="one_to_one", suffixes=[None, "_down"])

    # calculate ratio between coverage of upstream vs downstream region
    pas.loc[:, "up_down_ratio"] = pas["Score_up"] / pas["Score_down"]

    # replace nans if coverage 0 for both regions (leave otherwise)
    pas.loc[:, "up_down_ratio"] = np.where((pas["Score_up"].eq(0)) & (pas["Score_down"].eq(0)), 0.0, pas["up_down_ratio"])

    # output to BED-like file
    pas.to_csv(output_bed, sep="\t", index=False, header=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculate upstream vs downstream coverage ratio from BED files.")
    
    parser.add_argument("-p", "--pas_bed", required=True, help="BED file of original PAS (all PAS that were extended up/downstream)")
    parser.add_argument("-u", "--upstream_bed", required=True, help="BED file of PAS-upstream extended regions with mean coverage added to 'Score' field (produced by add_coverage_to_bed.py)")
    parser.add_argument("-d", "--downstream_bed", required=True, help="BED file of PAS-downstream extended regions with mean coverage added to 'Score' field (produced by add_coverage_to_bed.py)")
    parser.add_argument("-o", "--output_bed", required=True, help="Name of output BED-like containing regions from pas_bed with upstream, downstream and upstream:downstream coverage appended as 7th, 8th and 9th fields")
    parser.add_argument('--window_size', type=int, default=200, help='Size of the windows to extend used to generate upstream and downstream BEDs, default is 200.')
    parser.add_argument('--exclude_pas', action='store_true', help='Whether the pas nucleotide was removed from each window to prevent overlapping region between two windows, default is False')
    
    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    
    main(args.pas_bed, args.upstream_bed, args.downstream_bed, args.output_bed, args.window_size, args.exclude_pas)
