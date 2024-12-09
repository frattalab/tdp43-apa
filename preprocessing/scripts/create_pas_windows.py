#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import argparse
import sys


def main(bed_path: str,
         output_prefix: str,
         window_size: int = 200,
         return_merged: bool = False,
         exclude_pas: bool = False,
         upstream_key: str = "upstream",
         downstream_key: str = "downstream",
         name_sep: str = "|"):
    '''Create windows extended up and downstream of 3'end coordinates and output to BED file

    Parameters
    ----------
    bed_path : str
        _description_
    output_prefix : str
        _description_
    window_size : int, optional
        _description_, by default 200
    return_merged : bool, optional
        _description_, by default False
    exclude_pas : bool, optional
        _description_, by default False
    upstream_key : str, optional
        _description_, by default "upstream"
    downstream_key : str, optional
        _description_, by default "downstream"
    name_sep : str, optional
        _description_, by default "|"
    '''

    bed = pr.read_bed(bed_path)

    # get 3'ends
    bed_3p = bed.three_end()
    
    upstream = bed_3p.extend({"5": window_size})
    downstream = bed_3p.extend({"3": window_size})

    # return combined region if applicable (generate now so retain original Name field)
    if return_merged:
        merged = pr.concat([upstream, downstream]).merge(strand=True, by="Name").sort()

    # shift windows to exclude the PAS coordinate from up and downstream regions (so no shared nucleotides/coverage)
    if exclude_pas:
        # PAS is at end of interval, remove the final coord (strand-aware)
        upstream = upstream.subsequence(end=-1, strand=True)
        # PAS is at start of interval, remove the first coord (strand-aware)
        downstream = downstream.subsequence(start=1, strand=True)

    # create updated name field
    upstream.Name = upstream.Name + name_sep + upstream_key
    downstream.Name = downstream.Name + name_sep + downstream_key

    # write to file
    upstream.to_bed(".".join([output_prefix, upstream_key, "bed"]))
    downstream.to_bed(".".join([output_prefix, downstream_key, "bed"]))
    if return_merged:
        merged.to_bed(".".join([output_prefix, "merged", "bed"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create windows extended up and downstream of 3\'end coordinates and output to BED file.')
    parser.add_argument('bed_path', type=str, help='Path to the input BED file.')
    parser.add_argument('output_prefix', type=str, help='Prefix for the output BED files.')
    parser.add_argument('--window_size', type=int, default=200, help='Size of the windows to extend, default is 200.')
    parser.add_argument('--return_merged', action='store_true', help='Return merged regions, default is False.')
    parser.add_argument('--exclude_pas', action='store_true', help='Remove the pas nucleotide from each window to prevent overlapping region between two windows, default is False')
    parser.add_argument('--upstream_key', type=str, default='upstream', help='Key to append to Name field for upstream regions (and output file suffix), default is "upstream".')
    parser.add_argument('--downstream_key', type=str, default='downstream', help='Key to append to Name field for downstream regions (and output file suffix), default is "downstream".')
    parser.add_argument('--name_sep', type=str, default='|', help='Separator for the new name field, default is "|".')


    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(bed_path=args.bed_path,
         output_prefix=args.output_prefix,
         window_size=args.window_size,
         return_merged=args.return_merged,
         exclude_pas=args.exclude_pas,
         upstream_key=args.upstream_key,
         downstream_key=args.downstream_key,
         name_sep=args.name_sep)