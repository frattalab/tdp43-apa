#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
import sys
import argparse


default_bed_cols = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]

def main(cryptics_csv: str, output_file: str):

    df = pd.read_csv(cryptics_csv)

    pas = df.PAS_ID.str.split(":", expand = True)
    pas.columns = ["Chromosome", "Start", "Strand"]
    pas = pas.astype({"Start": np.int64})
    pas.loc[:, "End"] = pas.Start + 1

    # add in metadata
    pas = pas.merge(df, left_index=True, right_index=True)
    # print(pas.columns)
    # print(pas)

    # Create 'Name' field - combine PAS_ID with event type, gene_name, control usage and KD usage (after rounding)
    pas.loc[:, pas.select_dtypes(include=['float64', 'float32']).columns] = pas.select_dtypes(include=['float64', 'float32']).apply(lambda x: round(x, 5))
    # print(pas)
    pas.loc[:, "Name"] = pas["PAS_ID"].str.cat(pas[["Feature", "gene_name", "PolyA_usage_control", "PolyA_usage_TDP-43_knockdown"]].astype(str), sep="|")

    # add dummy Score column to satisfy BED format
    pas.loc[:, "Score"] = "."
    
    # convert to pyranges and output to BED
    pas = pas[default_bed_cols]
    bed = pr.PyRanges(pas)
    # print(bed)

    bed.to_bed(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert cryptic PAS table (Zeng et al. 2024 preprint, Supplementary Table 5, doi: 10.1101/2024.01.22.575730) to BED file.")
    parser.add_argument("cryptics_csv", help="Path to supplementary table containing cryptic PAS from Zeng et al. 2024 preprint (Supplementary Table 5)")
    parser.add_argument("output_file", help="Name of output BED file")
    
    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    
    main(args.cryptics_csv, args.output_file)
