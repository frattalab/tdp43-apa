#!/usr/bin/env python3


import pyranges as pr
import pandas as pd
import numpy as np
import sys
import argparse

def main(pas_bed: str,
         upstream_bed: str,
         downstream_bed: str,
         output_bed: str) -> None:
    '''
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

    # combine region BEDs with original BED file, adding Score column (mean coverage for given window)
    pas = pas.merge(up_bed[["Name", "Score"]], how="left", on="Name", validate="one_to_one", suffixes=[None, "_up"])
    pas = pas.merge(down_bed[["Name", "Score"]], how="left", on="Name", validate="one_to_one", suffixes=[None, "_down"])

    # calculate ratio between coverage of upstream vs downstream region
    pas.loc[:, "up_down_ratio"] = pas["Score_up"] / pas["Score_down"]

    # replace nans if coverage 0 for both regions (leave otherwise)
    pas.loc[:, "up_down_ratio"] = np.where((pas["Score_up"].eq(0)) & (pas["Score_down"].eq(0)), 0.0, pas["up_down_ratio"])

    # output to BED-like file
    pas.to_csv(output_bed, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
