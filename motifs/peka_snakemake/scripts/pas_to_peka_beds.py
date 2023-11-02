#!/usr/bin/env python3


import pandas as pd
import sys
import argparse

"""
Alla's instructions:
- BED6 files should be split by proximal and distal
- Create 1 BED file containing all sites
- Then a BED file containing the 'regulated' sites only

The output BED files should have following format:
- name column is a '.'
- score column is filled with arbitrary integer to mimic x-link count (e.g. 10)
"""


hlp = """usage: python pas_to_peka_beds.py PAS_BED OUTPUT_PREFIX

PAS_BED - BED file of pairs of proximal and distal sites for each last exon. The 'name' field (4th) must contain the string 'regulated' to define regulated sites (all others = background), and 'proximal' or 'distal' to define proximal and distal PAS respectively. the 'score' (5th) field must contain signed values in order to split into 'up and down' in treatment groups

OUTPUT_PREFIX - prefix to add to output BED files:
    - <OUTPUT_PREFIX>.<proximal/distal>.background_regulated.<up/down>.bed - combination of all background sites and all regulated sites (proximal/distal) with a given change in direction (up = score column has positive value, down = score column has negative value)
    - <OUTPUT_PREFIX>.<proximal/distal>.regulated.<up/down>.bed - all regulated sites (proximal/distal) with a given change in direction

requires pandas (satisfied with peka conda environment)
"""

parser = argparse.ArgumentParser(
           description="Split BED file by name into two groups for input to PEKA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # Add defaults to end of help strings
    )

parser.add_argument("-f", "--foreground-key",
                    dest="foreground_key",
                    default="regulated",
                    help="exact string in Name field of BED file denoting foreground regions")

parser.add_argument("-g1",
                    "--group1",
                    default=None,
                    help="Optional key to split input regions (by Name field) into separate groups (e.g. if have prox/distal PAS in same file)")

parser.add_argument("-g2",
                    "--group2",
                    default=None,
                    help="Optional key to split input regions (by Name field) into separate groups (e.g. if have prox/distal PAS in same file)")

parser.add_argument("input_bed",
                    type=str,
                    help="Path to input BED file")


parser.add_argument("output_prefix",
                    type=str,
                    help="Prefix to add to output BED files")

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()



bed_path = args.input_bed
foreground_key = args.foreground_key
grp1 = args.group1
grp2 = args.group2
out_prefix = args.output_prefix


bed = pd.read_csv(bed_path, sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand"])

bed_col_order = bed.columns.tolist()

# print(bed)

# define a mask to separate foreground and background regions
foreg_mask = bed["name"].str.contains(foreground_key)


# ensure matches
assert sum(foreg_mask) != 0, f"foreground key - {foreground_key} must be found in at least one interval"

# check if need to split into two groups
if grp1 is not None:
    if grp2 is None:
        raise Exception("group1 & group2 strings must be jointly specified")
    
    else:
        # mask to identify intervals belonging to two specified groups
        grp1_mask = bed["name"].str.contains(grp1)
        grp2_mask = bed["name"].str.contains(grp2)
        grp = True

elif grp2 is not None:
    raise Exception("group1 & group2 strings must be jointly specified")

else:
    # no need to split BED into two groups
    grp = False


# peka requires score & name columns in specific format
# now redefine the score and name columns to 10 & '.' respectively
bed.drop(columns=["score", "name"], inplace=True)
bed.loc[:, "score"] = 10
bed.loc[:, "name"] = "."

if grp:
    # split foreground & bg into two separate groups
    bed_fg_g1 = bed[foreg_mask & grp1_mask]
    bed_bg_g1 = bed[(~foreg_mask) & (grp1_mask)]
    bed_g1 = bed[grp1_mask]

    print(f"Number of foreground intervals in group1 - {len(bed_fg_g1)}")
    print(f"Number of background intervals in group1 - {len(bed_bg_g1)}")
    print(f"Number of total intervals in group1 - {len(bed_g1)}")
    
    bed_fg_g2 = bed[foreg_mask & grp2_mask]
    bed_bg_g2 = bed[(~foreg_mask) & (grp2_mask)]
    bed_g2 = bed[~foreg_mask & grp2_mask]

    print(f"Number of foreground intervals in group2 - {len(bed_fg_g2)}")
    print(f"Number of background intervals in group2 - {len(bed_bg_g2)}")
    print(f"Number of total intervals in group2 - {len(bed_g2)}")

    # write to file
    bed_fg_g1[bed_col_order].to_csv(".".join([out_prefix, grp1, "foreground.bed"]), index=False, sep="\t", header=False)
    bed_bg_g1[bed_col_order].to_csv(".".join([out_prefix, grp1, "background.bed"]), index=False, sep="\t", header=False)
    bed_g1[bed_col_order].to_csv(".".join([out_prefix, grp1, "all.bed"]), index=False, sep="\t", header=False)

    bed_fg_g2[bed_col_order].to_csv(".".join([out_prefix, grp2, "foreground.bed"]), index=False, sep="\t", header=False)
    bed_bg_g2[bed_col_order].to_csv(".".join([out_prefix, grp2, "background.bed"]), index=False, sep="\t", header=False)
    bed_g2[bed_col_order].to_csv(".".join([out_prefix, grp2, "all.bed"]), index=False, sep="\t", header=False)

else:
    # just simple forground & background
    bed_fg = bed[foreg_mask]
    bed_bg = bed[~foreg_mask]

    # write to file
    bed_fg[bed_col_order].to_csv(out_prefix + ".foreground.bed", index=False, sep="\t", header=False)
    bed_bg[bed_col_order].to_csv(out_prefix + ".background.bed", index=False, sep="\t", header=False)
    bed[bed_col_order].to_csv(out_prefix + ".all.bed", index=False, sep="\t", header=False)

