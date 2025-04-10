#!/usr/bin/env python3

import pandas as pd
import sys

hlp="""
usage: python polyadb_to_bed.py DB_TXT OUTPUT_BED

Convert a PolyADB text file to BED6 format
"""


if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
    print(hlp)
    sys.exit()


db = pd.read_csv(sys.argv[1], sep="\t")
out_bed = sys.argv[2]

# first construct BED coordinates. Assume that position follows UCSC 1-based formate, so subtract 1 to convert to representative site
db["Start"] = db["Position"] - 1

# construct BED coords
bed = db.rename(columns={"Position": "End", "PAS_ID": "Name", "Mean RPM": "Score"})
print(bed)

# select minimal cols for BED6 file
col_order = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
bed[col_order].to_csv(out_bed, sep="\t", index=False, header=False)