#!/usr/bin/env python3

import pyranges as pr
import sys

def main(in_bed,
         blacklist,
         out_bed):
    

    bed = pr.read_bed(in_bed)

    with open(blacklist, "r") as infile:
        bl_ids = [line.rstrip("\n") for line in infile]

    
    bed = bed.subset(lambda df: ~df.Name.isin(bl_ids))

    bed.to_bed(out_bed)


if __name__ == '__main__':
    
    
    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
        print("python3 filter_bed_name.py INPUT_BED BLACKLIST_NAMES_TXT OUTPUT_BED")
        sys.exit()

    main(sys.argv[1], sys.argv[2], sys.argv[3])
