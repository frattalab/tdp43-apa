#!/usr/bin/env python3

import pyranges as pr
import sys

def main(in_bed,
         blacklist,
         filter_str,
         out_bed):
    

    bed = pr.read_bed(in_bed)

    with open(blacklist, "r") as infile:
        bl_ids = [line.rstrip("\n") for line in infile]

    
    bed = bed.subset(lambda df: ~df.Name.isin(bl_ids))

    if filter_str != "None":
        bed = bed.subset(lambda df: df.Name.str.contains(filter_str, regex=False))

    bed.to_bed(out_bed)


if __name__ == '__main__':
    
    descrpn = """usage: python3 filter_bed_name.py INPUT_BED BLACKLIST_NAMES_TXT NAME_FILTER OUTPUT_BED

    INPUT_BED -
    BLACKLIST_NAMES_TXT - TXT file of IDs (1 per line) to remove from BED file (matched to 'Name' field in BED file)
    NAME_FILTER - str or 'None' filter the Name field for values that contain provided string. If don't wish to do this pass 'None'
    OUTPUT_BED - 
    
    """
    
    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
        print(descrpn)
        sys.exit()

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
