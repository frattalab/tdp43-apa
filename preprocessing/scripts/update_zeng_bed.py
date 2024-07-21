#!/usr/bin/env python3

import pyranges as pr
import sys

if len(sys.argv) == 1 or any(i in sys.argv for i in ["-h", "--help"]):
    print("Usage: python update_zeng_bed.py INPUT_BED OUTPUT_BED")
    sys.exit()

bed = pr.read_bed(sys.argv[1])

# convert to bed compliant coords by adding 1 to end
# (assuming coordinates lifted from LAPA, which should report in zero-based format (assuming so, as uses Python/PyRanges under the hood...))
bed = bed.assign("End", lambda df: df.End + 1)

bed.sort().to_bed(sys.argv[2])