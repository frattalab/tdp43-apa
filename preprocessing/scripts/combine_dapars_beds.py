#!/usr/bin/env python3

import argparse
import pyranges as pr
from pathlib import Path
import sys


def merge_bed_files(input_files, output_file):
    """
    Merge multiple BED files, collapse duplicates, and sort the results.
    """
    # Check if input files exist
    for file_path in input_files:
        if not Path(file_path).exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")

    # Read and merge all BED files
    ranges_list = [pr.read_bed(f) for f in input_files]

    # Concatenate all ranges
    merged = pr.concat(ranges_list)

    # Sort and drop duplicates
    merged = merged.sort()
    merged = merged.drop_duplicate_positions(strand=True)

    # Write to output file
    merged.to_bed(output_file)
    print(f"Successfully merged {len(input_files)} BED files to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Merge arbitrary number of BED files into a single sorted, duplicate-excluded BED file"
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Space-separated list of input BED files",
    )
    parser.add_argument("-o", "--output", required=True, help="Output BED file name")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    print("Provided input BED files: ", args.input)
    merge_bed_files(args.input, args.output)


if __name__ == "__main__":
    main()
