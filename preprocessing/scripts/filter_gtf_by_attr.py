#!/usr/bin/env python3

# Python script to filter a GTF file using pyranges

import pyranges as pr
import argparse
import logging
import sys


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Filter GTF file by list of attribute values (e.g. gene names)")
    parser.add_argument("input_file", help="Path to the input GTF file")
    parser.add_argument("output_file", help="Path to the output GTF file")
    parser.add_argument("-a", "--attribute_name", required=True, help="Attribute name/key to filter by (e.g. gene_id,gene_name)")
    parser.add_argument("-v", "--values", required=True, nargs='+', help="Space-separated list of allowed values")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    # read in small part of GTF file - check that provided attribute_name is present as column in GTF
    assert args.attribute_name in pr.read_gtf(args.input_file, nrows=100).columns, f"The attribute '{args.attribute_name}' is not present in the GTF file."

    # Load the full GTF file
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Print the parsed command line argument values
    logging.info("Input args: %r", args)

    # Load the full GTF file
    logging.info("Loading the full GTF file...")
    gr = pr.read_gtf(args.input_file)

    # Filter the GTF file
    logging.info("Filtering the GTF file...")
    in_values = set(args.values)
    filtered_gr = gr.subset(lambda df: df[args.attribute_name].isin(in_values))

    # Save the filtered GTF file
    logging.info("Saving the filtered GTF file...")
    filtered_gr.to_gtf(args.output_file)

    logging.info("Filtering complete.")


if __name__ == "__main__":
    main()