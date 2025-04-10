#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import argparse
import sys
from typing import Tuple

def parse_coordinates_from_name(name: str) -> Tuple[str, int, int, str]:
    """Extract coordinates from the name field of format chr:strand:start:end"""
    try:
        chrom, strand, start, end = name.split(':')
        return chrom, int(start), int(end), strand
    except ValueError as e:
        print(f"Error parsing coordinates from name: {name}")
        print(f"Expected format: chr:strand:start:end")
        sys.exit(1)

def process_bed_file(input_file: str, output_file: str) -> None:
    """
    Read BED file, update coordinates based on name field, and write output
    
    Args:
        input_file: Path to input BED file
        output_file: Path to output BED file
    """
    try:
        # Read BED file using PyRanges
        bed = pr.read_bed(input_file)
        
        # Convert to pandas DataFrame for easier manipulation
        df = bed.as_df()
        
        # Extract coordinates from Name column and update Start/End
        for idx, row in df.iterrows():
            _, start, end, _ = parse_coordinates_from_name(row['Name'])
            df.at[idx, 'Start'] = start
            df.at[idx, 'End'] = end
        
        # Convert back to PyRanges and write output
        result = pr.PyRanges(df)
        result.to_bed(output_file)
        
    except Exception as e:
        print(f"Error processing BED file: {str(e)}")
        sys.exit(1)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Update PATR clusters BED file to representative positions from the name field'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input BED file path'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output BED file path'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process the file
    process_bed_file(args.input, args.output)

if __name__ == "__main__":
    main()
