#!/usr/bin/env python3
"""
Standardise and clean PAPA GTF attribute values for output

This script processes multiple GTF files, standardizes attribute fields,
and combines them into a single sorted GTF file.
"""

import argparse
import pyranges as pr
import pandas as pd
import sys
from typing import List

def process_gtf_file(file_path: str) -> pr.PyRanges:
    """
    Process a single GTF file according to specified transformations.
    
    Args:
        file_path: Path to the GTF file
        
    Returns:
        PyRanges object with transformed attributes
    """

    # Read the GTF file into a PyRanges object
    gr = pr.read_gtf(file_path)
    
    # Select the required columns
    required_cols = ["Source", "Feature", "Score", 
                    "ref_gene_id", "transcript_id", 
                    "ref_gene_name", "le_id", "event_type"]
        
   # Check if all required columns exist in the GTF
    missing_cols = [col for col in required_cols if col not in gr.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in {file_path}: {', '.join(missing_cols)}")
        
    # Subset to required columns
    gr = gr[required_cols]
    
    # Remove duplicates while maintaining order of appearance
    gr.ref_gene_name = gr.ref_gene_name.apply(
        lambda x: ",".join(list(dict.fromkeys((x).split(","))))
    )
    return gr
        

def combine_and_process_gtfs(input_files: List[str], output_file: str) -> None:
    """
    Process multiple GTF files and combine them into a single output.
    
    Args:
        input_files: List of paths to input GTF files
        output_file: Path to the output GTF file
    """
    processed_grs = []
    
    # Process each input file
    for file_path in input_files:
        print(f"Processing {file_path}...")
        gr = process_gtf_file(file_path)
        if gr is not None:
            processed_grs.append(gr)
    
    if not processed_grs:
        print("Error: No valid GTF files were processed.", file=sys.stderr)
        sys.exit(1)
    
    # Combine all processed PyRanges objects
    print("Combining GTF files...")
    combined_gr = pr.concat(processed_grs)
    
    # Rename columns to remove 'ref_' prefix
    combined_gr = combined_gr.apply(lambda df: df.rename(
        columns={"ref_gene_id": "gene_id", "ref_gene_name": "gene_name"}
    ))
    
    # drop duplicate coordinates by le_id
    print("Dropping duplicate intervals by le_id...")
    combined_gr = combined_gr.apply(lambda df: df.drop_duplicates(subset=["le_id", "Start", "End"]))
    
    # Sort the GTF and write to output file
    print(f"Sorting the combined GTF and writing to {output_file}...")
    combined_gr.sort().to_gtf(output_file)
    
    print(f"Successfully processed {len(processed_grs)} GTF files.")


def main():
    parser = argparse.ArgumentParser(
        description="Standardise and clean PAPA GTF attribute values for simpler downstream analysis"
    )
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        nargs="+", 
        help="Paths to input GTF files"
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path to output combined GTF file"
    )
    
    args = parser.parse_args()
    
    combine_and_process_gtfs(args.input, args.output)


if __name__ == "__main__":
    main()