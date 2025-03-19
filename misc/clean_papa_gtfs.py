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
from typing import List, Optional

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
        

def combine_and_process_gtfs(input_files: List[str], output_file: str, 
                            event_type_table: Optional[str] = None,
                            event_type_col: str = "event_type") -> None:
    """
    Process multiple GTF files and combine them into a single output.
    
    Args:
        input_files: List of paths to input GTF files
        output_file: Path to the output GTF file
        event_type_table: Optional path to a TSV file containing event type information
        event_type_col: Column name for event type in the event_type_table (default: "event_type")
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
    
    # Add event type information if provided
    if event_type_table:
        print(f"Adding 'cleaned' event type information from {event_type_table}...")
        # Read the event type table
        event_df = pd.read_csv(event_type_table, sep='\t')
        
        # Subset to required columns
        event_df = event_df[['le_id', event_type_col]]
        
        # Rename event_type_col to 'simple_event_type'
        event_df = event_df.rename(columns={event_type_col: 'simple_event_type'})
        
        # Merge event type information with the combined GTF
        combined_gr = combined_gr.apply(lambda df: pd.merge(
            df, event_df, on='le_id', how='left'
        ))
    
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

    parser.add_argument(
        "-e", "--event_type_table",
        required=False,
        help="Path to a TSV file containing event type information"
    )
    parser.add_argument(
        "-c", "--event_type_col",
        default="event_type",
        help="Column name for event type in the event_type_table (default: 'event_type')"
    )
    
    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()
    
    args = parser.parse_args()
    
    combine_and_process_gtfs(args.input, args.output, args.event_type_table, args.event_type_col)
    

if __name__ == "__main__":
    main()