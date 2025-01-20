#!/usr/bin/env python3

import argparse
import pandas as pd
import pyranges as pr
from functools import reduce

import sys

def get_sequences(bed_file, fasta_file):
    """
    Extract sequences for intervals in BED file from FASTA file using PyRanges
    """
    try:
        # Read BED file into PyRanges
        gr = pr.read_bed(bed_file)
        
        # Check if all intervals have same length
        lengths = gr.lengths()
        if lengths.nunique() != 1:
            raise ValueError("All intervals must have the same length")
            
        # Extract sequences using PyRanges
        seq = pr.get_sequence(gr, fasta_file)
        gr.seq = seq
        
        return gr
        
    except Exception as e:
        print(f"Error processing input files: {str(e)}", file=sys.stderr)
        sys.exit(1)

def _df_per_bp_nucleotide_freq(df, seq_col):
    '''
    Calculate nucleotide frequencies for each position in a set of sequences
    '''
    # Expand sequence in seq_col into single column per nucleotide
    seq_df = pd.DataFrame(df[seq_col]
                         .apply(lambda x: list(x))
                         .values
                         .tolist())
           
    # Calculate frequency of each nucleotide at every position
    freq_df = seq_df.apply(pd.Series.value_counts).fillna(0)
    
    # Convert index to column
    freq_df = freq_df.reset_index().rename(columns={'index': 'nucleotide'})
              
    return freq_df

def per_bp_nucleotide_content(gr, seq_col="seq", return_counts=False):
    '''
    Calculate per position nucleotide content for set of constant length intervals
    '''
    assert seq_col in gr.columns
    assert gr.lengths().nunique() == 1
              
    # Calculate frequencies across all intervals
    freq_dfs = gr.apply(lambda df: _df_per_bp_nucleotide_freq(df, seq_col),
                       as_pyranges=False)
    
    # Sum counts across all chromosomes/strands
    freq_df = reduce(lambda a,b: a.add(b),
                    [df.set_index('nucleotide') for df in freq_dfs.values()])
    
    if return_counts:
        return freq_df.astype(int)
    
    # Calculate fractions
    freq_df = freq_df.divide(freq_df.sum(axis="index"), axis="columns")
              
    return freq_df

def adjust_position_labels(df, align_mode, interval_length):
    """
    Adjust position labels based on alignment mode
    """
    if align_mode == 'start':
        return df  # Keep original 0-based positions
        
    elif align_mode == 'end':
        # Rename columns to count backwards from 0
        new_cols = range(-(interval_length-1), 1)
        df.columns = new_cols
        
    elif align_mode == 'center':
        if interval_length % 2 == 0:
            raise ValueError("Center alignment requires odd-length intervals")
        
        # Calculate positions relative to center (0)
        mid = interval_length // 2
        new_cols = range(-mid, mid + 1)
        df.columns = new_cols
        
    return df

def main():
    parser = argparse.ArgumentParser(description='Generate per-position nucleotide frequency matrix from fixed-length intervals in a BED file')
    parser.add_argument('bed_file', help='Input BED file of fixed-length intervals')
    parser.add_argument('fasta_file', help='Input FASTA file (preferably indexed with pyfaidx, but will be generated if not provided)')
    parser.add_argument('output_file', help='Output TSV file for per-position frequency matrix')
    parser.add_argument('--align', choices=['start', 'end', 'center'], default='start',
                       help='How to name output columns with (zero-based) positions. If start, 1st column begins at 0. If end, last position column is 0. If center, then the middle position is set to 0 (default: start)')
    
    args = parser.parse_args()
    
    # Load intervals, check uniformity of interval length and extract nucleotide sequences
    gr = get_sequences(args.bed_file, args.fasta_file)
    
    # Calculate frequency matrix
    freq_matrix = per_bp_nucleotide_content(gr, return_counts=True)
    
    # Adjust position labels based on alignment mode
    interval_length = gr.lengths().iloc[0]
    freq_matrix = adjust_position_labels(freq_matrix, args.align, interval_length)
    
    # Save to TSV
    freq_matrix.to_csv(args.output_file, sep='\t', header=True, index=True)

if __name__ == '__main__':
    main()