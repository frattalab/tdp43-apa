#!/usr/bin/env python3

import pyranges as pr
import argparse

def extend_intervals(bed_file, extend_length, direction):
    """
    Extend intervals in a BED file by a specified length in the given direction.
    
    Args:
        bed_file (str): Path to input BED file
        extend_length (int): Length to extend intervals by
        direction (str): Direction to extend ('upstream', 'downstream', or 'both')
    
    Returns:
        PyRanges object with extended intervals
    """
    # Read the BED file
    gr = pr.read_bed(bed_file)
    
    # Define extension based on direction
    if direction == 'upstream':
        ext = {"5": extend_length}
    elif direction == 'downstream':
        ext = {"3": extend_length}
    elif direction == 'both':
        ext = extend_length  # When passing an int, PyRanges extends both ends
    
    # Extend the intervals
    gr = gr.extend(ext)
    
    return gr


def main():
    parser = argparse.ArgumentParser(description='Extend intervals in a BED file by a specified distance and in a specified direction.')

    parser.add_argument('input_bed',
                        help='Input BED file')

    parser.add_argument('output_bed',
                        help='Output BED file')

    parser.add_argument('-l', '--length',
                        type=int,
                        required=True,
                        help='Length to extend intervals by')

    parser.add_argument('-d', '--direction',
                        choices=['upstream', 'downstream', 'both'],
                        default='both',
                        help='Direction to extend intervals (default: both)')

    args = parser.parse_args()

    # Extend the intervals
    extended_gr = extend_intervals(args.input_bed, args.length, args.direction)

    # Write to output file
    extended_gr.to_bed(args.output_bed)

if __name__ == '__main__':
    main()
