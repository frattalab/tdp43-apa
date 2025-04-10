#!/bin/bash

usage() {
   echo "Usage: $0 [-w window_size] [-o out_bed_dir] <bed_files_directory> <reference.bed> <output.tsv>"
   echo "Options:"
   echo "  -w  Window size in base pairs (default: 0)"
   echo "  -o  Output directory for overlapping BED files (optional)"
   exit 1
}

# set default window-size (to be updated)
WINDOW_SIZE=0
OUT_BED_DIR=""

while getopts "w:o:h" opt; do
   case $opt in
       w) WINDOW_SIZE=$OPTARG ;;
       o) OUT_BED_DIR=$OPTARG ;;
       h) usage ;;
       ?) usage ;;
   esac
done

shift $((OPTIND-1))

if [ "$#" -ne 3 ]; then
   usage
fi

BED_DIR="$1"
REFERENCE="$2"
OUTPUT="$3"

if [ ! -d "$BED_DIR" ]; then
   echo "Error: Directory $BED_DIR does not exist"
   exit 1
fi

if [ ! -f "$REFERENCE" ]; then
   echo "Error: Reference file $REFERENCE does not exist"
   exit 1
fi

if [ -n "$OUT_BED_DIR" ] && [ ! -d "$OUT_BED_DIR" ]; then
   mkdir -p "$OUT_BED_DIR"
fi

# Create header for output file
echo -e "file_name\ttotal_unique\toverlapping\tpercent_overlapping\toverlapping_score_sum" > "$OUTPUT"

# Get total number of bed files for progress
total_files=$(ls -1 "$BED_DIR"/*.bed | wc -l)
current_file=0

# Process each BED file in the directory
for bed_file in "$BED_DIR"/*.bed; do
    if [ -f "$bed_file" ]; then
        # Update progress
        ((current_file++))
        echo "Processing file $current_file of $total_files: $(basename "$bed_file")"

        # Get filename without path
        filename=$(basename "$bed_file")

        # Count total unique intervals in file A (after removing duplicates)
        total_unique=$(sort -k1,1 -k2,2n -k3,3n -k6,6 "$bed_file" | uniq | wc -l)

        # Count overlapping intervals (after removing duplicates from A)
        # Store overlapping regions in temporary file/output directory if saving
        if [ -n "$OUT_BED_DIR" ]; then
            overlapping_regions="$OUT_BED_DIR/${filename%.bed}.overlapping.bed"
        else
            overlapping_regions=$(mktemp)
        fi

        sort -k1,1 -k2,2n -k3,3n -k6,6 "$bed_file" | uniq | \
        bedtools window -a stdin -b "$REFERENCE" -sm -u -w "$WINDOW_SIZE" > "$overlapping_regions"

        overlapping=$(wc -l < "$overlapping_regions")
        score_sum=$(awk '{sum += $5} END {print sum}' "$overlapping_regions")

        percent=$(awk "BEGIN {printf \"%.2f\", ($overlapping/$total_unique)*100}")

        echo -e "${filename}\t${total_unique}\t${overlapping}\t${percent}\t${score_sum}" >> "$OUTPUT"

        if [ -z "$OUT_BED_DIR" ]; then
            rm "$overlapping_regions"
        fi

    fi
done

echo "Analysis complete. Results written to $OUTPUT"