#!/bin/bash

# Check if the output prefix is passed as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <output_prefix>"
    echo "output_prefix - path to/name of prefix for output files. For transcript assignment script, 'tx_assignment.cryptic_<ids/genes>' is added"
    echo "Manually edit input files if wish to change"
    echo "Ensure 'pybioinfo' conda environment is active"
    exit 1
fi

# Assign the first positional argument to output_prefix
OUTPUT_PREFIX=$1

# Define input file variables
GTF_FILE="processed/2024-09-27_cryptic_last_exons.gtf"
REFERENCE_GTF="data/reference_filtered.gtf"
DECOYS_MODE="full"
METADATA_FILE="data/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv"
TX2LE_IDS="data/novel_ref_combined.cryptic_ids.tx2le.tsv"
TX2LE_GENE="data/novel_ref_combined.cryptic_genes.tx2le.tsv"
OUTPUT_TX2LE="${OUTPUT_PREFIX}.tx2le.tsv"

# skip running decoys if already exists
if [ -f "$OUTPUT_TX2LE" ]; then
    echo "Output file '$OUTPUT_TX2LE' already exists. Skipping add_decoys_to_gtf.py."
else
    # Run the first Python script (add_decoys_to_gtf.py)
    echo 'Running decoys script'
    time python scripts/add_decoys_to_gtf.py \
        -i "$GTF_FILE" \
        -r "$REFERENCE_GTF" \
        -d "$DECOYS_MODE" \
        -o "$OUTPUT_PREFIX"
fi

echo 'Checking assignment of cryptic transcripts to le_ids (just cryptic IDs (.cryptic_ids) and all IDs of cryptic genes)'
# Run the second Python script twice with different output prefixes
# First run: cryptic_ids
echo
echo 'Checking transcript to le_id assignment, focusing on cryptic IDs only'
echo
python scripts/test_decoys_tx_assignment.py \
    --output_prefix "${OUTPUT_PREFIX}.tx_assignment.cryptic_ids" \
    --metadata_file "$METADATA_FILE" \
    "$TX2LE_IDS" \
    "$OUTPUT_TX2LE"

# Second run: cryptic_genes
echo 
echo 'Checking transcript to le_id assignment, focusing on all IDs of cryptic genes'
echo
python scripts/test_decoys_tx_assignment.py \
    --output_prefix "${OUTPUT_PREFIX}.tx_assignment.cryptic_genes" \
    --metadata_file "$METADATA_FILE" \
    "$TX2LE_GENE" \
    "$OUTPUT_TX2LE"

