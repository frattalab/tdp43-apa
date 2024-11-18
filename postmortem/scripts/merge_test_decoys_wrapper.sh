#!/bin/bash

# Usage: ./merge_test_decoys_wrapper.sh <output_prefix>

OUTPUT_PREFIX=$1

# Paths to the input GTF files
NON_CRYPTICS_GTF="processed/decoys/novel_ref_combined.quant.non_cryptics.all.all.gtf"
DECOYS_GTF="processed/decoys/playground/2024-11-12_decoys.combined.gtf"
CRYPTIC_OTHERS_GTF="processed/decoys/novel_ref_combined.quant.cryptics.all.non_cryptic_ids.gtf"

# Paths to the metadata and output files
METADATA_FILE="data/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv"
TX2LE_IDS="data/novel_ref_combined.cryptic_ids.tx2le.tsv"
TX2LE_GENE="data/novel_ref_combined.cryptic_genes.tx2le.tsv"
OUTPUT_TX2LE="${OUTPUT_PREFIX}.tx2le.tsv"

# Call the merge_decoy_gtfs.py script
echo 'MERGING INPUT GTFS'
echo
python scripts/merge_decoy_gtfs.py --non_cryptics $NON_CRYPTICS_GTF --decoys $DECOYS_GTF --cryptic_others $CRYPTIC_OTHERS_GTF --output_prefix $OUTPUT_PREFIX

# Call the test_decoys_tx_assignment.py script with TX2LE_IDS
echo
echo 'CHECKING TRANSCRIPT TO LE_ID ASSIGNMENT, FOCUSING ON CRYPTIC IDS ONLY'
echo
python scripts/test_decoys_tx_assignment.py \
    --metadata $METADATA_FILE \
    --output_prefix "${OUTPUT_PREFIX}.tx_assignment.cryptic_ids" \
    "$TX2LE_IDS" \
    "$OUTPUT_TX2LE"

# Call the test_decoys_tx_assignment.py script with TX2LE_GENE
echo 
echo 'CHECKING TRANSCRIPT TO LE_ID ASSIGNMENT, FOCUSING ON ALL IDS OF CRYPTIC GENES'
echo
python scripts/test_decoys_tx_assignment.py \
    --metadata $METADATA_FILE \
    --output_prefix "${OUTPUT_PREFIX}.tx_assignment.cryptic_genes" \
    "$TX2LE_GENE" \
    "$OUTPUT_TX2LE"