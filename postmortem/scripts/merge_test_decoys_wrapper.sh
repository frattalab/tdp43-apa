#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Usage: ./merge_test_decoys_wrapper.sh <output_prefix>

OUTPUT_PREFIX=$1

# Paths to the input GTF files
# different permutations of cryptic events, non-cryptic genes and non-cryptic events within cryptic-containing genes
# i.e. output of scripts/split_gtf_by_cryptic_and_event_status.py
SPLIT_GTF_PREFIX="processed/decoys/novel_ref_combined.quant"
NON_CRYPTICS_GTF="${SPLIT_GTF_PREFIX}.non_cryptics.all.all.gtf"
CRYPTIC_OTHERS_GTF="${SPLIT_GTF_PREFIX}.cryptics.all.non_cryptic_ids.gtf"
CRYPTIC_3EXT_GTF="${SPLIT_GTF_PREFIX}.cryptics.ext3.ids.gtf"
CRYPTIC_PROX3EXT_GTF="${SPLIT_GTF_PREFIX}.cryptics.proxext3.ids.gtf"
CRYPTIC_COMPLEX_GTF="${SPLIT_GTF_PREFIX}.cryptics.complex.ids.gtf"

# output of scripts/get_decoy_tx.py - contains IPA, ALEs & their decoy transcripts
DECOYS_GTF="processed/decoys/playground/2024-11-12_decoys.combined.gtf"

# Paths to the metadata and output files
METADATA_FILE="data/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv"
TX2LE_IDS="data/novel_ref_combined.cryptic_ids.tx2le.tsv"
TX2LE_GENE="data/novel_ref_combined.cryptic_genes.tx2le.tsv"
OUTPUT_TX2LE="${OUTPUT_PREFIX}.tx2le.tsv"

# Call the merge_decoy_gtfs.py script
echo 'MERGING INPUT GTFS'
echo
# python scripts/merge_decoy_gtfs.py --non_cryptics $NON_CRYPTICS_GTF --decoys $DECOYS_GTF --cryptic_others $CRYPTIC_OTHERS_GTF --output_prefix $OUTPUT_PREFIX
python scripts/merge_decoy_gtfs.py \
    --non_cryptics $NON_CRYPTICS_GTF \
    --decoys $DECOYS_GTF \
    --cryptic_others $CRYPTIC_OTHERS_GTF \
    --cryptic_3exts $CRYPTIC_3EXT_GTF \
    --cryptic_prox3exts $CRYPTIC_PROX3EXT_GTF \
    --cryptic_complex $CRYPTIC_COMPLEX_GTF \
    --output_prefix $OUTPUT_PREFIX

# Call the test_decoys_tx_assignment.py script with TX2LE_IDS
echo
echo 'CHECKING TRANSCRIPT TO LE_ID ASSIGNMENT, FOCUSING ON CRYPTIC IDS ONLY'
echo
python scripts/test_decoys_tx_assignment.py \
    --metadata_file $METADATA_FILE \
    --output_prefix "${OUTPUT_PREFIX}.tx_assignment.cryptic_ids" \
    "$TX2LE_IDS" \
    "$OUTPUT_TX2LE"

# Call the test_decoys_tx_assignment.py script with TX2LE_GENE
echo 
echo 'CHECKING TRANSCRIPT TO LE_ID ASSIGNMENT, FOCUSING ON ALL IDS OF CRYPTIC GENES'
echo
python scripts/test_decoys_tx_assignment.py \
    --metadata_file $METADATA_FILE \
    --output_prefix "${OUTPUT_PREFIX}.tx_assignment.cryptic_genes" \
    "$TX2LE_GENE" \
    "$OUTPUT_TX2LE"