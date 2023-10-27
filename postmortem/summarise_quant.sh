#!/usr/bin/env bash

# Need this conda env active - envs/papa_r.yaml - papa_pipeline_r

script="scripts/tx_to_polya_quant.R"
sample_tbl="data/liu_facs/liu_facs_sample_sheet_dexseq_meta.csv"
salmon_dir="data/liu_facs/2023-08-23_salmon_reference_filtered_cryptic_le_combined/salmon_quant"
tx2le="processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.tx2le.tsv"
le2gene="processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.tx2gene.tsv"
out_prefix="processed/liu_facs/2023-10-27_summarised_pas_quantification_cryptics"

Rscript $script \
-s $sample_tbl \
-d $salmon_dir \
-t $tx2le \
-g $le2gene \
-o $out_prefix &> ${out_prefix}.tx_to_polya_quant.log
