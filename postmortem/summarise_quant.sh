#!/usr/bin/env bash

# Need this conda env active - envs/papa_r.yaml

script="scripts/tx_to_polya_quant.R"
sample_tbl="data/liu_facs_papa_sample_sheet.csv"
salmon_dir="data/liu_facs/2023-06-21_gencode_v40_cryptics_decoys_full/salmon_quant"
tx2le="processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.tx2le.tsv"
le2gene="processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.tx2gene.tsv"
out_prefix="processed/2023-06-21_summarised_pas_quantification_cryptics_fix_tx2le"

Rscript $script \
-s $sample_tbl \
-d $salmon_dir \
-t $tx2le \
-g $le2gene \
-o $out_prefix &> ${out_prefix}.tx_to_polya_quant.log
