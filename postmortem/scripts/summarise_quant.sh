#!/usr/bin/env bash

# Need this conda env active - envs/papa_r.yaml - papa_pipeline_r

script="scripts/tx_to_polya_quant.R"
sample_tbl="data/liu_facs/liu_facs_papa_sample_sheet.tdppos_first.csv"
salmon_dir="data/liu_facs/2024-11-20_decoys/salmon_quant/"
tx2le="processed/decoys/2024-11-20_decoys_novel_ref_combined.quant.tx2le.tsv"
le2gene="processed/decoys/2024-11-20_decoys_novel_ref_combined.quant.tx2gene.tsv"
out_prefix="processed/liu_facs/2024-11-20_summarised_pas_quantification_cryptics"

Rscript --vanilla $script \
-s $sample_tbl \
-d $salmon_dir \
-t $tx2le \
-g $le2gene \
-o $out_prefix &> ${out_prefix}.tx_to_polya_quant.log
