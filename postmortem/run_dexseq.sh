#!/usr/bin/env bash

# Need this conda env active - envs/papa_r.yaml - papa_pipeline_r

sample_tbl="data/liu_facs/liu_facs_sample_sheet_dexseq_meta.csv"
formulas="data/liu_facs/liu_facs_dexseq_formulas.txt"
counts="processed/liu_facs/2023-10-27_summarised_pas_quantification_cryptics.counts.tsv"
le2gene="processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.le2gene.tsv"
out_prefix="processed/liu_facs/2023-10-27_liu_facs_covariate_TDPneg_TDPpos"
base_key="TDPpos"
contrast_key="TDPneg"
contrast_name="TDPneg_vs_TDPpos"


Rscript scripts/run_dexseq.R \
-i $counts \
-g $le2gene \
-s $sample_tbl \
-f $formulas \
--base-key $base_key \
--contrast-key $contrast_key \
--constrast-name $contrast_name \
-o $out_prefix &> ${out_prefix}.run_dexseq.log

echo 'running simpler formula with just patient covariate'

sample_tbl="data/liu_facs/liu_facs_sample_sheet_dexseq_meta.csv"
formulas="data/liu_facs/liu_facs_dexseq_formulas_patient_covariate.txt"
counts="processed/liu_facs/2023-10-27_summarised_pas_quantification_cryptics.counts.tsv"
le2gene="processed/2023-06-22_cryptics_plus_decoys.decoys_full_fix_tx2le.le2gene.tsv"
out_prefix="processed/liu_facs/2023-11-21_liu_facs_covariate_patient_only_TDPneg_TDPpos"
base_key="TDPpos"
contrast_key="TDPneg"
contrast_name="TDPneg_vs_TDPpos"


Rscript scripts/run_dexseq.R \
-i $counts \
-g $le2gene \
-s $sample_tbl \
-f $formulas \
--base-key $base_key \
--contrast-key $contrast_key \
--constrast-name $contrast_name \
-o $out_prefix &> ${out_prefix}.run_dexseq.log
