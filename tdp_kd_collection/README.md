# Analysis of PAPA results for in-vitro TDP-43 depletion compendium

Scripts to analyse and visualise cryptic polyA detection across compendium of neuronal TDP-43 KD RNA-seq datasets

## plot_kd_collection.R

Generates scatter plot of median cryptic usage across TDP-43 KD datasets (Fig 1)

Takes combined output of `../preprocessing/scripts/clean_papa_tbls.R` and summary of cryptic events (output by `../preprocessing/scripts/manual_validation_summary.R`).

## plot_cryptic_consistency.R

Exploration and visualisation of consistency of cryptic APA expression across compendium of KD datasets. Generates the heatmap (Supplementary Fig 1) and detection/cryptic status bin plot (Supplementary Figure 2A).

Takes combined output of `../preprocessing/scripts/clean_papa_tbls.R` and summary of cryptic events (output by `../preprocessing/scripts/manual_validation_summary.R`).

## seddighi_de_volcano.R

Generates volcano plot of differential expression results from 'Seddighi i3Neuron' dataset, highlighting cryptic APA genes (Fig 3, Supplementary Fig 4).

Requires DESeq2 output results table (available at `data/seddighi.ipscCortical_neuron.DESEQ2_results.csv`) and summary of cryptic events (output by `../preprocessing/scripts/manual_validation_summary.R`).
