# Analysis of ELK1 TF activity

Notes:

- Contains a mix of finalised and WIP/unfinished scripts. Anything not referenced below is not required to reproduce analysis in scripts, and is retained purely for archival purposes. There is no guarantee of correct installation, running or validity of outputs.

## Prerequisites

- If it's your first time running R code from this project, run `renv::restore()` in your R console to install the required dependencies.

## Export Dorothea/Progeny targets to disk

scripts: `scripts/dl_collectri.R` and  `scripts/dl_progeny.R`

Prerequisite to GSEA analysis. Loads in the target tables from the respective packages and exports them to TSVs

## Perform GSEA analysis of HeLa TDP-43 KO differential expression

scripts: `script/hela_ko_tf_activity_gsea.R`

Using ELK1 HeLa ChIP-seq target genes as gene lists, evalutate whether they are overrepresented among differentially expressed genes using GSEA. Outputs analysis results.

Requires ELK1 and ELK4 ChIP-seq target lists (available under `data`) and the DESeq2 differential expression results table for the Ferguson HeLa dataset.

## Visualise GSEA results

script: `plot_hela_ko_tf_activity_gsea.R`

Generates GSEA running sum enrichment plot with p-value annotations as in Fig 3.

## Output combined lists of target genes from ELK1 & ELK4 HeLa ChIP-seq

script: `scripts/write_lists_hela_ko.R`

Makes a TSV with target genes and their intersection/'uniqueness' in each column. Generates table included in Supplementary Material