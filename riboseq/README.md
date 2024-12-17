# Ribo-seq and Frac-seq analysis

- [Prerequisites](#prerequisites)
- [Ribo-seq Analysis](#ribo-seq-analysis)
- [Frac-seq Analysis](#frac-seq-analysis)

## Prerequisites

- If it's your first time running R code from this project, run `renv::restore()` in your R console to install the required dependencies.

## Ribo-seq Analysis

### Ribo-seq differential expression analysis

script: `deseq2_riboseq.R`

The script performs a differential expression analysis of the i3Neuron Ribo-seq data between TDP-43 knockdown and control with DESeq2. Requires count matrices of reads assigned to annotated coding regions by featureCounts, as computed using our ['feature_counts' Snakemake pipeline](https://github.com/frattalab/rna_seq_single_steps).

### Cryptic APA set enrichment in Ribo-seq differential expression results (GSEA)

script: `gsea_riboseq.R`

Script to perform Gene Set Enrichment Analysis (GSEA) on Ribo-seq differential expression results, using cryptic APA event types as gene sets. Requires output of `deseq2_riboseq.R`

### Ribo-seq analysis plotting

script: `plot_riboseq.R`

Generate volcano plot and GSEA dot plot from Ribo-seq analyses.

## Frac-seq Analysis

### Compute TPM, count and % PAS usage matrices from Salmon output

script: `fracseq_tx_to_polya_quant.R`

Uses tximport to aggregate transcript TPMs across samples to compute isoform & gene-level TPM matrices. Requires a directory containing Salmon quantifications as computed using the ['single_steps' Salmon pipeline](https://github.com/frattalab/rna_seq_single_steps), a 'dummy' sample table (just mapping sample names to experimental fraction and cell line) and transcript to PAS/gene assignment tables (i.e. 'tx2gene', 'tx2le', 'le2gene') for the PAPA quantification GTF.

Essentially a copy of the 'tx_to_polya_quant.R' script from the main PAPA repository/pipeline modified to remove the command line interface.

### Visualise cryptic 3'Exts in Frac-seq data

script: `plot_fracseq.R`

Transformation and plotting code to produce Frac-seq plot (Supplementary Figures)
