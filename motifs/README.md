# Scripts to perform motif and iCLIP analysis

Scripts to perform de-novo motif enrichment., motif coveage and iCLIP analysis/visualisation around genomic landmarks

## define_iclip_regions_ale.ipynb and define_iclip_regions.ipynb

Scripts to generate BED files of representative cryptic and background intervals for each event. For details on the collapsing criteria, please see the methods in the manuscript.

The scripts require as input:

- GTF of input novel last exons used to generate combined reference of last exons (full last exon sequence, not just subtracted region) - produced by PAPA
- DEXSeq results df produced by PAPA - use the cleaned TSV produced by `../preprocessing/scripts/clean_papa_tbls.R`
- reference GTF used to identify novel last exons, quantify vs ref (typically the `reference_filtered.gtf` file produced by PAPA)
- GTF of last exon quantification regions used to generate Salmon index (typically named `novel_ref_combined.quant.last_exons.gtf`)
- transcript to le_id mapping table produced by PAPA (typically named `novel_ref_combined.tx2le.tsv`)
- summary df of cryptic events following manual validation of bleedthrough events and TXT file of failed ids - produced by `../preprocessing/scripts/manual_validation_summary.R`

## Creating combined BED file of iCLIP peaks

```bash
conda activate pybioinfo
cd data/iCLIP
cat tardbp-shsy5y-1-20210701-mh_mapped_to_genome_single_peaks.bed tardbp-shsy5y-2-20210701-mh_mapped_to_genome_single_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 5 -o sum > tardbp-shsy5y.concat.sort.chr.bed
```

## spliced_iclip_coverage.sh, bleedthrough_iclip_coverage.sh and d3utr_iclip_coverage.sh

Wrapper scripts to run bedtools coverage on input BED files to get iCLIP coverage for positions relative to genomic landmarks. Requires BED files produced by jupyter notebooks and the combined BED file of iCLIP peaks as input. Calls `scripts/cl_iclip_coverage.R` under the hood.

## iclip_plot.R

Script that imports coverage files produced by wrapper shell scripts, performs necessary preprocessing for visualisation of iCLIP binding profiles (Fig 1). Plotting functions are defined in `scripts/fncs_plot_iclip.R` and imported.

## peka_snakemake - Snakemake pipeline to run PEKA (de-novo kmer enrichment) and cv_coverage (pre-defined motif coverage)

See README.md within directory for run/usage instructions. Works generically, but to replicate workflow in manuscript it requires BED files produced by jupyter notebooks

## papa_process_peka.R

Cleans and merges a set of PEKA output tables present in a directory. Takes output of peka_snakemake pipeline as input.

## papa_gsea_peka.R

Runs a gene set enrichment analysis on motifs ranked by PEKA score using pre-defined groups of motifs (here the 'Halleger et al, 2021 TDP-43 enriched 6mers). Produces GSEA dot plot (Supplementary Figure).

Takes as input 'simplified', merged PEKA output table produced by `papa_process_peka.R`.

## papa_process_cv_coverage.R

Cleans and merges a set of per-motif cv_coverage output tables present in a directory (as produced by cv_coverage.smk in the peka_snakemake subdirectory). Motifs to extract are pre-defined in script ('Halleger et al, 2021 TDP-43 enriched 6mers') and kept per motif (i.e. no group-wise summarisation is performed).

## papa_plot_cv_coverage.R

Produces motif coverage maps around genomic landmarks for specified motif groups (Fig 1). Takes output of `scripts/papa_process_cv_coverage.R` as input. Plotting functions are defined in `scripts/fncs_plot_iclip.R` and `scripts/fncs_plot_iclip.R` and imported.
