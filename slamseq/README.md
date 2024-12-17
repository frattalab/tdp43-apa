# SLAM-seq analysis

Scripts to run GrandR and visualise stability of cryptic 3'UTR extension-containing genes.

- [Prerequisites](#prerequisities)
- [gene-level RNA half life estimation](#gene-level-rna-half-life-estimation)
- [isoform-specific RNA half life estimation](#isoform-specific-rna-half-life-estimation)

## Prerequisities

- If it's your first time running R code from this project, run `renv::restore()` in your R console to install the required dependencies.
- Generating regions for isoform-specific SLAM-seq analysis will require the general `pybioinfo` conda environment

## gene-level RNA half life estimation

### Estimating and plotting RNA half-lives between conditions

Script: `estimate_half_life.R`

All-in-one script to run GrandR on Grand-SLAM derived counts and plot estimated half lives for selected genes (Fig 3).

Requires TSV of Grand-SLAM inferred counts, available at `data/grandslam_exp1_grandR.tsv.gz`. You will need to decompress prior to running the script e.g.

```bash
gzip -d -k data/grandslam_exp1_grandR.tsv.gz
```

## isoform-specific RNA half life estimation

### Construct shared and isoform-specific regions for 3'Ext (and others) quantification

scripts: `cryptic_papa_to_regions.ipynb` and `background_papa_to_regions.ipynb`

Constructs featureCounts-compatible GTF consisting of isoform-specific last exon intervals for assigning conversion counts to isoforms. For 3'Exts, the regions are split into two bins where the proximal is the 'shared' region and the distal the cryptic specific region. Each is output with its own gene and exon. Events from Zeng et al. have also been added to the reference (many of which were present in annotation), but were not evaluated further. Requires as input:

- PAPA quantification GTF of last exons
- GENCODE v40 annotation (for Zeng et al. events not evaluated by my pipeline)

### Isoform-specific quantification

Performed using the external [fastq2EZbakR pipeline](https://github.com/isaacvock/fastq2EZbakR).

### Isoform-specific half-life estimation and visualisation

Script: `isoform_specific_analysis.Rmd`

Computes isoform-specific half-lives and associated confidence intervals. ELK1 half-life plot used in manuscript is also generated here. Feeds on output of fastq2EZbakR pipeline.
