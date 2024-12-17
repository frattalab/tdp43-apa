# Preprocessing PAPA results tables and in-silico PAPA evaluation

Scripts and analysis workflows to perform additional cleaning of PAPA outputs across in vitro datasets, evaluate cryptic PAS support by poly(A)-tail containing reads and run DaPars2 with cryptic 3'Ext coordinates.

NOTE: this directory contains many incomplete/experimental analyses that are retained only for archival purposes. Please use with caution. If a script is not described in README then assume it is not relevant to the main manuscript.

- [Prerequesites](#prerequesites)
- [Generating cleaned cryptic summary tables](#generating-cleaned-cryptic-summary-tables)
- [poly(A)-tail containing read (PATR) evaluation](#polya-tail-containing-read-patr-evaluation)
- [DaPars2 comparison](#dapars2-comparison)
- [Incomplete analysis and miscellaneous commands not relevant to main manuscript](#incomplete-analysis-and-miscellaneous-commands-not-relevant-to-main-manuscript)

## Prerequesites

- If its your first time running R code from this project, run `renv::restore()` in your R console to install the required dependencies
- Python/command line scripts can be executed using the general `pybioinfo` conda environment. Please see base repository README for installation instructions

## Generating cleaned cryptic summary tables

### clean_papa_tbls.R

Script to re-process, add additional output and combine PAPA DEXSeq results tables

1. Combined df containing results for all experiments (alogn with per-experiment)
2. De-duplicating gene name values (due to a stupid design decision on my part)
3. Re-annotating bleedthrough events if called as novel extension of annotated bleedthrough event

The script requires as input:

- A directory containing subdirectories of individual PAPA runs. The subdirectories should be structured as output by PAPA
  - Minimal files are the DEXSeq results tables output by PAPA (`dexseq_apa.results.processed.tsv`)
  - Script assumes you have different comparisons under different directories. The directory structure should also match that of a PAPA run e.g. <main_output_dir>/differential_apa/dexseq_apa.results.processed.tsv. The script grabs the experiment_name from 'main_output_dir', so make sure it is named appropriately.

### get_cryptics.R

Takes combined output of clean_papa_tbls.R and filters for cryptic events in each experiment, producing a simplified table with summary statistics on cryptic status.

Cryptic critera are hardcoded:

- padj < 0.05
- mean_PPAU_base < 0.1 (mean usage in control cells < 10 %)
- delta_PPAU_treatment_control > 0.1 (change in usage (treatment - control) > 10 %)

### manual_validation_summary.R

Filters out IPA/bleedthrough events that fail manual validation (i.e. intron retention artefacts) and outputs filtered cryptic summary table. Works with output of `get_cryptics.R` along with mv files provided under `data`

## poly(A)-tail containing read (PATR) evaluation

An external [Snakemake pipeline](https://github.com/SamBryce-Smith/bulk_polyatail_reads) is used to extract PATRs from BAM files and generate BEDs of PAS clusters. Please see linked repository for usage/running instructions.

### Count number of unique 3'end coordinates for each last exon ID (le_id)

`scripts/get_num_pas.py`. Currently, script applies only to any le_id evaluated for differential usage across experiments. Outputs number of unique 3'end coordinates per le_id with and without a overlap slack of 12nt (i.e. intervals within 12nt are considered as one). Takes as input:

- PAPA quantification GTF
- Combined TSV of DEXSeq results tables across experiments - produced by `scripts/clean_papa_tbls.R`
- Summary TSV of cryptic events - produced by `scripts/manual_validation_summary.R`

Output file: `processed/curation/cryptic_annot_comparison/2024-09-03_le_id_pas_counts.tsv`

```bash
conda activate pybioinfo # if applicable
python scripts/get_num_pas.py
```

### Calculate KD sample summarised PAS expression (TPM) for each dataset

`scripts/get_papa_summarised_tpm.R` calculates mean and median TPM values across KD samples for each dataset and le_id. Requires as input:

- Matrices of per le_id TPMs as computed by PAPA e.g. `differential_apa/summarised_pas_quantification.tpm.tsv` for each dataset
- Sample tables for each dataset as used by PAPA

### Get (samples of) annotated PAS matched to cryptic PAS for distribution of average expression and the number of unique PAS

Script: `scripts/get_matched_annot_pas.R`

Uses matchRanges to generate covariate matched samples of annotated events (n = 1000). Requires outputs of `scripts/get_papa_summarised_tpm.R` and `scripts/get_num_pas.py` to define covariates. Median TPM matrix is used to represent expression for each le_id (modelled using log2(TPM + 1)). 

Outputs a directory of text files (1 per sample) containing a list of covariate-matched annotated PAS/le_ids.

### Evaluate the efficacy of covariate matching

`scripts/analyse_matched_annot_pas.R` evaluates the annotated PAS samples for the covariate balance and their uniqueness. Not reported directly in manuscript, although essential to show that the balancing is excellent for all 1k samples and pairwise comparisons indicate most samples have approximately 50 % unique events (i.e. not just sampling the same events each time). Used as basis to go ahead with retaining all 1k samples to define an empirical distribution.

Warning: some of the steps are a little slow and are computed sequentially/without parallel processing.

### Get PATR PAS overlap with cryptic and (covariate-matched) annotated le_ids across range of distance thresholds

script: `scripts/get_patr_matches.py`

Feeds on output of `scripts/get_num_pas.py` and `get_matched_annot_pas.R`. Takes directory of different samples of covariate-balanced annotated le_ids (each stored in a separate text file), and one-by-one filters the PAPA GTF for these events and performs an overlap with PAS defined by polyA-tail reads, reporting as a binary 1/0 for overlap at a series of overlap thresholds.

Warning - not particularly well optimised. Takes approx 10 mins to run for 1k annotated ID samples and 7 distance thresholds

```bash
conda activate pybioinfo # if applicable
python scripts/get_patr_matches.py --papa-gtf data/novel_ref_combined.last_exons.gtf --patr-bed data/bulk_polya_reads/tdp_ko_collection/pas_clusters/condition__TDP43KD/two_class_simple/polya_clusters.bed --outdir processed/curation/cryptic_annot_comparison/ --pas-counts processed/curation/cryptic_annot_comparison/2024-09-03_le_id_pas_counts.tsv --annot-ids-dir processed/curation/cryptic_annot_comparison/ids/ --distance-thresholds 0,10,25,50,100,200,500 --nproc 8
```

### Visualise differences in overlap between annotated and cryptic PAS at each distance threshold

`scripts/compare_cryptic_annot_patr_matches.R` visualises the % overlap reported for cryptic relative to the 1k annotated PAS samples. Also computes an empirical p-value at each distance threshold assessing the null that cryptic and annotated PAS originate from the same distribution.

## DaPars2 comparison

### Download RefSeq GTF

Downloaded RefSeq transcript annotations from the UCSC genome browser FTP website (NCBI RefSeq uses different chromosome naming convention). Used the 110 release archive (believe it's the last using hg38?)

```bash
# assuming already in tdp43-apa/preprocessing
cd data/dapars_comparison
curl https://hgdownload.soe.ucsc.edu/goldenPath/archive/hg38/ncbiRefSeq/110/hg38.110.ncbiRefSeq.gtf.gz -o hg38.110.ncbiRefSeq.gtf.gz
```

### Filtering GTF to target genes

Filter GTF for entries corresponding to specified gene names (ELK1, SIX3 and TLX1)

```bash
conda activate pybioinfo # if applicable
python scripts/filter_gtf_by_attr.py -v ELK1 SIX3 TLX1 -a gene_name data/dapars_comparison/hg38.110.ncbiRefSeq.gtf.gz data/dapars_comparison/hg38.110.ncbiRefSeq.filtered.elk1_six3_tlx1.gtf
```

### Constructing extended transcript models for ELK1, SIX3 and TLX1

`scripts/get_dapars_ext_tx.py` overlaps input 3'Ext coordinates with annotated transcripts. Generates an extra extended event where the 3'end coordinates of matched annotated last exons are updated to the cryptic 3'Ext 3'end coordinate. The updated last exons are then returned to the original full transcript model. Outputs 'proximal', original annotated PAS (and distal PAS) alongside the updated models transcript.

```bash
conda activate pybioinfo # if applicable
python scripts/get_dapars_ext_tx.py -g data/dapars_comparison/hg38.110.ncbiRefSeq.filtered.elk1_six3_tlx1.gtf -b ../motifs/processed/iclip_regions/2023-12-15_3ext.last_exons.cryptic.bed -o processed/dapars_comparison/hg38.110.ncbiRefSeq.filtered.elk1_six3_tlx1
```

### Running DaPars2

Analysed Seddighi i3 cortical samples using the [APAeval snakemake pipeline](https://github.com/iRNA-COSI/APAeval/tree/main/method_workflows/DaPars2/Dapars2_snakemake) (commit hash: d7831b6). Ran separately using original annotated transcript models and the extended transcript models. The APAeval pipeline parses the DaPars2 output into BED files. I used the '01.bed' output files for each sample for futher analysis (single nucleotide PAS predictions, score field filled with placeholder value ('.')).

### Combining predictions across replicates (condition and annotation-wise)

Wrote a simple script to combine the BED files. Duplicate predictions are collapsed, arbitrarily selecting a representative interval (will only differ by 'Name' field, which contains the corresponding transcript ID for predictions). Origin of intervals/predictions is not tracked/reported. Commands used reported below (executed from `preprocessing` subdirectory with general conda environment active)

```bash
$ python scripts/combine_dapars_beds.py -i data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43*01.bed -o data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.exts.kd.bed 
Provided input BED files:  ['data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065403_S23_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065405_S46_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065407_S25_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065409_S29_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065411_S54_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065413_S19_DaPars2_01.bed']
Successfully merged 6 BED files to: data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.exts.kd.bed
$ python scripts/combine_dapars_beds.py -i data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43*01.bed -o data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.orig.kd.bed 
Provided input BED files:  ['data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065403_S23_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065405_S46_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065407_S25_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065409_S29_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065411_S54_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065413_S19_DaPars2_01.bed']
Successfully merged 6 BED files to: data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.orig.kd.bed
$ python scripts/combine_dapars_beds.py -i data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT*01.bed -o data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.orig.ctrl.bed 
Provided input BED files:  ['data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074709_S20_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074711_S17_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074713_S59_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074715_S37_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074717_S21_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074719_S22_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074721_S50_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074723_S36_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074725_S27_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074727_S44_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074729_S12_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074731_S13_DaPars2_01.bed']
Successfully merged 12 BED files to: data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.orig.ctrl.bed
```

### Visualise DaPars2-inferred cryptic PAS usage

`plot_dapars_pdui.R` takes the combined BED files and plots the DaPars2-inferred cryptic PAS (distal) usage split by condition for each gene.

## Incomplete analysis and miscellaneous commands not relevant to main manuscript

Here lies incomplete analysis I started, but did not finish during the revision process. Keeping here for archival purposes in case I need to return. Includes qualitative comparisons with Zeng et al. i3Neuron Quant-seq and attempts to refine reported cryptic PAS coordinates using PATRs.


### Putative updated 3'ends using polyA junction reads as input

```bash
python scripts/create_putative_pas.py -g data/novel_ref_combined.last_exons.gtf -b data/bulk_polya_reads/tdp_ko_collection/pas_clusters/condition__TDP43KD/two_class_simple/polya_clusters.bed --max_distance 10000 --use-bed-name -o processed/curation/2024-06-20_last_exons.max_distance_10000.rep
```

### Get BED file of cryptic PAS from Zeng et al Supplementary Table 5

```bash
mkdir -p processed/zeng_2024/ && python scripts/zeng_cryptics_to_bed.py data/zeng_2024/supplementary_table_s5.csv processed/zeng_2024/supplementary_s5.cryptic_pas.bed
```

### Combined BED file of original + updated putative 3'ends

Used as input to coverage.smk

```bash
$ basename $PWD
curation
```

```python
import pyranges as pr
>>> orig = pr.read_bed("2024-06-20_last_exons.max_distance_10000.rep.original.bed")
>>> upd = pr.read_bed("2024-06-20_last_exons.max_distance_10000.rep.updated.bed")
>>> comb = pr.concat([orig, upd])
>>> comb_3p = comb.three_end()
# Doing this just in case, drops ~ 25k rows
>>> comb_3p = comb_3p.apply(lambda df: df.drop_duplicates())
>>> comb_3p.sort().to_bed("2024-06-20_pas.max_distance_10000.original_updated_combined.bed")
```