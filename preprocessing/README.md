# Preprocessing PAPA results tables for downstream analysis

Scripts to perform additional cleaning of PAPA outputs across in vitro datasets

## clean_papa_tbls.R

Script to re-process, add additional output and combine PAPA DEXSeq results tables

1. Combined df containing results for all experiments (alogn with per-experiment)
2. De-duplicating gene name values (due to a stupid design decision on my part)
3. Re-annotating bleedthrough events if called as novel extension of annotated bleedthrough event

The script requires as input:

- A directory containing subdirectories of individual PAPA runs. The subdirectories should be structured as output by PAPA
  - Minimal files are the DEXSeq results tables output by PAPA (`dexseq_apa.results.processed.tsv`)
  - Script assumes you have different comparisons under different directories. The directory structure should also match that of a PAPA run e.g. <main_output_dir>/differential_apa/dexseq_apa.results.processed.tsv. The script grabs the experiment_name from 'main_output_dir', so make sure it is named appropriately.

## get_cryptics.R

Takes combined output of clean_papa_tbls.R and filters for cryptic events in each experiment, producing a simplified table with summary statistics on cryptic status.

Cryptic critera are hardcoded:

- padj < 0.05
- mean_PPAU_base < 0.1 (mean usage in control cells < 10 %)
- delta_PPAU_treatment_control > 0.1 (change in usage (treatment - control) > 10 %)

## manual_validation_summary.R

Filters out IPA/bleedthrough events that fail manual validation (i.e. intron retention artefacts) and outputs filtered cryptic summary table. Works with output of `get_cryptics.R` along with mv files provided under `data`

# Miscellaneous commands

## Putative updated 3'ends using polyA junction reads as input

```bash
python scripts/create_putative_pas.py -g data/novel_ref_combined.last_exons.gtf -b data/bulk_polya_reads/tdp_ko_collection/pas_clusters/condition__TDP43KD/two_class_simple/polya_clusters.bed --max_distance 10000 --use-bed-name -o processed/curation/2024-06-20_last_exons.max_distance_10000.rep
```

## Get BED file of cryptic PAS from Zeng et al Supplementary Table 5

```bash
mkdir -p processed/zeng_2024/ && python scripts/zeng_cryptics_to_bed.py data/zeng_2024/supplementary_table_s5.csv processed/zeng_2024/supplementary_s5.cryptic_pas.bed
```

## Combined BED file of original + updated putative 3'ends

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

## Count number of unique 3'end coordinates for each last exon ID (le_id)

Currently, script applies only to any le_id evaluated for differential usage across experiments. Outputs number of unique 3'end coordinates per le_id with and without a overlap slack of 12nt.

Output file: `processed/curation/cryptic_annot_comparison/2024-09-03_le_id_pas_counts.tsv`

```bash
python scripts/get_num_pas.py
```

## Get PATR PAS overlap with cryptic and (covariate-matched) annotated le_ids across range of distance thresholds

Feeds on output of `scripts/get_num_pas.py` and `get_matched_annot_pas.R`. Takes directory of different samples of covariate-balanced annotated le_ids (each stored in a separate text file), and one-by-one performs an overlap with PAS defined by polyA-tail reads at a series of overlap thresholds.

Warning - not particularly well optimised. Takes approx 10 mins to run for 1k annotated ID samples and 7 distance thresholds

```bash
python scripts/get_patr_matches.py --papa-gtf data/novel_ref_combined.last_exons.gtf --patr-bed data/bulk_polya_reads/tdp_ko_collection/pas_clusters/condition__TDP43KD/two_class_simple/polya_clusters.bed --outdir processed/curation/cryptic_annot_comparison/ --pas-counts processed/curation/cryptic_annot_comparison/2024-09-03_le_id_pas_counts.tsv --annot-ids-dir processed/curation/cryptic_annot_comparison/ids/ --distance-thresholds 0,10,25,50,100,200,500 --nproc 8
```
