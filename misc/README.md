# Miscellaneous scripts

## combine_supp_data.R

Script to generate supplemental data file. Requires as input:

- table of dataset characteristics (provided with repo)
- cryptic ALE SJ quantification and selectivity classification (produced by `postmortem/scripts/clean_ale_junction_counts.R`)
- Riboseq differential expression results for cryptic APAs (produced by `postmortem/scripts/plot_riboseq.R`)
- BED files of representative last exon coordinates for each APA category (produced by `motifs/notebooks/define_iclip_regions_<ale|d3utr>.ipynb` - both scripts used)
- Summary dataframe of cryptic event annotation and expression across in-vitro datasets (produced by `preprocessing/scripts/manual_validation_summary.R`)
- List of target genes for ELK1 & ELK4 in HeLa ChIP-seq data (produced by `tf_activity/scripts/write_lists_hela_ko.R`)