# Snakemake pipeline to run PEKA

Note: this pipeline uses a forked version of PEKA to allow running on arbitrary coordinates and report additional outputs. The underlying algorithms are not changed in any way.

To run on provided data:

- download a copy of genome sequence (FASTA format) & generate a fasta index (e.g. samtools faidx)
- update the 'genome' and 'genome_idx' values in config.yaml to corresponding paths

submit a dry run with
`snakemake -n -p --configfile config/config.yaml`

If everything looks good, run the pipeline with following command. Can take up to 5 min for PEKA to run. conda env can also take some time to install on first run
`snakemake -p --configfile config/config.yaml --use-conda --cores 1`

For own runs, the pipeline requires the following inputs:

- 'sample table' - see 'config/sample_table.csv' for an example, and config.yaml for description of column names/fields
- BED file of single nucleotide coordinates containing two groups of interest. The 'foreground' group are identified by an extract string present in the corresponding Name field. All entries not containing the exact string are considered as 'background' regions.
- Genome fasta and FAI index
- Regions GTF genome segmentation file produced as output of "iCount segment" function. For purposes of arbitrary groups of coordinates, the preprovided test GTF file (grabbed from peka repo, provided here in `data` is sufficient)



### cv_coverage.smk

Uses the same config file, but this is useful if you're just interested in getting % coverage for a set of motifs (aggregated & optionally per-motif) in foreground and background regions

submit a dry run with
`snakemake -n -p --configfile config/config.yaml -s cv_coverage.smk`

If everything looks good, run the pipeline with following command. 
`snakemake -p --configfile config/config.yaml -s cv_coverage.smk --use-conda --cores 1`

In contrast to when cv_coverage is run under the main Snakefile, the outputs of these runs are stored under `cv_coverage` in the main output dir specified in the config file. If per_motif, the aggregate results are stored under `cv_coverage_<foreground/background>`. Otherwise, the per-motif coverages are stored under `<main_output_dir/cv_coverage/<kmer>/cv_coverage_<foreground/background>`