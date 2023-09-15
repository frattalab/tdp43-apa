
`med_scatter` - scatter plot of any le_id with padj < 0.05 across datasets, cryptics in top left quadrant labeleld in orange if satisfy cryptic criteria using the medians of deltas & base usage. A select few IDs are labelled with gene name (`c("ELK1", "ARHGAP32", "STMN2", "CNPY3", "TLX1", "ANKRD27")`).

`plot_med_df` - dataframe used to generate med_scatter plot. Saved to TSV (with fewer columns) at `processed/2023-09-15_cryptics_scatter_standard_plot_tbl.tsv`

`med_scatter_lab_1dataset` - scatter plot of any le_id with padj < 0.05 across datasets, cryptics in top left quafrant labeleld in orange if satisfy cryptic criteria using the medians of delta changes. IDs that satisfy cryptic criteria in at least 1 dataset, but when using median base usage & delta usage they do not pass the ctryptic cutoff are labelled in purple.

`TODO - output table for avove plot`

`cryp_any_not_med_df` - dataframe storing IDs cryptic in at least 1 dataset but failing the cryptic expression criteria when use median base & delta usages across significant datasets. Saved to TSV at `processed/2023-09-15_cryptics_1dataset_not_median_base_delta_tbl.tsv`
