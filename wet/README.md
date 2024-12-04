# Visualisation and statistical analysis of wet lab data

## fish_process_counts.R

Cleans provided table of FISH quantification results in format suitable for downstream analysis + visualisation. Calculates mean foci counts across images in each replicate and subcellular foci count ratio.

takes as input:

- Counts table `fish_counts_cleaned.tsv` provided by Puja (available under `data`)

## fish_plot.R

Performs normalisation and statistical test for foci counts and ratios. Produces plots of foci counts for main figure (Fig 3). Takes input from `fish_process_counts.R`

## wb_plot.R

Produces scatter plot of quantification of ELK1 Western Blot intensities in control vs knockdown i3Neurons.

takes as input:

- Precomputed intensities in each replicate normalised to tubulin intensities (available under `data`)