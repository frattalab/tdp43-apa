# Code for cryptic polydenylation manuscript

Analysis code for cryptic polyadenylation manuscript, currently on biorxiv:

> **TDP-43 loss induces extensive cryptic polyadenylation in ALS/FTD**
>
>*Sam Bryce-Smith, Anna-Leigh Brown, Puja R. Mehta, Francesca Mattedi, Alla Mikheenko, Simone Barattucci, Matteo Zanovello, Dario Dattilo, Matthew Yome, Sarah E. Hill, Yue A. Qi, Oscar G. Wilkins, Kai Sun, Eugeni Ryadnov, Yixuan Wan, NYGC ALS Consortium, Jose Norberto S. Vargas, Nicol Birsa, Towfique Raj, Jack Humphrey, Matthew Keuss, Michael Ward, Maria Secrier, Pietro Fratta*
>
>bioRxiv 2024.01.22.576625; doi: https://doi.org/10.1101/2024.01.22.576625

## Dependencies/Installation

All code has been tested/run in Linux based environments (local = Ubuntu 22.04.1 LTS via Windows Subsystem for Linux 2, remote = UCL Computer Science SGE cluster). We cannot guarantee correct function/installation in other environments.

For R-based analysis, RStudio is strongly recommended to make proper use of R projects and Renv files

The minimal pre-requesites are:

- git
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)/[mamba](https://mamba.readthedocs.io/en/latest/installation.html) (mamba recommended as much quicker!)
- R 4.3.2
- RStudio (you may be able to use another IDE compatible with R projects, but I have no experience of these)

Once these are satisfied, clone and enter the repo locally using the following commands:

```bash
git clone https://github.com/frattalab/tdp43-apa.git
cd tdp43-apa
```

### Python code/general packages

Assuming you have conda/mamba available on your system, you can install all the non-R based dependencies with the following command:

```bash
<conda/mamba> env create -f py_bioinfo_full.yaml
```

Once installation is complete, you can activate the environment with the following command:

```bash
conda activate pybioinfo
```

### R code

Every individual subdirectory has its own [RStudio project](https://support.posit.co/hc/en-us/articles/200526207-Using-RStudio-Projects) and [Renv](https://rstudio.github.io/renv/articles/renv.html) files.

- Open R project of interest (i.e. subdirectory here in repo)

  - [See guide here](https://support.posit.co/hc/en-us/articles/200526207-Using-RStudio-Projects), can open via RStudio's 'Open Project' command or by opening the '.Rproject' file in your system's file browser.
- The correct version renv should automatically be downloaded and installed upon opening the project (if not already installed). You can then run `renv::restore()` in the console to install the required packages for the given project.

## Assorted notes

- Subdirectories contain a mix of essential and WIP/experimental scripts. Check the **README** in each subdirectory for a description of the scripts required to reproduce analyses in the main manuscript.
- The input dependencies of each script are not well documented (apologies). Generally, should assume that will need to run `preprocessing` scripts first to generate minimal outputs for other steps, and `misc` subdirectory contains scripts that should be run last (to generate supplementary tables). Eventual plan is to produce a Snakemake pipeline to automate running different steps (and document the required running order).
- Minimal data to reproduce analysis is still a work in progress (again apologies). Eventual plan is to integrate with Zenodo (possibly using SciDataFlow).
- Scripts have varying levels of generalisability. If you'd like to adapt some of the code to your inputs, please do get in touch and I'll do my best to give you some advice/help out.