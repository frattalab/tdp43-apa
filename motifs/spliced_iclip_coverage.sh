#!/usr/bin/env bash

echo 'iCLIP coverage for spliced exon starts, SH-SY5Y expressed background'

Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-12-14_papa_cryptic_spliced.background_shsy5y.le_start.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_shsy5y/spliced/le_start -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools

echo 'iCLIP coverage for spliced pas, SH-SY5Y expressed background'

Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-12-14_papa_cryptic_spliced.background_shsy5y.pas.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_shsy5y/spliced/pas -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools
