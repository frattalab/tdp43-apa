#!/usr/bin/env bash

echo 'iCLIP coverage for bleedthrough exon starts, SH-SY5Y expressed background'

Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-12-14_papa_cryptic_bleedthrough_uniq.background_shsy5y.le_start.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_shsy5y/bleedthrough_uniq/le_start -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools

echo 'iCLIP coverage for bleedthrough pas, SH-SY5Y expressed background'

Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-12-14_papa_cryptic_bleedthrough_uniq.background_shsy5y.pas.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_shsy5y/bleedthrough_uniq/pas -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools

# echo 'iCLIP coverage for bleedthrough exon starts, all background'
# 
# Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-11-24_papa_cryptic_bleedthrough_uniq.background_all.le_start.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_all/bleedthrough_uniq/le_start -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools
# 
# echo 'iCLIP coverage for bleedthrough pas, all background'
# 
# Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-11-24_papa_cryptic_bleedthrough_uniq.background_all.pas.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_all/bleedthrough_uniq/pas -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools