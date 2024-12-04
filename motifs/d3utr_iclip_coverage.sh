#!/usr/bin/env bash

echo 'iCLIP coverage for 3pUTR proximal PAS, SH-SY5Y expressed background'

Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-12-15_papa_cryptic_d3utr.background_shsy5y.pas.proximal.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_shsy5y/d3utr/proximal -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools

echo 'iCLIP coverage for 3pUTR distal PAS, SH-SY5Y expressed background'

Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-12-15_papa_cryptic_d3utr.background_shsy5y.pas.distal.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_shsy5y/d3utr/distal -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools

echo 'iCLIP coverage for shortening 3pUTR proximal PAS, SH-SY5Y expressed background'

Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-12-15_papa_cryptic_d3utr_proximal.background_shsy5y.pas.proximal.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_shsy5y/d3utr_proximal/proximal -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools

echo 'iCLIP coverage for shortening 3pUTR distal PAS, SH-SY5Y expressed background'

Rscript scripts/cl_iclip_coverage.R -r processed/iclip_regions/2023-12-15_papa_cryptic_d3utr_proximal.background_shsy5y.pas.distal.bed -i data/iCLIP/tardbp-shsy5y.concat.sort.chr.bed -c data/GRCh38.primary_assembly.genome.chromsizes.txt -w 500 -o processed/iclip_maps/coverage/background_shsy5y/d3utr_proximal/distal -p /home/sam/mambaforge-pypy3/envs/pybioinfo/bin/bedtools
