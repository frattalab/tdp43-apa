# SLAM-seq analysis

Scripts to run GrandR and visualise stability of cryptic 3'UTR extension-containing genes.

## estimate_half_life.R

All-in-one script to run GrandR on Grand-SLAM derived counts and plot estimated half lives for selected genes (Fig 3).

Requires TSV of Grand-SLAM inferred counts, available at `data/grandslam_exp1_grandR.tsv.gz`. You will need to decompress prior to running the script e.g.

```bash
gzip -d -k data/grandslam_exp1_grandR.tsv.gz
```
