# DaPars comparison for 3'Ext and 3'UTR shortening events


## Datasets

Downloaded RefSeq transcript annotations from the UCSC genome browser FTP website (NCBI RefSeq uses different chromosome naming convention). Used the 110 release archive (believe it's the last using hg38?)

```bash
# assuming already in tdp43-apa/preprocessing
cd data/dapars_comparison
curl https://hgdownload.soe.ucsc.edu/goldenPath/archive/hg38/ncbiRefSeq/110/hg38.110.ncbiRefSeq.gtf.gz -o hg38.110.ncbiRefSeq.gtf.gz
```

## Preprocessing step - filtering GTF to target genes

Filter GTF for specified gene names - will first test out code using ELK1, SIX3 and TLX1

```bash
python scripts/filter_gtf_by_attr.py -v ELK1 SIX3 TLX1 -a gene_name data/dapars_comparison/hg38.110.ncbiRefSeq.gtf.gz data/dapars_comparison/hg38.110.ncbiRefSeq.filtered.elk1_six3_tlx1.gtf
```

## Constructing extended transcript models

Idea:

- Extend existing transcript models to the inferred cryptic 3'Ext
- Output 'proximal', original annotated PAS (and distal PAS) alongside the updated models transcript
- Use GTF as input to DaPars2 snakemake pipeline, take BED file of predicted PAS output by pipeline and compare to annotated

### Constructing extended transcript models for ELK1, SIX3 and TLX1

```bash
python scripts/get_dapars_ext_tx.py -g data/dapars_comparison/hg38.110.ncbiRefSeq.filtered.elk1_six3_tlx1.gtf -b ../motifs/processed/iclip_regions/2023-12-15_3ext.last_exons.cryptic.bed -o processed/dapars_comparison/hg38.110.ncbiRefSeq.filtered.elk1_six3_tlx1
```

## Running DaPars2

Used the APAeval snakemake pipeline to do this, using Seddighi i3 cortical samples. Ran separately using original annotated transcript models and the extended transcript models. The APAeval pipeline parses the DaPars2 output into BED files. I used the '01.bed' output files for each sample for futher analysis (single nucleotide PAS predictions, score field filled with placeholder value ('.')).

### Combining predictions across replicates (condition and annotation-wise)

Wrote a simple script to combine the BED files. Duplicate predictions are collapsed, arbitrarily selecting a representative interval (will only differ by 'Name' field, which contains the corresponding transcript ID for predictions). Origin of intervals/predictions is not tracked/reported. Commands used reported below (executed from `preprocessing` subdirectory with general conda environment active)

```bash
$ python scripts/combine_dapars_beds.py -i data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43*01.bed -o data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.exts.kd.bed 
Provided input BED files:  ['data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065403_S23_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065405_S46_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065407_S25_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065409_S29_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065411_S54_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/exts_elk1_six3_tlx1/TDP43_19065413_S19_DaPars2_01.bed']
Successfully merged 6 BED files to: data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.exts.kd.bed
$ python scripts/combine_dapars_beds.py -i data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43*01.bed -o data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.orig.kd.bed 
Provided input BED files:  ['data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065403_S23_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065405_S46_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065407_S25_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065409_S29_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065411_S54_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/TDP43_19065413_S19_DaPars2_01.bed']
Successfully merged 6 BED files to: data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.orig.kd.bed
$ python scripts/combine_dapars_beds.py -i data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT*01.bed -o data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.orig.ctrl.bed 
Provided input BED files:  ['data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074709_S20_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074711_S17_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074713_S59_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074715_S37_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074717_S21_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074719_S22_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074721_S50_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074723_S36_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074725_S27_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074727_S44_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074729_S12_DaPars2_01.bed', 'data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/orig_elk1_six3_tlx1/NT_19074731_S13_DaPars2_01.bed']
Successfully merged 12 BED files to: data/dapars_comparison/apaeval_dapars2_snakemake/seddighi_i3_cortical/merged.orig.ctrl.bed
```
