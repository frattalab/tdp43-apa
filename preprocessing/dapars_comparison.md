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