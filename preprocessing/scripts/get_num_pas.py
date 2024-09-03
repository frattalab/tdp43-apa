


# %%
import pyranges as pr
import pandas as pd
import numpy as np
import os

gtf = pr.read_gtf("data/novel_ref_combined.last_exons.gtf")
# df containing cryptic events across any datasets
cryptics_df = pd.read_csv("processed/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv", sep="\t")
# combined DEXSeq output of all datasets (used to get cleaned event type and annotation status)
dexseq_df = pd.read_csv("processed/PAPA/2023-12-10_i3_cortical_zanovello.all_datasets.dexseq_apa.results.processed.cleaned.tsv",
                        sep="\t", 
                        usecols=["le_id", "simple_event_type", "annot_status"])
outdir = "processed/curation/cryptic_annot_comparison"

# duplicated once per dataset (remove)
dexseq_df.drop_duplicates(inplace=True)

gtf

# %%
# subset GTF to minimal cols

# Extract all PAS, drop duplicate coordinates for each le_id
gtf = gtf[["le_id", "transcript_id", "ref_gene_id", "ref_gene_name", "event_type"]].three_end()
gtf = gtf.apply(lambda df: df.drop_duplicates(subset=["Start", "End", "Strand", "le_id"]))

# All remaining entries should be unique PAS coordinates - count per le_id
le_pas_counts = gtf.le_id.value_counts().reset_index()
le_pas_counts

# %%

# Repeat unique PAS counting, instead collapsing closely spaced PAS to a single interval
gtf_clpsd = gtf.cluster(strand=True, by="le_id", slack=12).apply(lambda df: df.drop_duplicates(subset=["Cluster"]))
le_pas_counts_clpsd = gtf_clpsd.le_id.value_counts().reset_index()
le_pas_counts_clpsd

# %%
# combine PAS count dfs
le_pas_counts = le_pas_counts.merge(le_pas_counts_clpsd, on="le_id", suffixes=[None, "_clpsd"])
le_pas_counts
# %%
# add annotation information to df (first collapsing duplicated/distinct annots to a single string)
dexseq_df = dexseq_df.groupby("le_id").agg(lambda x: ','.join(sorted(set(x))))
le_pas_counts = dexseq_df.merge(le_pas_counts, on="le_id")
le_pas_counts

# %%
# add cryptic info
le_pas_counts.loc[:, "cryptic_status"] = np.where(le_pas_counts.le_id.isin(set(cryptics_df.le_id)), 1, 0)
le_pas_counts = le_pas_counts.reindex(columns=["le_id", "annot_status", "cryptic_status", "simple_event_type", "count", "count_clpsd12"])
le_pas_counts

# %%
# double-check le_ids only represented once
print(le_pas_counts.le_id.value_counts().describe())

# %%
# output to disk
if not os.path.exists(outdir):
        os.makedirs(outdir)

le_pas_counts.to_csv(os.path.join(outdir, "2024-09-03_le_id_pas_counts.tsv"), sep="\t", index=False, header=True)


# %%