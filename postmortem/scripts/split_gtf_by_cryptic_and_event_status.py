#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import sys

help = """Usage: python split_gtf_by_cryptic_event_status.py GTF CRYPTICS_SUMMARY_DF OUTPUT_PREFIX
GTF - input PAPA GTF
CRYPTICS_SUMMARY_DF - Summary df containing cryptic IDs and event_type annotations
OUTPUT_PREFIX - <output_prefix>.<cryptic_status>.<event_type>.<ids_or_genes>.gtf
"""

if len(sys.argv) == 1 or any(h in sys.argv for h in ["-h", "--help"]):
    print(help)
    sys.exit(0)

print("Reading in input GTF...")
in_gtf = pr.read_gtf(sys.argv[1])
cryptics_df = pd.read_table(sys.argv[2], sep="\t")
output_prefix = sys.argv[3]

# subset for cryptic containing genes
print("Filtering input GTF for cryptics/non-cryptics")
cryptics_df.loc[:, "gene_id"] = cryptics_df.le_id.str.split("_", expand=True)[0]
cryptic_genes = set(cryptics_df.gene_id)
gtf_cryptic = in_gtf.subset(lambda df: df.ref_gene_id.isin(cryptic_genes))
gtf_other = in_gtf.subset(lambda df: ~df.ref_gene_id.isin(cryptic_genes))

# subset cryptic for different event types
ipa_le_ids = set(cryptics_df.loc[cryptics_df.simple_event_type == "bleedthrough", "le_id"])
ale_le_ids = set(cryptics_df.loc[cryptics_df.simple_event_type == "spliced", "le_id"])
ext3_le_ids = set(cryptics_df.loc[cryptics_df.simple_event_type == "distal_3utr_extension", "le_id"])
proxext3_le_ids = set(cryptics_df.loc[cryptics_df.simple_event_type == "proximal_3utr_extension", "le_id"])
complex_le_ids = set(cryptics_df.loc[cryptics_df.simple_event_type.str.contains(",", regex=False), "le_id"])

print("Outputting non-cryptics GTF")
# output non-cryptics to a single GTF)
gtf_other.to_gtf(output_prefix + ".non_cryptics.all.all.gtf")

print("Outputting cryptic le_ids to GTF (split by event type, 1 per file)")
# Output per-event type GTFs of just the cryptic IDs
for id_key, le_ids in zip(["ipa", "ale", "ext3", "proxext3", "complex"],
                    [ipa_le_ids, ale_le_ids, ext3_le_ids, proxext3_le_ids, complex_le_ids]):
    gtf_cryptic.subset(lambda df: df.le_id.isin(le_ids)).to_gtf(".".join([output_prefix, "cryptics", id_key, "ids", "gtf"]))

print("Outputting GTF of cryptic-containing genes (with cryptic IDs excluded)")
#  Output GTF of cryptic genes (with the cryptic IDs excluded)
gtf_cryptic.subset(lambda df: df.le_id.isin(set(cryptics_df.le_id))).to_gtf(output_prefix + ".cryptics.all.non_cryptic_ids.gtf")



