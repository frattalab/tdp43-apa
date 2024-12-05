#!/usr/bin/env python3

# %%
import pyranges as pr
import pandas as pd


# Load in cryptic information
cryptics_df = pd.read_csv("../data/2023-12-10_cryptics_summary_all_events_bleedthrough_manual_validation.tsv", sep="\t")
cryptic_le_ids = set(cryptics_df.le_id)
papa_gtf = pr.read_gtf("../data/PAPA/novel_ref_combined.last_exons.gtf").subset(lambda df: df.le_id.isin(cryptic_le_ids))
papa_gtf

# %%
# Need to add 'event_type_simple' to attribute field (for generate event-specific decoys regions)
papa_gtf = papa_gtf.apply(lambda df: df.merge(cryptics_df[["le_id", "simple_event_type", "annot_status"]],
                                   on="le_id",
                                   how="left",
                                   suffixes=[None, "_simple"]).rename(columns={"simple_event_type": "event_type_simple"}))


# %%
# double check event type assignment
papa_gtf.as_df().drop_duplicates(subset=["le_id", "event_type_simple"]).event_type_simple.value_counts()

# %%
# output to GTF
papa_gtf.to_gtf("../processed/2024-09-27_cryptic_last_exons.gtf")
# %%
