import pyranges as pr
import pandas as pd


original_bed = pr.read_bed("processed/curation/2024-05-23_last_exons.max_distance_10000.rep.original.bed")
updated_bed = pr.read_bed("processed/curation/2024-05-23_last_exons.max_distance_10000.rep.updated.bed")

original_bed = original_bed.three_end()
updated_bed = updated_bed.three_end()

comb_bed = pr.concat([original_bed, updated_bed])
print(len(comb_bed))
print(comb_bed)
# original intervals have score 0, joined = non-zero
# sort ascending
comb_bed_reduced = comb_bed.apply(lambda df: df.sort_values(by=["Start", "End", "Score"], ascending=True)).drop_duplicate_positions(keep="first")
print(len(comb_bed_reduced))
print(comb_bed_reduced)
comb_bed_reduced = comb_bed_reduced.sort()
print(comb_bed_reduced)

comb_bed_reduced.to_bed("processed/curation/2024-05-23_pas.max_distance_10000.original_updated_combined.bed")
