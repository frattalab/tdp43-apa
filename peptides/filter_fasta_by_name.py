
from Bio import SeqIO
import pandas as pd
import gzip
import sys


def read_fasta(file):
    with gzip.open(file, 'rt') as handle:
        sequences = {record.id.split("|")[1]: str(record.seq) for record in SeqIO.parse(handle, 'fasta') 
                     if record.id.split("|")[1] in ids}
    return sequences

fa = sys.argv[1]
pred_peps = sys.argv[2]
# USP31, CNPY3, ANKRD27
ids = ["ENST00000219689.12", "ENST00000372836.5", "ENST00000587352.5"]

pred = pd.read_csv(pred_peps, sep="\t")
pred["transcript_id"] = pred.cryptic_transcript_id.str.split(";", expand=True)[0]
print(pred)



# Usage

sequences = read_fasta(fa)

# print(sequences)

seq_df = pd.DataFrame.from_dict(sequences,orient="index", columns=["peptide_seq_full"]).reset_index().rename(columns={"index": "transcript_id"})

seq_df = seq_df.merge(pred, how="left", on="transcript_id")
print(seq_df)

# for _, row in seq_df.iterrows():
#     print(row["transcript_id"])
#     print(row["peptide_seq_full_x"])
#     print(row["peptide_seq_full_y"])

# do they match?
print(seq_df["peptide_seq_full_x"] == seq_df["peptide_seq_full_y"].str.rstrip("*"))