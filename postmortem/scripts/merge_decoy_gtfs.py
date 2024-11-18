import argparse
import pyranges as pr
import sys

def main(non_cryptics_path: str, decoys_path: str, cryptic_others_path: str, output_prefix: str,
         event_tx_col: str = 'transcript_id', event_id_col: str = 'le_id', gene_id_col: str = 'gene_id', gene_name_col: str = 'gene_name') -> None:
    # Read the GTF files
    non_cryptics = pr.read_gtf(non_cryptics_path)
    decoys = pr.read_gtf(decoys_path)
    cryptic_others = pr.read_gtf(cryptic_others_path)

    # Combine the GTF files into a single object
    combined = pr.concat([non_cryptics, decoys, cryptic_others])

    # Remove duplicate entries
    combined = combined.drop_duplicate_positions(strand=True)

    # Output the merged GTF to a file
    print(f"Writing combined GTF to file - {output_prefix + '.merged.gtf'}")
    combined.to_gtf(output_prefix + ".merged.gtf")

    # Produce 'tx2le', 'tx2gene', 'le2gene', and 'le2name' files
    print(f"Writing 'tx2le' (transcript_id | le_id) to TSV... - {output_prefix + '.tx2le.tsv'}")
    (combined.as_df()
     [[event_tx_col, event_id_col]]
     .rename(columns={event_tx_col: "transcript_id", event_id_col: "le_id"})
     .drop_duplicates()
     .sort_values(by="transcript_id")
     .to_csv(output_prefix + ".tx2le.tsv",
             sep="\t",
             index=False,
             header=True)
     )

    print(f"Writing 'tx2gene' (transcript_id | gene_id) to TSV... - {output_prefix + '.tx2gene.tsv'}")
    (combined.as_df()
     [[event_tx_col, gene_id_col]]
     .rename(columns={event_tx_col: "transcript_id", gene_id_col: "gene_id"})
     .drop_duplicates()
     .sort_values(by="transcript_id")
     .to_csv(output_prefix + ".tx2gene.tsv",
             sep="\t",
             index=False,
             header=True)
     )

    print(f"Writing 'le2gene' (le_id | gene_id) to TSV... - {output_prefix + '.le2gene.tsv'}")
    (combined.as_df()
     [[event_id_col, gene_id_col]]
     .rename(columns={event_id_col: "le_id", gene_id_col: "gene_id"})
     .drop_duplicates()
     .sort_values(by="le_id")
     .to_csv(output_prefix + ".le2gene.tsv",
             sep="\t",
             index=False,
             header=True)
     )

    print(f"Writing 'le2name' (le_id | gene_name) to TSV... - {output_prefix + '.le2name.tsv'}")
    (combined.as_df()
     [[event_id_col, gene_name_col]]
     .rename(columns={event_id_col: "le_id", gene_name_col: "gene_name"})
     .drop_duplicates()
     .dropna()
     .sort_values(by="le_id")
     .to_csv(output_prefix + ".le2name.tsv",
             sep="\t",
             index=False,
             header=True)
     )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge decoy-augmented cryptic GTF with original quantification GTF for downstream analysis")
    parser.add_argument("--non_cryptics", required=True, help="Path to the non-cryptics GTF file.")
    parser.add_argument("--decoys", required=True, help="Path to the decoys GTF file.")
    parser.add_argument("--cryptic_others", required=True, help="Path to the GTF file containing non-cryptic isoforms from cryptic genes")
    parser.add_argument("--output_prefix", required=True, help="Prefix for the output files.")
    parser.add_argument("--event_tx_col", default="transcript_id", help="Column name for event transcript ID.")
    parser.add_argument("--event_id_col", default="le_id", help="Column name for event ID.")
    parser.add_argument("--gene_id_col", default="ref_gene_id", help="Column name for gene ID.")
    parser.add_argument("--gene_name_col", default="ref_gene_name", help="Column name for gene name.")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    
    main(args.non_cryptics, args.decoys, args.cryptic_others, args.output_prefix,
         args.event_tx_col, args.event_id_col, args.gene_id_col, args.gene_name_col)