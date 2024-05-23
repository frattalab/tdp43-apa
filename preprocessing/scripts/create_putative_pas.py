#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
from typing import Literal
import argparse
import sys


def _df_swap_coord(df: pd.DataFrame,
                   change: Literal['Start', 'End'],
                   replace_col: str
                   ):
    '''Swap one end of a coordinate range (i.e. either start or end) with a joined value

    Note: intended to be applied to internal df of pyranges 0.x/1.x dataframe i.e. all strand values are the same

    Adapted from pr.methods.new_position._new_position (to only swap a single coordinate)

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    change : Literal[&#39;Start&#39;, &#39;End&#39;]
        _description_
    replace_col : str
        _description_
    old_out_suffix : str
        _description_

    Returns
    -------
    _type_
        _description_
    '''    ''''''
    
    # '''
    # Swap values in Start/End coordinates with a provided column
    # Adapted from pr.methods.new_position._new_position (to only swap a single coordinate)
    # '''
    assert isinstance(df, pd.DataFrame)
    assert change in ["Start", "End"]
    assert replace_col in df.columns

    # copy of original    
    to_change = df[change].copy()
    
    df.loc[:, change] = df[replace_col]
    df.loc[:, replace_col] = to_change

    return df


def _df_update_3p(df: pd.DataFrame, replace_suffix: str = "_b"):
    '''_summary_

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    '''

    assert "Strand" in df.columns
    assert (df["Strand"] == "+").all() or (df["Strand"] == "-").all()

    if (df["Strand"] == "+").all():
        # Update End col to update 3'end of interval
        out_col = "End" + replace_suffix
        assert out_col in df.columns
        out = _df_swap_coord(df, "End", out_col)
        

    else:
        # Update Start col to update 3'end of interval
        out_col = "Start" + replace_suffix
        assert out_col in df.columns
        out = _df_swap_coord(df, "Start", out_col)
    
    out_msk = out["End"] - out["Start"] >= 0
    print(f"Number of negative putative intervals (to be dropped) - {len(out_msk)}")
    out = out[out_msk]

    return out



minimal_gtf_cols = ["le_id", "ref_gene_name"]

def main(gtf_path: str,
         pas_bed_path: str,
         max_distance: int = 5000,
         output_prefix: str = "last_exons.updated_pas"):
    '''_summary_

    Parameters
    ----------
    gtf_path : str
        _description_
    bed_path : str
        _description_
    max_distance : int, optional
        _description_, by default 5000
    '''

    # read in PAPA GTF
    print("reading in PAPA GTF...")
    gtf = pr.read_gtf(gtf_path)[minimal_gtf_cols]



    # read in PAS BED
    print("Reading in PAS BED...")
    pas_bed = pr.read_bed(pas_bed_path)

    # Join last exons with all PAS, dropping those that have no matches/cannot be updated
    print(f"Joining last exons to all PAS within specified window (either side) - {max_distance}")
    gtf_upd = gtf.join(pas_bed, strandedness="same", slack=max_distance,how=None, suffix="_pas")

    # print(gtf_upd)

    # Generate all putative updated 3'ends
    gtf_upd = gtf_upd.apply(lambda df: _df_update_3p(df, replace_suffix="_pas"))

    # print(gtf_upd)
    # print(gtf_upd.columns)

    # Remove duplicates
    gtf_upd = gtf_upd.drop_duplicate_positions(strand=True)

    # print(gtf_upd)

    # construct 'updated' identifier, consisting of isoform_id (le_id), pas_id (Name field) and gene_name (why not)
    gtf_upd = gtf_upd.assign("NameUpdated", lambda df: df["le_id"].str.cat(df[["ref_gene_name", "Name"]], sep="|"))
    # print(gtf_upd[["NameUpdated", "le_id"]])

    # get set of all updated IDs
    upd_ids = set(gtf_upd.as_df()["le_id"])

    # Clean up to be ready as a BED file
    gtf_upd = gtf_upd.drop(like="_pas$").drop(["le_id", "ref_gene_name", "Name"])
    gtf_upd = gtf_upd.apply(lambda df: df.rename(columns={"NameUpdated": "Name"}))
    


    gtf_nupd = gtf.subset(lambda df: ~df["le_id"].isin(upd_ids))
    gtf_nupd = gtf_nupd.drop_duplicate_positions()
    nupd_ids = set(gtf_nupd.as_df()["le_id"])

    print(f"Number of IDs that were updated - {len(upd_ids)}")
    print(f"Number of IDs that weren't updated - {len(nupd_ids)}")

    # as before - create updated identifier, adding dummy value for 3rd field/updated PAS
    gtf_nupd = gtf_nupd.assign("NameDummy", lambda df: pd.Series(["not_updated"]*len(df), index=df.index))
    gtf_nupd = gtf_nupd.assign("NameUpdated", lambda df: df["le_id"].str.cat(df[["ref_gene_name", "NameDummy"]], sep="|"))
    gtf_nupd.Score = 0
    gtf_nupd = gtf_nupd.apply(lambda df: df.rename(columns={"NameUpdated": "Name"}))[["Name", "Score"]]

    # print(gtf_nupd)

    gtf_orig = gtf.drop_duplicate_positions()
    gtf_orig = gtf_orig.assign("NameDummy", lambda df: pd.Series(["not_updated"]*len(df), index=df.index))
    gtf_orig = gtf_orig.assign("NameUpdated", lambda df: df["le_id"].str.cat(df[["ref_gene_name", "NameDummy"]], sep="|"))
    gtf_orig.Score = 0
    gtf_orig = gtf_orig.apply(lambda df: df.rename(columns={"NameUpdated": "Name"}))[["Name", "Score"]]

    print("Outputting to file...")
    gtf_upd.to_bed(output_prefix + ".updated.bed")
    gtf_nupd.to_bed(output_prefix + ".not_updated.bed")
    gtf_orig.to_bed(output_prefix + ".original.bed")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Generate all putative updated 3'ends of PAPA predictions with BED file of PAS")
    parser.add_argument('-g', '--gtf_path', type=str, required=True, help="Path to the PAPA novel GTF file")
    parser.add_argument('-b', '--pas_bed_path', type=str, required=True, help="Path to the PAS BED file")
    parser.add_argument('--max_distance', type=int, default=5000, help="Maximum distance between novel PAS and BED PAS to retain as an updated end (default: 5000)")
    parser.add_argument('-o', '--output_prefix', type=str, default="last_exons.updated_pas", help="Output prefix (default: last_exons.updated_pas)")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.gtf_path, args.pas_bed_path, args.max_distance, args.output_prefix)