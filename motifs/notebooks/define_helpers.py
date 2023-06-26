#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np

def n_uniq_coords(df: pd.DataFrame) -> int:
    '''count unique 3'end coordinates (strand aware)

    Assumes df follows pr.PyRanges format (i.e. minimally has Start, End and Strand columns)

    Parameters
    ----------
    df : pd.DataFrame
        _description_

    Returns
    -------
    int
        _description_
    '''


    if df["Strand"].unique()[0] == "+":
        n_coords = df["End"].nunique()
    
    elif df["Strand"].unique()[0] == "-":
        n_coords = df["Start"].nunique()
    
    else:
        raise ValueError("Strand column must contain '+' or '-'")

    return n_coords


def _df_select_rep_site(df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    '''Select a representative polyA site for a given last exon isoform

    intended to be applied to a grouped dataframe (e.g. by "le_id")

    Parameters
    ----------
    df : pd.DataFrame
        _description_

    Returns
    -------
    tuple[pd.DataFrame, str]
        _description_
    '''

    # First select sites of le_id that match atlas site (if applicable)
    if df["atlas_filter"].astype(int).sum() > 0:

        out_df = df[df["atlas_filter"].eq("1")]

        # 2. Within these, see if multiple different 3'ends predicted
        ### TODO: use Name (atlas ID) instead of 3' coordinate to match
        n_atlas = out_df["Name"].nunique()

        if n_atlas > 1:
                # multiple different atlas sites predicted
                # see if somes sites were called in > 1 dataset - select those with max number of datasets called
                end3_c = out_df.groupby("Name")["experiment_id"].nunique()
                # have a series with Name (PolyASite ID) as index labels, count of experimnts as values
                # get the IDs with largest number of unqiue experiments (as as series)
                end3_c = end3_c[end3_c == end3_c.max()]

                if len(end3_c) > 1:
                    # multiple 3'ends/ atlas sites predicted in same number of datasets
                    # raise Exception(f"{out_df['le_id_quant'].drop_duplicates().iloc[0]} - le_id has multiple predicted atlas sites")
                    # prefer the shortest extension in this case (conservative, defines smaller window in)

                    # calculate length
                    out_df = out_df[out_df["Name"].isin(end3_c.index)]
                    out_df["reg_len"] = out_df.End - out_df.Start

                    # select shortest extension
                    out_df = out_df[out_df.reg_len == out_df.reg_len.min()]
                    # Select a single row for le_id corresponding to selected PAS
                    out_df = out_df.sort_values(by="experiment_id").drop_duplicates(subset=["Start", "End", "Strand"]).drop(columns="reg_len")

                    return out_df, "atlas_max_datasets_shortest"
                                    
                else:
                    # 1 pas remaining
                    end3_coord = end3_c.index[0]

                    # Select a single row for le_id corresponding to selected PAS
                    out_df = out_df[out_df["Name"] == end3_coord].sort_values(by="experiment_id").drop_duplicates(subset=["Start", "End", "Strand"])

                    decision_str = "atlas_max_datasets"


        else:
            # only single atlas site predicted - just select a single row
            end3_coord = out_df["Name"].drop_duplicates().index[0]
            # Select a single row for le_id corresponding to selected PAS
            out_df = out_df[out_df["Name"] == end3_coord].sort_values(by="experiment_id").drop_duplicates(subset=["Start", "End", "Strand"])

            decision_str = "atlas_1_pred"
        

        return out_df, decision_str

    else:
        # TODO: implement filtering for events called with motifs only
        out_df = df[df["atlas_filter"].eq("0")]

        # select predicted site with motif that is closest to expected distance (21nt upstream)
        # -74_TATAAA,-14_TATAAA
        # [-74_TATAAA, -14_TATAAA]
        # [-74, -14]
        # absolute distance from expected position (21nt upstream)
        # select min deviating motif for each predicted PAS
        pas_m = out_df["pas_motifs"].apply(lambda x: min([abs(21 - abs(int(i.split("_")[0]))) for i in x.split(",")]))

        # select predicted pas with minimum deviation from 21nt
        pas_m = pas_m[pas_m == pas_m.min()]

        if len(pas_m) > 1:
            # multiple predicted PAS with motif same distance from expected

            # check if exactly the same coordinate predicted for these events

            # (if true, just sort alphabetically by experiment id and pick first)
            # otherwise, select shortest
            n_pas = n_uniq_coords(out_df)

            if n_pas > 1:
                # select shortest extension
                # calculate length
                out_df = out_df.loc[pas_m.index, :]

                out_df["reg_len"] = out_df.End - out_df.Start

                # select shortest extension
                out_df = out_df[out_df.reg_len == out_df.reg_len.min()]
                # Select a single row for le_id corresponding to selected PAS
                out_df = out_df.sort_values(by="experiment_id").drop_duplicates(subset=["Start", "End", "Strand"]).drop(columns="reg_len")

            else:
                # just a single coord predicted in multiple datasets - 
                out_df = out_df.sort_values(by="experiment_id").drop_duplicates(subset=["Start", "End", "Strand"])

            return out_df, "motif_shortest_min"
        
        else:
            out_df = out_df.loc[pas_m.index, :]
            

            decision_str = "motif_1_min"


        return out_df, decision_str
    

def select_rep_site(gr: pr.PyRanges, id_col: str = "le_id") -> tuple[pr.PyRanges, dict]:
    '''_summary_

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    id_col : str, optional
        _description_, by default "le_id"

    Returns
    -------
    tuple[pr.PyRanges, dict]
        _description_
    '''


    # store selected dfs for each isoform
    sel_dfs = []
    # store selection decision for each isoform
    sel_dict = {"atlas_1_pred": [],
                  "atlas_max_datasets": [],
                  "atlas_max_datasets_shortest": [],
                  "motif_1_min": [],
                  "motif_shortest_min": []}

        
    # select representative pas according to prioritisation scheme
    for le_id, df in gr.as_df().groupby(id_col):
        out_df, decision_str = _df_select_rep_site(df)
        sel_dfs.append(out_df)
        sel_dict[decision_str].append(le_id)

    # combine per-id dfs and output pyranges object
    sel_df = pd.concat(sel_dfs)  
    out_gr = pr.PyRanges(sel_df)     

    return out_gr, sel_dict
          
    
