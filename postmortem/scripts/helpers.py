#!/usr/bin/env python3

from __future__ import print_function
import pyranges as pr
import pandas as pd
import numpy as np
import sys


def eprint(*args, **kwargs):
    '''
    Nice lightweight function to print to STDERR (saves typing, I'm lazy)
    Credit: https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python (MarcH)
    '''
    print(*args, file=sys.stderr, **kwargs)


def _n_ids(gr, id_col, is_df=False):

    assert id_col in gr.columns

    if not is_df:
        return len(set(gr.as_df()[id_col]))
    else:
        return gr[id_col].nunique()


def _pd_merge_gr(df, df_to_merge, how, on, suffixes, to_merge_cols):
    '''
    Perform a pd.merge inside a pr.apply to add columns from gr_to_merge based on metadata in gr_df
    Here, will check the chromosome and strand of provided gr (gr_df)
    and subset gr_to_merge to chr/strand before converting to df
    This should cut down on memory requirements when convert to df (which requires pd.concat() x chrom & strand)
    For this kind of merge, only expect joins between regions on same chr/strand
    '''
    #chromsomes returns a list of chr names (always 1 val)
    assert isinstance(to_merge_cols, list)

    if df_to_merge.empty:
        eprint("df_to_merge for chr/strand pair {} is empty - returning to_merge cols filled with NaNs".format(",".join([df.Chromosome.iloc[0], df.Strand.iloc[0]])))

        df_cols = df.columns.tolist()
        # on won't have suffix added - need to remove as a target column
        to_merge_cols = [col for col in to_merge_cols if col != on]

        # eprint(df_cols)
        # eprint(to_merge_cols)

        # list of cols shared between dfs - need to add suffixes[1]
        # list of cols only in df_to_merge - these stay the same in a merge
        only_cols = [col for col in to_merge_cols if col not in df_cols]
        shared_cols = [col for col in to_merge_cols if col in df_cols]

        # eprint("only_cols - {}".format(only_cols))
        # eprint("shared_cols - {}".format(shared_cols))

        out_shared = [col + suffixes[1] for col in shared_cols]
        target_cols = out_shared + only_cols

        nrows = len(df.index)

        out_cols = {col: pd.Series([np.nan]*nrows) for col in target_cols}

        return df.assign(**out_cols)

    else:
        return df.merge(df_to_merge,
                        how=how,
                        on=on,
                        suffixes=suffixes)


def check_concat(gr):
    '''
    Validate that underlying dfs of (concatenated) PyRanges object have same number of columns

    If dfs do not have same shape, methods like gr.assign will fail due to inconsistent index to add column
    This heuristic will find the first chr/str with the most columns
    and update all other chr/str dfs with missing cols from the max filled wih NaNs (making up to same number of cols)
    '''

    col_lens = {chr_str: len(df.columns) for chr_str, df in gr}

    if len(set(col_lens.values())) == 1:
        # All columns are of same length
        return gr

    else:
        # Some dfs have missing columns
        # Want to make sure each df has same number (and labels) of columns
        default_cols = ["Chromosome", "Start", "End", "Strand"]

        # Get all metadata cols from every df in gr
        # Will be lots of duplicates in here, but this way will avoid massively changing order of columns (from the first df in gr)
        extra_cols = [col for df in gr.values() for col in df.columns if col not in default_cols]

        # Generate list of target cols without duplicates (preserving order of appearance)
        all_cols = list(dict.fromkeys(default_cols + extra_cols))

        # Update all dfs so they have same columns (and order)
        # Col will be filled with NaNs if df didn't have before
        out_gr = gr.apply(lambda df: df.reindex(columns=all_cols))

        return out_gr


def get_internal_regions(gr,
                         feature_col="Feature",
                         feature_key="exon",
                         id_col="transcript_id",
                         region_number_col="exon_number",
                         ):
    '''
    Return gr of internal exons for each transcript_id
    In process, exon_number_col will be converted to type 'int'
    '''

    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)


    # Pull out exons, convert exon_number to int
    exons_gr = gr.assign(region_number_col,
                         lambda df: df[region_number_col].astype(float).astype("Int64"),
                         nb_cpu=1)

    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    exons_gr = exons_gr.apply(lambda df: df.sort_values(by=[id_col,
                                                            region_number_col
                                                            ],
                                                        ascending=True),
                              nb_cpu=1)

    # Filter out 1st + last exons for each ID
    # first exons for each transcript (.ne(1))
    # keep="last" sets last dup value to 'False' & all others True
    # This will filter out last exons

    out_gr = (exons_gr.subset(lambda df: (df[region_number_col].ne(1).astype(bool)) &
                     (df.duplicated(subset=["transcript_id"], keep="last")),
                     nb_cpu=1
                    )
             )

    return out_gr


def get_terminal_regions(gr,
                         feature_col = "Feature",
                         feature_key = "exon",
                         id_col = "transcript_id",
                         region_number_col = "exon_number",
                         source = None,
                         which_region="last",
                         filter_single = False,
                         ):
    '''
    Return gr of last exons for each transcript_id
    In process, region_number_col will be converted to type 'int'
    StringTie merged GTFs (or whatever tool single_steps/stringtie_longreads.smk is using)
    reports exon_number that DOES NOT RESPECT STRAND (from browsing in IGV)
    i.e. for minus-strand - largest exon_number for transcript corresponds to FIRST EXON, not last
    Annotated (i.e. Ensembl) reported exon_numbers DO RESPECT STRAND (i.e. max always = last exon)

    if Do respect strand, put source = None (default)
    if Don't respect strand, put source = "stringtie" (i.e. plus strand = max, minus strand = min)
    '''

    assert source in [None, "stringtie"]
    assert which_region in ["first", "last"]
    assert region_number_col in gr.columns.tolist()
    assert feature_col in gr.columns.tolist()
    assert id_col in gr.columns.tolist()

    # Make sure only 'exon' features are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)

    # Make sure region_number_col is int
    try:
        mod_gr = (gr.assign(region_number_col,
                            lambda df: df[region_number_col].astype(float).astype(int),
                            nb_cpu=1)
                  )
    except KeyError:
        # Currently getting weird KeyError with assign for certain chromosome
        # Mostly non-std chrom names
        # No error if do '.<exon_number>' to assign, but this makes inflexible to colname
        # Also no error if gr -> df assign -> gr
        eprint("pr.assign returned KeyError. Converting {} to int via pandas df conversion".format(region_number_col))

        mod_gr = gr.as_df()
        mod_gr[region_number_col] = mod_gr[region_number_col].astype(float).astype(int)
        mod_gr = pr.PyRanges(mod_gr)


    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    mod_gr = mod_gr.apply(lambda df: df.sort_values(by=[id_col, region_number_col], ascending=True),
                          nb_cpu=1)


    # Filter out single-exon transcripts
    if filter_single:
        eprint("Filtering for multi-exon transcripts...")
        eprint("Before: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))

        # Setting to 'False' marks all duplicates as True (so keep these)
        mod_gr = mod_gr.subset(lambda df: df.duplicated(subset=[id_col], keep=False), nb_cpu=1)

        eprint("After: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))


    if source is None:
        # source = None means that 1 = first region of group regardless of strand
        # Pick last region entry by max region number for each transcript (id_col)
        # Pick first region entry by min region number for each transcript (id_col)

        # keep="last" sets last in ID to 'False' and all others true (negate to keep last only)
        # keep="first" sets first in ID to 'False'

        out_gr = mod_gr.subset(lambda df: ~(df.duplicated(subset=[id_col], keep=which_region)),
                               nb_cpu=1
                              )


    elif source == "stringtie":
        # Numbering Doesn't respect strand
        # Need to flip selecting first/last in group depending on strand
        # minus strand - pick min if Minus strand, max if plus strand

        if which_region == "first":
            # + strand - pick first in group, - strand - pick last in group

            out_gr = (mod_gr.subset(lambda df:
                                    #1. plus strand & first in group/ID
                                    (df["Strand"] == "+") & ~(df.duplicated(subset=[id_col],
                                                                            keep="first")) |
                                    #2. minus strand & last in group/ID
                                    (df["Strand"] == "-") & ~(df.duplicated(subset=[id_col],
                                                                            keep="last")),
                                    nb_cpu=1)
                     )

        elif which_region == "last":
            # + strand - pick last in group/ID
            # - strand - pick first in group/ID
            out_gr = (mod_gr.subset(lambda df:
                                    #1. plus strand & last in group/ID
                                    (df["Strand"] == "+") & ~(df.duplicated(subset=[id_col],
                                                                            keep="last")) |
                                    #2. minus strand & first in group/ID
                                    (df["Strand"] == "-") & ~(df.duplicated(subset=[id_col],
                                                                            keep="first")),
                                    nb_cpu=1)
                     )


    return out_gr


def _df_add_region_number(df,id_col,sort_col="Start"):
    '''
    Return a Series of strand-aware region numbers (5'-3' in 1..n)
    Function to be used internally in a pr.assign (mainly by add_region_number)
    '''
    if id_col not in df.columns.tolist():
        raise KeyError(f"id_col - {id_col} - is not present in df for chr/strand pair {','.join([df.Chromosome.iloc[0], df.Strand.iloc[0]])}")

    elif (df.Strand == "+").all():
        # Start position smallest to largest = 5'-3'

        return df.groupby(id_col)[sort_col].rank(method="min", ascending=True)

    elif (df.Strand == "-").all():
        # Start position largest to smallest = 5'-3'

        return df.groupby(id_col)[sort_col].rank(method="min", ascending=False)

    elif df.empty:
        eprint("df is empty - returning empty pd.Series")
        return pd.Series()


def add_region_number(gr,
                      id_col="transcript_id",
                      feature_key="intron",
                      out_col="intron_number",
                      feature_col="Feature",
                      nb_cpu=1):
    '''
    Adds column to gr containing a strand aware region number column,
    ordered 5'-3' 1..n by a group of features (e.g. transcript)
    '''

    # Make sure only 'feature_key' rows are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)

    # Make sure sorted by position first.
    gr = gr.sort()

    # Add in region number column in strand aware manner, so 1 = most 5', n = most 3'

    gr = gr.assign(out_col, lambda df: _df_add_region_number(df, id_col), nb_cpu=nb_cpu)

    return gr


def cluster_to_region_number(gr, group_id_col, out_col="le_number", cluster_col="Cluster"):
    '''
    Returns gr with 'out_col' column added
    where out_col is leftmost to rightmost cluster_col converted to a
    strand-aware 1..n order by group_id_col
    1 = most 5' site in group_id_col
    n = most 3' site in group_id_col
    '''

    # For groupby.rank to work appropriately, need a single row per 'Cluster'
    c2p = (gr[[group_id_col, cluster_col]]
           .apply(lambda df: df.drop_duplicates(subset=cluster_col)
                  )
           )

    # eprint(c2p)

    # Add 1..n 5'-3' region number as a column
    c2p = c2p.assign(out_col,
                     lambda df: _df_add_region_number(df, group_id_col, cluster_col))

    # eprint(c2p)

    # Return 'out_col' to original gr
    c2p_cols = c2p.columns.tolist()
    out_gr = gr.apply_pair(c2p, lambda df, df_to_merge: _pd_merge_gr(df,
                                                                     df_to_merge,how="left",
                                                                     on=cluster_col,
                                                                     suffixes=[None, "_match"],
                                                                     to_merge_cols=c2p_cols)
                           )

    # eprint(out_gr[["Cluster", "le_number"]])
    # eprint(out_gr.columns)

    # avoid adding extra 'PyRanges' cols (Chromosome etc.) from c2p
    return out_gr.drop(like="_match$")


def annotate_le_ids(combined, id_col="gene_id", le_id_outcol="le_id"):
    '''
    '''

    # Make sure all dfs of gr have same columns (number & labels)
    # This *should* be the case, but if chr/strand df is unique to one of the concatenated grs then can see different num of cols
    combined = check_concat(combined)

    combined = combined.cluster(by=id_col)

    # eprint(combined)

    # Assign 5'-3' 1..n 'last_exon number' for each gene
    # Group together overlapping exons with a common identifier
    # .cluster(strand=None) groups as ('I') expect i.e. only overlapping intervals on the same strand can be merged
    combined = cluster_to_region_number(combined, id_col)

    eprint("assigning 'le_id' (last exon ID) for each gene...")

    # eprint(combined[[id_col, "transcript_id", "Feature", "le_number", "Cluster"]].print(n=20))

    combined = combined.assign(le_id_outcol,
                               lambda df: df[id_col] + "_" + df["le_number"].astype(int).astype(str)
                               )

    return combined



def _df_add_region_rank(df,
                        id_col,
                        region_number_col,
                        first_key,
                        internal_key,
                        last_key):
    '''
    '''

    conditions = [df[region_number_col] == 1,
                  # keep="last" sets last in ID to 'False' and all others true (negate to keep last only)
                  # Safe as gr is sorted by tx_id and region_number_col prior
                  ~df.duplicated(subset=[id_col], keep="last")]

    choices = [first_key, last_key]

    decisions = np.select(conditions, choices, default=internal_key)

    return pd.Series(decisions, index=df.index)


def add_region_rank(gr,
                    id_col="transcript_id",
                    region_number_col="exon_number",
                    out_col="region_rank",
                    first_key="first",
                    internal_key="internal",
                    last_key="last"):
    '''
    Add a column specifying whether region corresponds to the most 3' region in the group (e.g. transcript) or not (1/0)
    '''

    # Sort regions by id_col & 5'-3' position (region_number_col must be strand_aware)

    gr = gr.apply(lambda df: df.sort_values(by=[id_col,
                                                region_number_col],
                                            ascending=True),
                  nb_cpu=1)

    gr = gr.assign(out_col,
                   lambda df: _df_add_region_rank(df,
                                                  id_col,
                                                  region_number_col,
                                                  first_key,
                                                  internal_key,
                                                  last_key)
                   )

    return gr

def _df_collapse_metadata(df, id_col, standard_cols, collapse_cols, collapse_uniq_cols, collapse_sep):
    '''
    Intended to be applied to internal dfs of PyRanges objects
    '''

    found_collapsed = [col for col in collapse_cols if col in df.columns]

    not_found_collapsed = set(collapse_cols) - set(found_collapsed)

    if len(not_found_collapsed) > 0:
        chr_strand = f"{df.Chromosome.drop_duplicates()[0]},{df.Strand.drop_duplicates()[0]}"
        eprint(f"following 'collapse_cols' columns not found in df (chr/strand) - {chr_strand} - {', '.join(not_found_collapsed)}")

    grouped = df.groupby(id_col)

    # Pick first entry for all standard_cols, these should be same for all rows of id_col
    # Leaves a df with id_col values as index
    std_collapsed = grouped[standard_cols].first()

    # For collapse cols, collapse to single row of delimited strings for each column
    # Again leave a df with id_col values as index labels
    clp_collapsed = grouped[found_collapsed].agg(lambda col: collapse_sep.join(col.astype(str)))

    if collapse_uniq_cols is not None:
        # Collapse these cols to single row of delimited strings whilst dropping duplicates
        # Again leave a df with id_col values as index labels
        clp_uniq_collapsed = grouped[collapse_uniq_cols].agg(lambda col: collapse_sep.join(list(dict.fromkeys(col.astype(str)))))

        int_collapsed = clp_collapsed.merge(clp_uniq_collapsed, left_index=True, right_index=True)

        collapsed = std_collapsed.merge(int_collapsed, left_index=True, right_index=True)

    else:
        # combine by id_col
        collapsed = std_collapsed.merge(clp_collapsed, left_index=True, right_index=True)

    return collapsed


def collapse_metadata(gr,
                      id_col="transcript_id",
                      standard_cols=["Chromosome", "Start", "End", "Strand"],
                      collapse_cols=None,
                      collapse_uniq_cols=None,
                      collapse_sep=","):
    '''
    Collapse to a single entry/row per ID entry whilst retaining/collapsing metadata on duplicate rows
    standard_cols: list of column labels that have the same value for all entries of id_col and do not need to be collapsed.
    This is essential for PyRanges standard columns, as you do not want to be changing their dtypes to string. All columns labels in this list retain their dtype, and the first value is retained
    collapse_cols: list of column labels containing metadata you'd like to collapse to a single row (separated by collapse_sep)
        If None, then all columns in gr except for standard_cols, id_col & collapse_uniq_cols will be collapsed
    collapse_uniq_cols: list of column labels containing metadata you'd like to collapse to a single row whilst dropping duplicate values. Values will maintain order of appearance in df
    '''

    assert all([True if col in gr.columns else False for col in standard_cols])

    if collapse_uniq_cols is not None:
        # Here just checking the columns are found in df
        assert all([True if col in gr.columns else False for col in collapse_uniq_cols])

    if collapse_cols is None:
        if collapse_uniq_cols is not None:
            def_cols = standard_cols + [id_col] + collapse_uniq_cols

        else:
            def_cols = standard_cols + [id_col]

        collapse_cols = [col for col in gr.columns if col not in def_cols]

    else:
        assert all([True if col in gr.columns else False for col in collapse_cols])


    return gr.apply(lambda df: _df_collapse_metadata(df,
                                                     id_col,
                                                     standard_cols,
                                                     collapse_cols,
                                                     collapse_uniq_cols,
                                                     collapse_sep
                                                     )
                    )

