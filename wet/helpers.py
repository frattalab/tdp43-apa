import pyranges as pr
import pandas as pd
from typing import Literal


def get_terminal_regions(gr: pr.PyRanges,
                         feature_col = "Feature",
                         feature_key = "exon",
                         id_col = "transcript_id",
                         region_number_col = "exon_number",
                         number_type: Literal["stranded", "unstranded"] = "stranded",
                         which_region: str = "last",
                         filter_single = False,
                         ):
    '''Return the first/last interval in group of intervals

    Requires a column that provides a 1..n numbering of intervals within each group (can be generated by add_region_number). Extraction will always be with respect to strand (can handle strand-aware/non-strand aware ranking in region_number_col) 

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    feature_col : str, optional
        _description_, by default "Feature"
    feature_key : str, optional
        _description_, by default "exon"
    id_col : str, optional
        _description_, by default "transcript_id"
    region_number_col : str, optional
        _description_, by default "exon_number"
    number_type : Literal[&quot;stranded&quot;, &quot;unstranded&quot;], optional
        _description_, by default "stranded"
    which_region : str, optional
        _description_, by default "last"
    filter_single : bool, optional
        _description_, by default False

    Returns
    -------
    _type_
        _description_
    '''

    assert number_type in ["stranded", "unstranded"]
    assert which_region in ["first", "last"]
    assert region_number_col in gr.columns.tolist()
    assert feature_col in gr.columns.tolist()
    assert id_col in gr.columns.tolist()

    # Make sure only 'exon' features are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)
    # Make sure region_number_col is int
    # assert gr.as_df()[region_number_col].dtype

    # TODO: Make this into a separate function
    # Make sure region_number_col is int
    # try:
    #     mod_gr = (gr.assign(region_number_col,
    #                         lambda df: df[region_number_col].astype(float).astype(int),
    #                         nb_cpu=1)
    #               )
    # except KeyError:
    #     # Currently getting weird KeyError with assign for certain chromosome
    #     # Mostly non-std chrom names
    #     # No error if do '.<exon_number>' to assign, but this makes inflexible to colname
    #     # Also no error if gr -> df assign -> gr
    #     print("pr.assign returned KeyError. Converting {} to int via pandas df conversion".format(region_number_col))

    #     mod_gr = gr.as_df()
    #     mod_gr[region_number_col] = mod_gr[region_number_col].astype(float).astype(int)
    #     mod_gr = pr.PyRanges(mod_gr)


    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    mod_gr = gr.apply(lambda df: df.sort_values(by=[id_col, region_number_col], ascending=True),
                          nb_cpu=1)


    # Filter out single-exon transcripts
    if filter_single:
        print("Filtering for multi-exon transcripts...")
        print("Before: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))

        # Setting to 'False' marks all duplicates as True (so keep these)
        mod_gr = mod_gr.subset(lambda df: df.duplicated(subset=[id_col], keep=False), nb_cpu=1)

        print("After: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))


    if number_type == "stranded":
        # source = None means that 1 = first region of group regardless of strand
        # Pick last region entry by max region number for each transcript (id_col)
        # Pick first region entry by min region number for each transcript (id_col)

        # keep="last" sets last in ID to 'False' and all others true (negate to keep last only)
        # keep="first" sets first in ID to 'False'

        out_gr = mod_gr.subset(lambda df: ~(df.duplicated(subset=[id_col], keep=which_region)),
                               nb_cpu=1
                              )


    else:
        # Numbering doesn't respect strand
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


def _df_collapse_metadata(df, id_col, standard_cols, collapse_cols, collapse_uniq_cols, collapse_sep):
    '''
    Intended to be applied to internal dfs of PyRanges objects
    '''

    found_collapsed = [col for col in collapse_cols if col in df.columns]

    not_found_collapsed = set(collapse_cols) - set(found_collapsed)

    if len(not_found_collapsed) > 0:
        chr_strand = f"{df.Chromosome.drop_duplicates()[0]},{df.Strand.drop_duplicates()[0]}"
        print(f"following 'collapse_cols' columns not found in df (chr/strand) - {chr_strand} - {', '.join(not_found_collapsed)}")

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
                                                     ))