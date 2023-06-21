#!/usr/bin/env python3

from __future__ import print_function
import pyranges as pr
import pandas as pd
import numpy as np
from helpers import *
import sys
import argparse

'''
Quick and dirty script to augment last exon GTF with dummy introns/exons for quantification
Input:
- Reference GTF
- Input last exons GTF
- Output prefix

Outputs:
- GTF with last exons + overlapping ref exons/introns
- tx2le, le2gene etc. assignment TSVs excluding the dummy introns/exons
'''

# sin3b_tx = "ENST00000596802.5"

def _n_ids(gr, id_col, is_df=False):

    assert id_col in gr.columns

    if not is_df:
        return len(set(gr.as_df()[id_col]))
    else:
        return gr[id_col].nunique()


def _df_grp_update_le_number(df,
                             le_num_col="le_number",
                             is_ext_col="is_extension"):
    '''
    # TODO: update to cumulative sum based approach
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.cumsum.html?highlight=cumsum
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.core.groupby.DataFrameGroupBy.cumsum.html?highlight=cumsum
    - Basically be, cumulative sum of 'is_extension' (modified so same le_id = 1 extension only (1st))
    - Should have same index, so add df.groupby(id_col).apply(lambda df: df[le_number].add(cumsum))
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.add.html?highlight=add#pandas.Series.add

    Returns: pd.Series
    '''

    # Track how many extension IDs for gene starting from 5'end
    # Once encounter an ext, all IDs downstream need to be shifted down by n of preceeding extensions
    ext_count = 0

    # Track le_id of most recently found ext event
    # If gene has alt extension isoforms of the same last exon
    # Will group together for simplification's sake
    prev_ext_le_id = ""

    out = []

    for _, row in df.iterrows():
        if row[is_ext_col] == 1:
            # Is an extension
            # Check if same LE as before
            if row[le_num_col] == prev_ext_le_id:
                # Tx is a diff extension of same le - group together
                out.append(row[le_num_col] + ext_count)

            else:
                # New distinct extension event
                ext_count += 1
                out.append(row[le_num_col] + ext_count)

                # Track le_id in case another for the same last exon
                prev_ext_le_id = row[le_num_col]

        elif row[is_ext_col] == 0:
            # Not an extension
            # If no previous extension, le number is unchaged
            # Otherwise shifted by n extensions
            out.append(row[le_num_col] + ext_count)

    return pd.Series(out, index=df.index)


def update_extension_le_number(df, id_col="gene_id", out_col="le_number_ext"):
    '''
    '''

    # assert

    df = (df.groupby(id_col)
          .apply(lambda grp: grp.assign(**{out_col: _df_grp_update_le_number(grp)
                                           }
                                        )
                 )
          )


    return df


def update_ext_le_nums(gr,
                      type_col="event_type",
                      le_ext_key="last_exon_extension",
                      id_col="gene_id",
                      number_col="le_number"
                      ):
    '''
    '''

    gr = gr.assign("is_extension",
                   lambda df: pd.Series(np.where(df[type_col].str.contains(le_ext_key,
                                                                           na=False,
                                                                           regex=False),
                                                 1,
                                                 0
                                                 )
                                        )
                   )


    # Extra sort so extension events come after non-extensions of the same le_id
    gr = gr.apply(lambda df:
                  df.sort_values(by=[id_col, number_col, "is_extension"])
                  )

    # Adds column 'le_number_ext'
    gr = gr.apply(lambda df: update_extension_le_number(df))

    # eprint(combined_ext[["ref_gene_id", "le_id", "le_number", "event_type", "is_extension", "le_number_ext"]].print(n=50))

    # Reassign le_number (number_col) with updated number, return to same shape as input gr
    gr = (gr.drop([number_col, "is_extension"])
            .apply(lambda df: df.rename(columns={"le_number_ext": number_col}))
          )

    return gr




def annotate_le_ids_ext(combined, id_col="gene_id", event_type_col="event_type", ext_key="last_exon_extension", le_id_outcol="le_id"):
    '''
    '''

    # Make sure all dfs of gr have same columns (number & labels)
    # This *should* be the case, but if chr/strand df is unique to one of the concatenated grs then can see different num of cols
    combined = check_concat(combined)

    # Identify last exon Ids containing novel extensions - these need to be handled separately
    # set of cluster IDs for each chr, strand tuple
    # {(chr, strand): {cluster_IDs}}
    et_mask = combined.as_df()[event_type_col].str.contains(ext_key, regex=False)
    # eprint(et_mask)
    # eprint(et_mask[et_mask.isna()])
    # eprint(combined.as_df().loc[et_mask[et_mask.isna()].index, ["gene_id", "gene_name", "transcript_id", "event_type","annot_status","event_type_temp"]])
    # eprint(combined.as_df().loc[et_mask[et_mask.isna().index], :])

    d_ext_gene_ids = combined.apply(lambda df: set(df[df[event_type_col].str.contains(ext_key, regex=False)][id_col]) # should be ref gene ID in novel events (found)
                                    ,
                                    as_pyranges=False)

    # Combine into a single set
    # https://blog.finxter.com/union-multiple-sets-in-python/
    ext_gene_ids = set().union(*d_ext_gene_ids.values())

    # Separate objects for extension-containing & non-extension containing genes
    combined_ext = combined.subset(lambda df: df[id_col].isin(ext_gene_ids))
    combined_n_ext = combined.subset(lambda df: ~df[id_col].isin(ext_gene_ids))

    if len(combined_ext) > 0:
        # Are extensions, need to update as separate le numbers

        # First steps are the same - group together overlapping with common identifier
        # Then convert to gene-wise 5'-3' last exon number
        combined_ext = combined_ext.cluster(by=id_col)
        combined_ext = cluster_to_region_number(combined_ext, id_col)

        # Update le_number for extension event
        combined_ext = update_ext_le_nums(combined_ext,
                                          type_col=event_type_col,
                                          le_ext_key=ext_key,
                                          id_col=id_col,
                                          )

        combined_ext = combined_ext.assign(le_id_outcol,
                                           lambda df: df[id_col] + "_" + df["le_number"].astype(int).astype(str)
                                           )

    # Assign last exon numbers for all genes without novel 3'UTR extensions

    # Assign 5'-3' 1..n 'last_exon number' for each gene
    # Group together overlapping exons with a common identifier
    # .cluster(strand=None) groups as ('I') expect i.e. only overlapping intervals on the same strand can be merged
    combined_n_ext = combined_n_ext.cluster(by=id_col)
    combined_n_ext = cluster_to_region_number(combined_n_ext, id_col)

    combined_n_ext = combined_n_ext.assign(le_id_outcol,
                                           lambda df: df[id_col] + "_" + df["le_number"].astype(int).astype(str)
                                           )

    if len(combined_ext) > 0:
        # GTF containing defined last exons (with le_id etc. defined)
        combined_f = pr.concat([combined_ext, combined_n_ext])

    else:
        combined_f = combined_n_ext


    return combined_f


def _df_match_5p_adj(df, suffix):
    '''
    Subset of rows where 5'end is bookended to overlapping feature
    i.e. 5'end of main feature is directly adjacent to 3'end of everlapping feature

    Match directly adjacent/bookended features by the 5'end of the target feature
    '''

    if (df.Strand == "+").all():
        # Compare Start to suffixed end
        decisions = np.where(df['Start'] == df['End' + suffix],
                             True,
                             False)

    elif (df.Strand == "-").all():
        # End (5' end of feature) to Start (3' end of overlapping feature)
        decisions = np.where(df['End'] == df['Start' + suffix],
                             True,
                             False)

    return pd.Series(decisions, index=df.index)


def _df_match_3p_adj(df, suffix):
    '''
    Match directly adjacent/bookended features by the 3'end of the target feature
    '''

    if (df.Strand == "+").all():
        # Compare End (3'end of feature) to Start (5' end of overlapping feature)
        decisions = np.where(df['End'] == df['Start' + suffix],
                             True,
                             False)

    elif (df.Strand == "-").all():
        # Compare Start (3'end of feature) to End (5' end of overlapping feature)
        decisions = np.where(df['Start'] == df['End' + suffix],
                             True,
                             False)

    return pd.Series(decisions, index=df.index)


def _collapse_straddled_exons(df, id_col, suffix):
    '''
    # # for each overlapping intron, subdivide into bookended intervals
    # # find a set of transcript IDs to represent each individual interval
    # # If duplicate tx for a given interval , sort alphabetically then take the first tx ID
    # # Note this is done separately for upstream and dowsntream exons
    # # So as long as tx is represented once, it means a unique combination of flanking exons
    '''

    # Grouping by intron ID then exon coordinates
    rep_tx = (set(df.groupby([id_col + suffix, "Start", "End"])
                  .apply(lambda df: df[id_col].sort_values().iloc[0])
                  .values
                  )
              )

    # eprint(rep_tx)


    # Now subset to representative/selected transcripts
    sub_df = df[df[id_col].isin(rep_tx)]

    return sub_df


def get_straddled_exons(exons, introns, id_col="transcript_id", suffix="_b"):
    '''
    Return exons straddling (upstream and downstream) the input introns

    Params:
    exons - pr.PyRanges() -
    introns - pr.PyRanges() -
    id_col - str - name of column containing group identifier for intervals
    '''

    exons_cols = exons.columns.tolist()
    introns_cols = introns.columns.tolist()

    # Remove Chromosome - is never suffixed/duplicated upon join
    introns_cols = [col for col in introns_cols if col != "Chromosome"]

    # list of cols only in introns - these stay the same in a merge/join
    added_cols = [col + suffix if col in exons_cols else col for col in introns_cols]

    # # list of cols shared between dfs - need to add suffixes[1]
    # shared_cols = [col for col in to_merge_cols if col in exons_cols]
    #
    #
    # # Define cols names added from introns by join that will want to remove at the end
    # out_shared = [col + suffixes[1] for col in shared_cols]
    # target_cols = out_shared + only_cols


    # ids_introns = set(introns.as_df()[id_col])

    exons_olap = exons.join(introns,
                            how=None,  # Only keep overlapping exons
                            strandedness="same",
                            slack=1,  # bookended intervals can overlap
                            suffix=suffix
                            )

    # TODO: subset instead for exact matches at 5'end (upstream exon) or 3'end (downstream exon) of intron. Would then be agnostic to group identifiers
    exons_olap = exons_olap.subset(lambda df: _df_match_5p_adj(df, suffix) | _df_match_3p_adj(df, suffix))

    # Now check that given transcript has two exons present (i.e. upstream & downstream exon)
    n_pre_nexons = _n_ids(exons_olap, id_col)
    eprint(f"Number of transcripts with at least 1 exon exactly matching either side of intron - {n_pre_nexons}")

    # Now check that given transcript/ identifier ('id_col') has two exons present (i.e. upstream & downstream exon)
    ids_2exons = set(exons_olap.as_df().groupby(id_col)["Start"].nunique().loc[lambda x: x == 2].index)
    eprint(f"Number of transcripts with exactly bookended exons flanking the intron - {len(ids_2exons)}")

    exons_olap = exons_olap.subset(lambda df: df[id_col + suffix].isin(ids_2exons))

    # Collapse duplicate decoy exon pairs to a single represenative entry
    exons_olap = exons_olap.apply(lambda df: _collapse_straddled_exons(df, id_col, suffix))

    return exons_olap.drop(added_cols)


def _df_extend_intron(df, id_col, suffix):
    '''
    Extend intron to boundaries of bookended/straddled exons
    I
    '''

    # for each intron, subdivide into bookended intervals
    # find a set of transcript IDs to represent each individual interval
    # If duplicate tx for a given interval , sort alphabetically then take the first tx ID
    # Note this is done separately for upstream and dowsntream exons
    # So as long as tx is represented once, it means a unique combination of flanking exons
    rep_tx = (set(df.groupby([id_col, "Start" + suffix, "End" + suffix])
                  .apply(lambda df: df[id_col + suffix].sort_values().iloc[0])
                  .values
                  )
              )

    # eprint(rep_tx)

    # return rep_tx

    # Now subset to representative/selected transcripts
    sub_df = df[df[id_col + suffix].isin(rep_tx)]

    # Now extend introns for every selected transcript
    grouped = sub_df.groupby([id_col, id_col + suffix])

    # Select start and end coords

    # Starts should be smallest value (left-most)
    starts = grouped["Start" + suffix].min().reset_index().rename(columns={"Start_b": "Start"})
    # End coord should be largest value (right-most)
    ends = grouped["End" + suffix].max().reset_index().rename(columns={"End_b": "End"})

    # For all other cols, values should be identical so just pick the first that appears in each case
    other_cols = [col for col in df.columns if col not in ["Start", "End"]]

    others = grouped[other_cols].first().reset_index(drop=True)
    # eprint(others)

    # Combine updated coords with other metadata cols
    combined_coords = starts.merge(ends, on=[id_col, id_col + suffix])
    # eprint(combined_coords)

    combined = combined_coords.merge(others, on=[id_col, id_col + suffix])

    return combined[df.columns.tolist()]


def get_retention_event(introns, exons, id_col="transcript_id", suffix="_b"):
    '''
    Merge introns with 'straddled' exons (upstream and downstream exons) to represent a complete intron retention event
    '''

    exons_cols = exons.columns.tolist()
    introns_cols = introns.columns.tolist()

    # Remove Chromosome - is never suffixed/duplicated upon join
    # want to keep transcript-id from exons so know which ref transcript this decoy originated from
    exons_cols = [col for col in introns_cols if col not in ["Chromosome", "transcript_id"]]

    # list of cols only in exons - these stay the same in a merge/join
    added_cols = [col + suffix if col in introns_cols else col for col in exons_cols]


    introns_olap = introns.join(exons,
                            how=None,  # Only keep introns with overlap
                            strandedness="same",
                            slack=1,  # bookended intervals can overlap (this picks up )
                            suffix=suffix
                            )

    # eprint(introns_olap)

    # subset to where exact matches of intron 5'end to exon 3'end | intron 3' end to exon 5'end
    introns_olap = introns_olap.subset(lambda df: _df_match_5p_adj(df, suffix) | _df_match_3p_adj(df, suffix))

    # eprint(introns_olap)

    n_pre_nexons = _n_ids(introns_olap, id_col + suffix)
    eprint(f"Number of transcripts with at least 1 exon exactly matching either side of intron - {n_pre_nexons}")

    # Now check that given transcript/ identifier ('id_col') has two exons present (i.e. upstream & downstream exon)
    ids_2exons = set(introns_olap.as_df().groupby(id_col + suffix)["Start" + suffix].nunique().loc[lambda x: x == 2].index)
    eprint(f"Number of transcripts with exactly bookended exons flanking the intron - {len(ids_2exons)}")

    introns_olap = introns_olap.subset(lambda df: df[id_col + suffix].isin(ids_2exons))

    # Extend introns to boundaries of flanking exons
    introns_ext = introns_olap.apply(lambda df: _df_extend_intron(df, id_col, suffix))

    return introns_ext.drop(added_cols)



def add_3p_extension_length(gr,
                            ref_gr,
                            id_col="transcript_id",
                            out_col="3p_extension_length",
                            suffix="_b",
                            ref_cols_to_keep=["gene_id",
                                              "transcript_id",
                                              "exon_id",
                                              "Start",
                                              "End",
                                              "region_rank"],
                            nb_cpu=1):
    '''
    Add column '3p_extension_length' reporting 3'end extension of overlapping regions in gr & ref_gr (distance relative to gr)
    Note that for each unique ID (id_col) in gr, the smallest extension will be reported
    Avoids cases where overlaps with short & long isoforms, but tx is just reassembly of long isoform
    (so it's just reassembly of longer isoform, but is an extension relative to shorter ref isoform)
    '''

    # Find columns unique to ref_gr, so can drop at the end (no suffix added)
    not_cols = gr.columns.tolist()



    # Set of reference cols want to drop at end (either suffixed or not)
    # Note that Chromosome is never copied + suffixed in a join
    ref_cols = list(set(ref_gr.columns) - set(["Chromosome"]) - set(ref_cols_to_keep))

    ref_to_drop = [col if col not in not_cols else col + suffix for col in ref_cols]

    # Only overlapping intervals kept
    joined = gr.join(ref_gr,
                     strandedness="same",
                     how=None,
                     suffix=suffix,
                     nb_cpu=nb_cpu)

    joined = joined.assign(out_col,
                           lambda df: df["End"] - df["End" + suffix] if (df["Strand"] == "+").all() else
                           df["Start" + suffix] - df["Start"],
                           nb_cpu=nb_cpu)


    # To avoid capturing extensions of shorter isoforms (that really just are the longer known isoform)
    # Pick the smallest extension for each transcripts
    joined = joined.apply(lambda df: df.sort_values([id_col, out_col],
                                                    ascending=True).drop_duplicates(subset=[id_col],
                                                                                    keep="first"),
                          nb_cpu=nb_cpu)


    # eprint(joined.columns)

    return joined.drop(ref_to_drop)


def _df_5p_end_tolerance(df, rank_col, first_key, first_5p_tolerance, other_5p_tolerance, suffix):
    '''
    '''

    if (df["Strand"] == "+").all():
        # 5'end = Start
        # conditions = [(df[rank_col] == first_key) & (abs(df["Start"] - df["Start_b"]) <= first_5p_tolerance),
        #               (df[rank_col] != first_key) & (abs(df["Start"] - df["Start_b"]) <= other_5p_tolerance)
        #               ]
        # choices = [True, True]
        decisions = np.where((df[rank_col] == first_key) & (abs(df["Start"] - df["Start" + suffix]) <= first_5p_tolerance) |
                             (df[rank_col] != first_key) & (abs(df["Start"] - df["Start" + suffix]) <= other_5p_tolerance),
                             True,
                             False
                             )

    elif (df["Strand"] == "-").all():
        # 5'end = End
        decisions = np.where((df[rank_col] == first_key) & (abs(df["End"] - df["End" + suffix]) <= first_5p_tolerance) |
                             (df[rank_col] != first_key) & (abs(df["End"] - df["End" + suffix]) <= other_5p_tolerance),
                             True,
                             False
                             )


    return pd.Series(decisions, index=df.index)


def _df_add_event_type(df, id_col, rank_col, rkey2key, collapse_by_id=True):
    '''
    '''

    assert isinstance(rkey2key, dict)

    # eprint(df[rank_col].drop_duplicates())
    # eprint(rkey2key)

    assert any([True if key in set(df[rank_col]) else False for key in rkey2key.keys()])

    conditions = [df[rank_col] == key for key in rkey2key.keys()]

    choices = list(rkey2key.values())

    decisions = np.select(conditions, choices)

    if not collapse_by_id:
        return pd.Series(decisions, index=df.index)

    else:
        # Some IDs could have multiple assignments based on overlapping exon
        # Collapse to comma-separated string if multiple

        if not df.empty:
            to_join = df.loc[:, [id_col]]
            to_join["event_type_decision"] = pd.Series(decisions, index=df.index)

            to_join = (to_join.groupby(id_col)
                       ["event_type_decision"]
                       .agg(lambda x: ",".join(sorted(set(x)))) # sorted converts set --> list
                       .reset_index() #return tx_id to column
                       )

            # Now join by tx_id to get collapsed ID in correct (original) index
            df2 = df.merge(to_join, on=id_col)

            return df2["event_type_decision"]

        else:
            return pd.Series(decisions, index=df.index)



# def find_extension_events(novel_le,
#                           ref_exons,
#                           min_extension_length,
#                           first_5p_tolerance,
#                           other_5p_tolerance,
#                           tolerance_filter=True,
#                           return_filtered_ids=True,
#                           id_col="transcript_id",
#                           rank_col="region_rank",
#                           event_type_outcol="event_type",
#                           first_key="first_exon_extension",
#                           internal_key="internal_exon_extension",
#                           last_key="last_exon_extension",
#                           suffix="_ref"):
#     '''
#     Return gr with last exons that extend a reference exon
#     Criteria:
#     1 - Overlap with reference exon, novel last exon 3'end is downstream of ref_exon 3'end
#     2 - Smallest extension length is >= min_extension_length
#         - Smallest to avoid reassembly of longer isofrm being called as extension of shorter isoform
#     3 - 5'ends of novel_le & ref_exons match within given tolerance (nt)
#         - first exons should be more lenient (imprecision of TSS annotation, reassembly not nt precise)
#     '''

#     assert isinstance(min_extension_length, int)
#     assert isinstance(first_5p_tolerance, int)
#     assert isinstance(other_5p_tolerance, int)

#     # Find events extending 3'ends of ref exon, add 3p_extension_length (nt) column
#     # 1 length per le returned (smallest)
#     # Events with no overlap with ref exons are dropped
#     novel_le_ext = add_3p_extension_length(novel_le, ref_exons, suffix=suffix)

#     pre_len_ids = set(novel_le_ext.as_df()[id_col])

#     eprint(f"Number of events with any 3' ref exon extension - {len(pre_len_ids)}")

#     # eprint(f"SIN3B - {sin3b_tx in pre_len_ids}")

#     ext_len_dist = (novel_le_ext.as_df()
#                     ["3p_extension_length"]
#                     .describe(percentiles=[i * 0.1 for i in range(1,11,1)])
#                     )

#     eprint(f"3' extension length distribution\n{ext_len_dist}")

#     # Subset for extensions of min length
#     novel_le_ext = novel_le_ext.subset(lambda df: df["3p_extension_length"] >= min_extension_length)
#     post_len_ids = set(novel_le_ext.as_df()[id_col])

#     eprint(f"After minimum length filter - {min_extension_length} - number of extension events - {len(post_len_ids)}")
#     # eprint(f"SIN3B - {sin3b_tx in post_len_ids}")


#     if return_filtered_ids:
#         # Track IDs which do not pass min length filter
#         # In ds steps (e.g. finding spliced events, they could be considered 'spliced' events because no ref overlap, & ref last intron matches novel last intron)
#         min_len_filt_ids = pre_len_ids - post_len_ids
#         # eprint(f"SIN3B fails min length filter - {sin3b_tx in min_len_filt_ids}")

#     if tolerance_filter:
#         # Check 5'ends overlap within given tolerance
#         novel_le_ext = novel_le_ext.subset(lambda df: _df_5p_end_tolerance(df,
#                                                                            rank_col,
#                                                                            first_key,
#                                                                            first_5p_tolerance,
#                                                                            other_5p_tolerance,
#                                                                            suffix
#                                                                            )
#                                            )
#         post_tol_ids = set(novel_le_ext.as_df()[id_col])

#         eprint(f"After 5'end match tolerance filter, number of events - {len(post_tol_ids)}")

#         if return_filtered_ids:
#             end_5p_filt_ids = post_tol_ids - post_len_ids

#     else:
#         eprint("No 5'end match tolerance filtering performed...")
#         if return_filtered_ids:
#             end_5p_filt_ids = set()

#     # assign a 'event_type' column based on overlapping exon 'rank'
#     novel_le_ext = novel_le_ext.assign(event_type_outcol,
#                                        lambda df: _df_add_event_type(df,
#                                                                      id_col,
#                                                                      rank_col,
#                                                                      rkey2key={"first": first_key,
#                                                                                "internal": internal_key,
#                                                                                "last": last_key}
#                                                                      )
#                                        )

#     ev_types = novel_le_ext.as_df()[[id_col,
#                                       event_type_outcol]
#                                     ].drop_duplicates().value_counts(subset=[event_type_outcol])

#     eprint(f"Number of events of each type- {ev_types}")

#     if return_filtered_ids:
#         return novel_le_ext, post_len_ids, post_tol_ids

#     else:
#         return novel_le_ext


def main(in_gtf,
         in_le_gtf,
         output_prefix):
    '''
    '''

    eprint("reading in reference GTF...")
    gtf = pr.read_gtf(in_gtf)

    eprint("Extracting reference introns...")
    introns = gtf.features.introns(by="transcript")

    # Collapse duplicated introns to a single representative entry
    # If an intron is identical between transcripts, sort in alphabetical order and pick the first transcript to represent
    # eprint(introns)
    eprint("Remove duplicate introns, selecting representative transcript_id in alphabetical order in case of multiple values...")
    introns = introns.apply(lambda df: df.groupby(["gene_id", "Start", "End"]).apply(lambda df: df.sort_values(by="transcript_id").head(n=1)).reset_index(drop=True))
    # eprint(introns)
    # Need 'Feature' column in GTF to be 'exon' so gffread will extract FASTA sequence for it
    introns.Feature = "exon"

    eprint("Extracting and annotating reference exons as first, internal or last...")
    exons = gtf.subset(lambda df: df.Feature == "exon")
    exons = add_region_number(exons, feature_key="exon",out_col="exon_number")
    exons = add_region_rank(exons)

    # Get reference last exons
    ref_le = get_terminal_regions(exons)


    # Get all non-last reference exons
    ref_e_nl = pr.concat([get_terminal_regions(exons, which_region="first"),
                          get_internal_regions(exons)]
                         )

    # Remove sections of ref exons overlapping first/internal exons
    # This is so matches exactly with how did previously for PAPA
    ref_le = ref_le.subtract(ref_e_nl, strandedness="same")
    # ref_le.event_type_temp = "NULL"

    eprint("Reading in input last exons")
    le = pr.read_gtf(in_le_gtf)
    # le.event_type_temp = "NULL"

    # Store input transcript IDs. WIll allow to track which LEs correspond to input txipts
    le_tx_ids = set(le.transcript_id)

    if "Cluster" in le.columns:
        le = le.apply(lambda df: df.rename(columns={"Cluster": "Cluster_orig"}))

    if "le_number" in le.columns:
        le = le.apply(lambda df: df.rename(columns={"le_number": "le_number_orig"}))

    # Find the introns that overlap with input last exons
    # These will serve as 'dummy introns'
    # introns_le_in = introns.overlap(le, strandedness="same")
    # eprint(le.subset(lambda df: df["event_type"].str.contains("extension", regex=False)))

    # NaNs don't return boolean with str.contains below (?!)

    # Reassign event type - if annotated txs in input GTF these are not annotated with event type
    # Bleedthrough events would otherwise be excluded from analysis
    # sin3b_tx = le[le.gene_name == "SIN3B"].transcript_id.iloc[0]



    # le_bld, missed_len, missed_tol = find_extension_events(le, ref_e_nl, min_extension_length=100, first_5p_tolerance=100, other_5p_tolerance=0, return_filtered_ids=True, event_type_outcol="event_type_temp")
    # eprint(le_bld[["gene_id", "gene_name", "event_type", "event_type_temp"]])
    # eprint(le_bld[le_bld.gene_name == "SIN3B"])
    # eprint(missed_len.intersection({sin3b_tx}))
    # eprint(missed_tol.intersection({sin3b_tx}))
    # le = le.assign("event_type",
    #                lambda df: pd.Series(np.where(df["event_type"].isna(), "NA", df["event_type"])))

    # extract bleedthrough events from input last exons
    eprint("identifying introns overlapping bleedthrough events...")
    le_bld = le.subset(lambda df: df["event_type_simple"].str.contains("bleedthrough", regex=False))
   

    # find introns overlapping bleedthrough last exons
    introns_le_in_bld = introns.overlap(le_bld,
                                        strandedness="same")
    # eprint(introns_le_in_bld[introns_le_in_bld.gene_id == "SIN3B"])



    eprint("Expanding overlapping intron to straddled exons to represent intron retention event...")
    # Extend introns to the boundaries of straddled exons for bleedthroughs (reflect an IR event)
    introns_le_in_bld_ext = get_retention_event(introns_le_in_bld, ref_e_nl)

    introns_le_in_bld_ext = (introns_le_in_bld_ext.assign("transcript_id",
                                                          lambda df: df["transcript_id"] + "_" + df["transcript_id_b"] + "_decoy_ir").drop("transcript_id_b")
                                                 .assign("gene_id", lambda df: df["gene_id"] + "_decoy_ir")
                                                 )
    

    # Now extract introns overlapping other event types (to act as decoys)
    eprint("identifying introns overlapping non-bleedthrough events to act as decoys...")
    le_n_bld = le.subset(lambda df: ~df["event_type_simple"].str.contains("bleedthrough", regex=False))

    # find introns that overlap non-bleedthroughs
    introns_le_in_n_bld = introns.overlap(le_n_bld, strandedness="same")

    # assign decoy tx id for introns overlapping non-bleedthroughs
    introns_le_in_n_bld = (introns_le_in_n_bld.assign("transcript_id",
                                                      lambda df: df["transcript_id"] + "_" + "NA" + "_decoy_ir")
                            .assign("gene_id", lambda df: df["gene_id"] + "_decoy_ir")
                                                     )

    introns_le = pr.concat([introns_le_in_bld_ext, introns_le_in_n_bld])

    # # if wanted to do the same for reference events would do so here...
    # introns_le_ref = pr.PyRanges()  # introns.overlap(ref_le, strandedness="same")


    # Now define 'decoy exons' that represent a splicing event that skips the bleedthrough into the downstream (internal) exon
    # Extract the ref exons either side of the dummy intron
    # Bit like having the 'spliced in' region of a splice junction
    # These will serve as 'dummy exons' - not included in PPAU calcs but prevents all reads being assigned to last exon
    eprint("Identifying bleedthrough skipping splicing events to act as decoys...")
    exons_le_in = get_straddled_exons(exons, introns_le_in_bld)

    # assign dummy gene & transcript ID
    exons_le_in = (exons_le_in.assign("transcript_id",
                           lambda df: df["transcript_id"] + "_decoy_exon")
                .assign("gene_id",
                        lambda df: df["gene_id"] + "_decoy_exon"))
    



    # Annotate ref & input last exons to last exon IDs, can generate tx2le, le2gene tables

    # before joining, need a common gene_id column containing the reference gene ID (if a nvoel transcript, StringTie's gene_id attribute is a self-constructed one)
    # so that last exons (ref & novel) can be considered under same group (and le_ids be assigned appropriately)

    ref_le.gene_id_common = ref_le.gene_id
    le.gene_id_common = le.ref_gene_id

    le_comb = pr.concat([ref_le, le])
    le_comb = check_concat(le_comb)
    
    # assign a common event type col (annotating le_ids requires special consideration for distal 3'UTR extensions (consider as unique isoform), reference last exons don't contain this annotation so just need to assign a dummy value)
    le_comb = le_comb.assign("event_type_common",
                            lambda df: pd.Series(np.where(df["event_type_simple"].isna(),"NULL", df["event_type_simple"]))
                            )
     

    eprint("Annotating le_ids for downstream quantification...")
    le_comb = annotate_le_ids_ext(le_comb, id_col="gene_id_common", event_type_col="event_type_common", ext_key="distal_3utr_extension").drop("event_type_common")


    # To line up with PAPA GTF, need to remove last exons completely contained w/in annotated internal exons
    # (After assigning le_ids (why I've done this I do not know...))

    le_comb_c = le_comb.overlap(ref_e_nl, strandedness="same", how="containment", invert=True)

    # Some le_ids can be dropped if they are completely contained within non-last exons
    le_ids_dropped = set(le_comb_c.le_id) - set(le_comb.le_id)
    eprint(f"Number of last exon IDs dropped due to complete containment inside ref overlapping exons - {len(le_ids_dropped)}")


    # Now want a full GTF with last exons (to quant & summarise) & dummy exons/introns
    quant_combined = pr.concat([le_comb_c, introns_le, exons_le_in])


    eprint(f"Writing quantification-ready GTF to file - {output_prefix + '.quant.gtf'}")

    quant_combined.to_gtf(output_prefix + ".quant.gtf")

    eprint(f"Writing 'tx2gene' (transcript_id | gene_id) to TSV... - {output_prefix + '.tx2gene.tsv'}")

    (le_comb.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False)) # remove single isoform genes (keep='False' marks all duplicates as True (so keep these))
     .as_df()
     [["transcript_id", "gene_id"]]
     .drop_duplicates()
     .to_csv(output_prefix + ".tx2gene.tsv",
             sep="\t",
             index=False,
             header=True)
     )

    eprint(f"Writing 'tx2le' (transcript_id | le_id) to TSV... - {output_prefix + '.tx2le.tsv'}")

    (le_comb.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False))
     .as_df()
     [["transcript_id", "le_id"]]
     .drop_duplicates()
     .sort_values(by="le_id")
     .to_csv(output_prefix + ".tx2le.tsv",
             sep="\t",
             index=False,
             header=True)
     )

    eprint(f"Writing 'le2gene' (le_id | gene_id) to TSV... - {output_prefix + '.le2gene.tsv'}")
    (le_comb.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False))
     .as_df()
     [["le_id", "gene_id"]]
     .drop_duplicates()
     .sort_values(by="le_id")
     .to_csv(output_prefix + ".le2gene.tsv",
             sep="\t",
             index=False,
             header=True)
     )

    eprint(f"Writing 'le2name' (le_id | gene_name) to TSV... - {output_prefix + '.le2name.tsv'}")
    (le_comb.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False))
     .as_df()
     [["le_id", "gene_name"]]
     .drop_duplicates()
     .dropna()
     .sort_values(by="le_id")
     .to_csv(output_prefix + ".le2name.tsv",
             sep="\t",
             index=False,
             header=True)
     )

    eprint(f"Writing 'input tx2le' (transcript-id | le_id)  to TSV... - {output_prefix + '.input_tx2le.tsv'}")

    (le_comb.subset(lambda df: df.duplicated(subset=["gene_id"], keep=False))
     .subset(lambda df: df["transcript_id"].isin(le_tx_ids))
     .as_df()
     [["transcript_id", "le_id"]]
     .drop_duplicates()
     .dropna()
     .sort_values(by="le_id")
     .to_csv(output_prefix + ".input_tx2le.tsv",
             sep="\t",
             index=False,
             header=True)
     )


if __name__ == '__main__':

    # usage = "python add_dummies_to_gtf.py REF_GTF LAST_EXONS_GTF OUTPUT_PREFIX"
    #
    # if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
    #     eprint(usage)
    #     sys.exit(0)


    descrpn = """Generate quantification ready GTF of last exons (with decoy transcripts corresponding to overlapping exons/introns) & group transcripts according to shared last exon"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-i",
                        "--input-gtf",
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to input GTF of last exons")

    parser.add_argument("-r",
                        "--ref-gtf",
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to reference GTF file containing gene/transcript annotations")

    parser.add_argument("-d",
                        "--decoy-type",
                        required=True,
                        default="full",
                        choices=["full", "unique"],
                        help="Type of decoy transcripts/introns to generate. 'full' keeps the full bleedthrough last exon, adding the overlapping internal exon and its downstream exon as a decoy transcript along with the overlapping intron. 'unique' takes the specific/unique region of the bleedthrough last exon (not overlapping with ref internal exons) and adds the full overlapping intron as a decoy transcript."
                        )

    parser.add_argument("-o",
                        "--output-prefix",
                        default="last_exons_w_decoys",
                        help="Output prefix for GTF, tx2le,tx2gene tables etc.")


    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.ref_gtf, args.input_gtf, args.output_prefix)
