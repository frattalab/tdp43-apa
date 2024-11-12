#!/usr/bin/env python3

from __future__ import print_function
import pyranges as pr
import pandas as pd
import numpy as np
from helpers import add_region_number, get_internal_regions, get_terminal_regions, add_region_rank
import sys
import argparse
from typing import Iterable, Set, List


'''
Generate decoy transcript models for IPA (bleedthrough) and ALE ('spliced') last exons
- Bleedthroughs
    - spliced event (combo of two flanking internal exons)
    - intron retention (IR) - flanking exons + full intron
- ALEs
    - intron (full intron)

Input:
- Reference GTF
- GTF of last exons of interest split by event type ()

Output:
- GTF combining decoy models with last exons of interest
(optionally extending bleedthrough models to boundary of upstream exon, if provided PAPA quant GTF)
'''


def _n_ids(gr, id_col, is_df=False):

    assert id_col in gr.columns

    if not is_df:
        return len(set(gr.as_df()[id_col]))
    else:
        return gr[id_col].nunique()


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

    # print(rep_tx)


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
    print(f"Number of transcripts with at least 1 exon exactly matching either side of intron - {n_pre_nexons}")

    # Now check that given transcript/ identifier ('id_col') has two exons present (i.e. upstream & downstream exon)
    ids_2exons = set(exons_olap.as_df().groupby(id_col)["Start"].nunique().loc[lambda x: x == 2].index)
    print(f"Number of transcripts with exactly bookended exons flanking the intron - {len(ids_2exons)}")

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

    # print(rep_tx)

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
    # print(others)

    # Combine updated coords with other metadata cols
    combined_coords = starts.merge(ends, on=[id_col, id_col + suffix])
    # print(combined_coords)

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

    # print(introns_olap)

    # subset to where exact matches of intron 5'end to exon 3'end | intron 3' end to exon 5'end
    introns_olap = introns_olap.subset(lambda df: _df_match_5p_adj(df, suffix) | _df_match_3p_adj(df, suffix))

    # print(introns_olap)

    n_pre_nexons = _n_ids(introns_olap, id_col + suffix)
    print(f"Number of transcripts with at least 1 exon exactly matching either side of intron - {n_pre_nexons}")

    # Now check that given transcript/ identifier ('id_col') has two exons present (i.e. upstream & downstream exon)
    ids_2exons = set(introns_olap.as_df().groupby(id_col + suffix)["Start" + suffix].nunique().loc[lambda x: x == 2].index)
    print(f"Number of transcripts with exactly bookended exons flanking the intron - {len(ids_2exons)}")

    introns_olap = introns_olap.subset(lambda df: df[id_col + suffix].isin(ids_2exons))

    # Extend introns to boundaries of flanking exons
    introns_ext = introns_olap.apply(lambda df: _df_extend_intron(df, id_col, suffix))

    return introns_ext.drop(added_cols)


def expand_ids(ids: Iterable[str], delimiter: str = ',') -> Set[str]:
    '''Expands an iterable of (potentially) delimited strings into a set of individual strings.

    Parameters
    ----------
    ids : Iterable[str]
        An iterable of strings, where each string may contain 
                             multiple items separated by the specified delimiter.
    delimiter : str, optional
        delimiter (str, optional): The delimiter used to split the strings, by default ','

    Returns
    -------
    Set[str]
        _description_
    '''
    expanded_ids = set()
    for id in ids:
        if delimiter in id:
            expanded_ids.update(id.split(delimiter))
        else:
            expanded_ids.add(id)
    return expanded_ids


def decoy_assignment_category(df: pd.DataFrame, id_col: str = "le_id", indicator_col: str = "le_id_in", out_col: str = "decoy_assignment", out_keys: List[str] = ["all", "some", "none"]) -> pd.DataFrame:
    '''_summary_

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    id_col : str, optional
        _description_, by default "le_id"
    indicator_col : str, optional
        _description_, by default "le_id_in"
    out_col : str, optional
        _description_, by default "decoy_assignment"
    out_keys : List[str], optional
        _description_, by default ["all", "some", "none"]

    Returns
    -------
    pd.DataFrame
        _description_
    '''
    
    assert len(out_keys) == 3, f"out_keys must contain 3 values, n = {len(out_keys)}"
    
    # Group by group_id and calculate the mean of le_id_in
    summ = df.groupby(id_col)[indicator_col].mean().reset_index()
    
    # Determine the decoy_assignment (category) based on the mean
    summ[out_col] = np.select(
        [summ[indicator_col] == 1, summ[indicator_col] == 0],
        [out_keys[0], out_keys[2]],
        default=out_keys[1]
    )
    
    # Return the group_id and decoy_assignment columns
    return summ[[id_col, out_col]]


def main(ref_gtf: str,
         ipa_gtf: str,
         ale_gtf: str,
         output_prefix: str,
         gene_id_col: str = "ref_gene_id",
         event_id_col: str = "le_id",
         event_tx_col: str = "transcript_id",
         ref_gene_id_col: str = "gene_id"
         ):
    '''Generate decoy transcript models representing alternative processing decisions for ALE and IPA events

    Parameters
    ----------
    ref_gtf : str
        Reference GTF
    ipa_gtf : str
        GTF containing bleedthrough/IPA events for which to generate decoy transcripts
    ale_gtf : str
        GTF containing spliced/ALE events for which to generate decoy transcripts
    output_prefix : str
        Prefix for output GTFs, suffixed with ".decoys.gtf" (decoys only) or ".combined.gtf" (decoys + IPA + ALE)
    gene_id_col : str
        name of column containing gene ID in ipa_gtf/ale_gtf, by default "ref_gene_id"
    event_id_col : str, optional
        _name of column containing event identifier (i.e. le_id), by default "le_id"
    event_tx_col : str, optional
        Name of column containing identifier one below --event_id_col (i.e. id column for overlapping intervals that are grouped into a single event), by default "transcript_id"
    ref_gene_id_col : str
        name of column containing gene ID in ref_gtf, by default "gene_id"
        '''
    
    ipa = pr.read_gtf(ipa_gtf)
    ale = pr.read_gtf(ale_gtf)

    # Extract gene IDs from input target events (only process GTF for these genes)
    ipa_gene_ids = expand_ids(ipa.as_df()[gene_id_col])
    ale_gene_ids = expand_ids(ale.as_df()[gene_id_col])
    comb_gene_ids = ipa_gene_ids.union(ale_gene_ids)
 
    # track event IDs for downstream checks
    ipa_ids = set(ipa.as_df()[event_id_col])
    ale_ids = set(ale.as_df()[event_id_col])
    
    # track 'sub' event IDs mapped to each event (used to track/report decoy assignment status)
    ipa_ids2tx = ipa.as_df()[[event_id_col, event_tx_col]].drop_duplicates()
    ale_ids2tx = ale.as_df()[[event_id_col, event_tx_col]].drop_duplicates()
   
    print("Reading in reference GTF & filtering for input genes...")
    gtf = pr.read_gtf(ref_gtf).subset(lambda df: df[ref_gene_id_col].isin(comb_gene_ids))

    print("Extracting reference introns...")
    introns = gtf.features.introns(by="transcript")

    # Collapse duplicated introns to a single representative entry
    # If an intron is identical between transcripts, sort in alphabetical order and pick the first transcript to represent
    # (actual tx doesn't matter)
    # print(introns)
    print("Remove duplicate introns, selecting representative transcript_id in alphabetical order in case of multiple values...")
    # introns = introns.apply(lambda df: df.groupby([ref_gene_id_col, "Start", "End"]).apply(lambda df: df.sort_values(by="transcript_id").head(n=1)).reset_index(drop=True)).sort()
    # print(introns)
    # Need 'Feature' column in GTF to be 'exon' so gffread will extract FASTA sequence for it
    introns.Feature = "exon"

    print("Extracting and annotating reference exons as first, internal or last...")
    exons = gtf.subset(lambda df: df.Feature == "exon")
    # Collapse duplicated exons to a single representative entry (as for introns)
    # exons = exons.apply(lambda df: df.groupby([ref_gene_id_col, "Start", "End"]).apply(lambda df: df.sort_values(by="transcript_id").head(n=1)).reset_index(drop=True)).sort()
    exons = add_region_number(exons, feature_key="exon",out_col="exon_number")
    exons = add_region_rank(exons)

    # Get all non-last reference exons
    # ref_e_nl = pr.concat([get_terminal_regions(exons, which_region="first"),
    #                       get_internal_regions(exons)]
    #                      )
    
    print("Generating decoys for IPA/bleedthrough events...")
    # find introns overlapping bleedthrough last exons
    introns_ipa = introns.overlap(ipa, strandedness="same")
    print("Expanding overlapping intron to straddled exons (IPA/bleedthrough intron retention decoy)...")
    # Extend introns to the boundaries of straddled exons for bleedthroughs (reflect an IR event)
    introns_ipa_ext = get_retention_event(introns_ipa, exons)

    # assign dummy IDs (gene, tx) for decoy events
    introns_ipa_ext = (introns_ipa_ext.assign("transcript_id",
                                              lambda df: df["transcript_id"] + "_" + df["transcript_id_b"] + "_decoy_ir")
                                              .drop("transcript_id_b")
                                              .assign("gene_id", lambda df: df["gene_id"] + "_decoy_ir")
                                                 )
    
    # Now define 'decoy exons' that represent a splicing event that skips the bleedthrough into the downstream (internal) exon
    # Extract the ref exons either side of the dummy intron
    # Bit like having the 'spliced in' region of a splice junction
    # These will serve as 'dummy exons' - not included in PPAU calcs but prevents all reads being assigned to last exon
    print("Identifying flanking exons for intron containing bleedthrough/IPA (IPA/bleedthrough splicing decoy)...")
    exons_ipa = get_straddled_exons(exons, introns_ipa)
    # assign dummy gene & transcript ID
    exons_ipa = (exons_ipa.assign("transcript_id",
                           lambda df: df["transcript_id"] + "_decoy_exon")
                .assign("gene_id",
                        lambda df: df["gene_id"] + "_decoy_exon"))
    
    print("Identify introns containing ALE events (ALE intron retention decoy)")
    introns_ale = introns.overlap(ale, strandedness="same")
    
    # get ALEs not overlapping annotated introns
    ale_nintrons = ale.overlap(introns, strandedness="same", invert=True)
    
    
    nintrons_ale_ids = set(ale_nintrons.as_df()[event_id_col])
    print(f"Number of ALEs not overlapping annotated introns - {len(nintrons_ale_ids)}")
    
    

    # assign decoy tx id for introns overlapping non-bleedthroughs
    introns_ale = (introns_ale.assign("transcript_id",
                                      lambda df: df["transcript_id"] + "_" + "NA" + "_decoy_ir")
                                      .assign("gene_id", lambda df: df["gene_id"] + "_decoy_ir")
                                      )

    print("Sanity checking for complete assignment of decoy events...")
    print("Bleedthrough/IPA events")
    # exons_IPA - every IPA event should have an overlapping (or bookended) exon
    ipa_ex_ids = set(ipa.join(exons_ipa, strandedness="same", how=None, slack=1).as_df()[event_id_col])
    missing_ipa_ex_ids = ipa_ids.difference(ipa_ex_ids)
    if len(missing_ipa_ex_ids) > 0:
        print(f"IPA events missing flanking exon splicing decoy event, n = {len(missing_ipa_ex_ids)}")
        if len(missing_ipa_ex_ids) <= 25:
            print(sorted(list(missing_ipa_ex_ids)))
    
    # Every IPA event should also have an overlapping intron
    ipa_in_ids = set(ipa.overlap(introns_ipa_ext, strandedness="same").as_df()[event_id_col])
    missing_ipa_in_ids = ipa_ids.difference(ipa_in_ids)
    if len(missing_ipa_in_ids) > 0:
        print(f"IPA events missing overlapping intron retention decoy event, n = {len(missing_ipa_in_ids)}")
        if len(missing_ipa_in_ids) <= 25:
            print(sorted(list(missing_ipa_in_ids)))

    print("Spliced/ALE events")
    # All ALEs should also have an overlapping intron
    ale_in_ids2tx = ale.overlap(introns_ale, strandedness="same").as_df()[[event_id_col, event_tx_col]].drop_duplicates()
    ale_in_ids = set(ale_in_ids2tx[event_id_col])
    missing_ale_in_ids = ale_ids.difference(ale_in_ids)
    
    # decoy assignment category - check whether all txs of a le_id are overlapping intron
    ale_in_ids2tx = ale_ids2tx.merge(ale_in_ids2tx, how="left", on=event_tx_col, suffixes=[None, "_in"])
    
    # Make join column an indicator columm (1 if decoy assigned, 0 if not)
    joined_col = event_id_col + "_in"
    ale_in_ids2tx[joined_col] = np.where(ale_in_ids2tx[joined_col].isna(), 0, 1)
    # Categorise le_ids according to whether all/some/none of its transcripts/sub-ids are assigned a decoy transcript
    ale_decoys_summ = decoy_assignment_category(ale_in_ids2tx, id_col=event_id_col, indicator_col=joined_col)
    
    print("Summary counts for ALEs and assignment of decoy transcripts")
    print(ale_decoys_summ["decoy_assignment"].value_counts())
    
    # sanity check - ALEs not overlapping introns should be excluded from this analysis
    ale_in_nintron_ids = ale_in_ids.intersection(nintrons_ale_ids)
    missing_ale_in_nintron_ids = missing_ale_in_ids.intersection(nintrons_ale_ids)
    if len(missing_ale_in_ids) > 0:
        print(f"ALE events missing overlapping intron retention decoy event, n = {len(missing_ale_in_ids)}")
        if len(missing_ale_in_ids) <= 25:
            print(missing_ale_in_ids) 

        print(f"Intersection between IDs found in output and those allegedly with intron decoys, n = {len(ale_in_nintron_ids)}")
        if len(ale_in_nintron_ids) <= 25:
            print(ale_in_nintron_ids) 
        print(f"Intersection between IDs dropped and those without an overlapping intron, n = {len(missing_ale_in_nintron_ids)}")
        if len(missing_ale_in_nintron_ids) <= 25:
            print(missing_ale_in_nintron_ids) 
    
    print("Generating and outputting combined files...")
    # combine decoy events into a single GTF
    decoys_comb = pr.concat([exons_ipa, introns_ipa_ext, introns_ale]).sort()

    # combine decoys with input events
    all_comb = pr.concat([ale, ipa, decoys_comb]).sort()

    # Output to file
    decoys_comb.to_gtf(output_prefix + ".decoys.gtf")
    all_comb.to_gtf(output_prefix + ".combined.gtf")
    ale_decoys_summ.sort_values(by="decoy_assignment").to_csv(output_prefix + ".decoy_assignment.summary.ale.tsv", sep="\t", index=False, header=True)
    ale_in_ids2tx.rename(columns={"le_id_in": "decoy"}).to_csv(output_prefix + ".decoy_assignment.all.ale.tsv", sep="\t", index=False, header=True)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate decoy transcript models representing alternative processing decisions for ALE and IPA events"
    )

    parser.add_argument(
        "-r", "--ref-gtf", type=str, required=True,
        help="Reference GTF"
    )

    parser.add_argument(
        "-i", "--ipa-gtf", type=str, required=True,
        help="GTF containing bleedthrough/IPA events for which to generate decoy transcripts"
    )

    parser.add_argument(
        "-a", "--ale-gtf", type=str, required=True,
        help="GTF containing spliced/ALE events for which to generate decoy transcripts"
    )

    parser.add_argument(
        "-o", "--output-prefix", type=str, required=True,
        help='Prefix for output GTFs, suffixed with ".decoys.gtf" (decoys only) or ".combined.gtf" (decoys + IPA + ALE)'
    )

    parser.add_argument(
        "--gene_id_col", type=str, default="ref_gene_id",
        help='Name of column containing gene ID in ipa_gtf/ale_gtf, by default "ref_gene_id"'
    )

    parser.add_argument(
        "--event_id_col", type=str, default="le_id",
        help='Name of column containing event identifier (i.e. le_id), by default "le_id"'
    )
    
    parser.add_argument(
        "--event_tx_col", type=str, default="transcript_id",
        help='Name of column containing identifier one below --event_id_col (i.e. id column for overlapping intervals that are grouped into a single event), by default "transcript_id"'
    )

    parser.add_argument(
        "--ref_gene_id_col", type=str, default="gene_id",
        help='Name of column containing gene ID in ref_gtf, by default "gene_id"'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(
        ref_gtf=args.ref_gtf,
        ipa_gtf=args.ipa_gtf,
        ale_gtf=args.ale_gtf,
        output_prefix=args.output_prefix,
        gene_id_col=args.gene_id_col,
        event_id_col=args.event_id_col,
        event_tx_col=args.event_tx_col,
        ref_gene_id_col=args.ref_gene_id_col
    )