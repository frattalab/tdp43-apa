#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
from typing import Literal, Optional
import sys

'''
Extract splice junctions for AS-ALE (/ last exon spliced) events.

Due to a pretty silly move on my part, my pipeline does not output splice junctions for novel spliced last exons (nor last exons)
Given want to use SJs to quantify these events in NYGC data, need to extract using only last exon coordinates.

In pipeline, approach to define last exons:
- 3'end doesn't overlap any annotated exon
- terminal splice junction must match an annotated splice junction at the 5'end. 

Can exploit these criteria to grab SJs 'in reverse'

Strategy for novel last exons:
- Extract introns (SJs) from annotated transcripts
- find introns in which last exons are completely contained
- (w/ some form of pr.join), construct a SJ with start coordinate of annotated intron and end coordinate = start of exon coordinate
- As a sanity check, compare with AL's provided list of cryptic SJs inferred with MAJIQ - for common genes, do we get the same SJ coordinates

Strategy for annotated last exons:
- Any spliced last exon that is not captured with novel LE strategy
- pr.join with slack=1 (allow bookended intervals to overlap) annotated SJs with last exon coordinates
- Subset for exact matches at SJ 3'end/last exon 5'end
- Retain the SJs for downstream analysis


### TODOs
- Switch to clean script for generating SJs
- Create a TXT file with masked IDs (4 for now)
- For last exons not within gene bodies ('gene' entry), select the terminal introns for these transcripts, then do .nearest() to attach to last exon. Then can reconstruct SJ with existing function. 
'''


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
        print("df is empty - returning empty pd.Series")
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
        print("pr.assign returned KeyError. Converting {} to int via pandas df conversion".format(region_number_col))

        mod_gr = gr.as_df()
        mod_gr[region_number_col] = mod_gr[region_number_col].astype(float).astype(int)
        mod_gr = pr.PyRanges(mod_gr)


    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    mod_gr = mod_gr.apply(lambda df: df.sort_values(by=[id_col, region_number_col], ascending=True),
                          nb_cpu=1)


    # Filter out single-exon transcripts
    if filter_single:
        print("Filtering for multi-exon transcripts...")
        print("Before: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))

        # Setting to 'False' marks all duplicates as True (so keep these)
        mod_gr = mod_gr.subset(lambda df: df.duplicated(subset=[id_col], keep=False), nb_cpu=1)

        print("After: {}".format(len(set(mod_gr.as_df()[id_col].tolist()))))


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




def _df_intron_exon_to_sj(df: pd.DataFrame) -> pd.DataFrame:
    '''Convert joined intron + overlapping exon df of PyRanges object to splice junction coordinates (the first bin between the two intervals)

    If + strand
    Start coord = intron start (i.e. 'Start' column), End coord = exon start (i.e. 'Start_b' column)
    If - strand
    Start coord = Exon start (i.e. 'End_b' column), End coord = intron start (i.e. 'End' column)
    
    Parameters
    ----------
    df : pd.DataFrame (internal to pr.PyRanges object)
        _description_

    Returns
    -------
    pd.DataFrame
        _description_

    Raises
    ------
    Exception
        Strand column of df doesn't contain all '+' or '-'
    '''

    df[["Start_intron", "End_intron"]] = df[["Start", "End"]]

    if (df.Strand == "+").all():
        df["End"] = df["Start_b"]

    elif (df.Strand == "-").all():
        df["Start"] = df["End_b"]

    else:
        raise Exception("Strand column must only contain all '+' or '-'")
    
    return df






def intron_le_to_sj(gr: pr.PyRanges,
                    introns_gr: pr.PyRanges) -> pr.PyRanges:
    '''_summary_

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    introns_gr : pr.PyRanges
        _description_

    Returns
    -------
    pr.PyRanges
        _description_
    '''

    # overlap-join introns with LEs
    ref_introns_le = introns_gr.join(gr, strandedness="same")

    # parse SJ coordinates from overlapping exon and intron (intron start - last exon start)
    # Note that reference introns can be duplicated across transcript IDs, so only retain one row per intron (don't care about overlapping tx ID)
    le_sj = ref_introns_le.apply(_df_intron_exon_to_sj).drop_duplicate_positions()

    return le_sj


def annot_le_to_sj(gr: pr.PyRanges,
                   introns_gr: pr.PyRanges) -> pr.PyRanges:
    '''_summary_

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    introns_gr : pr.PyRanges
        _description_

    Returns
    -------
    pr.PyRanges
        _description_
    '''

    # now to get SJs for annotated ALEs
    # overlap-join SJs and LEs, allowing bookended/directly adjacent to overlap
    # subset for direct match at intron 3'end, last exon 5'End
    # Note that reference introns can be duplicated across transcript IDs, so only retain one row per intron (don't care about tx ID)
    le_sj_ref = (introns_gr.join(gr, strandedness="same", slack=1)
                 .subset(lambda df: ((df.Strand == "+") & (df["End"] == df["Start_b"])) |
                         ((df.Strand == "-") & (df["Start"] == df["End_b"]))
                         )
                         .drop_duplicate_positions()
           )
    
    # start and end coordinates are already the SJ, so no need to extract

    return le_sj_ref


def distal_le_to_sj(gr: pr.PyRanges,
                    introns_gr: pr.PyRanges,
                    id_col: Optional[str] = "transcript_id",
                    region_num_col: Optional[str] = "intron_number") -> pr.PyRanges:
    '''_summary_

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    introns_gr : pr.PyRanges
        _description_
    id_col : Optional[str], optional
        _description_, by default "transcript_id"
    region_num_col : Optional[str], optional
        _description_, by default "intron_number"

    Returns
    -------
    pr.PyRanges
        _description_
    '''
    
    # extract terminal introns for each transcript/group of intervals
    last_introns = get_terminal_regions(introns_gr, id_col=id_col, region_number_col=region_num_col)

    # find .nearest() fownstream last exon for each terminal intron 
    ref_introns_le = gr.nearest(last_introns, strandedness="same", overlap=False, how="upstream", suffix="_int")

    # parse SJ coordinates from overlapping exon and intron (intron start - last exon start)
    # Note that reference introns can be duplicated across transcript IDs, so only retain one row per intron (don't care about overlapping tx ID)
    # function expects start cols to be intron coords, _b coords the exon coords
    ref_introns_le = ref_introns_le.apply(lambda df: df.rename(columns={"Start": "Start_b",
                                                                        "End": "End_b",
                                                                        "Start_int": "Start",
                                                                        "End_int": "End"})
                                                                        )    
    
    le_sj = ref_introns_le.apply(_df_intron_exon_to_sj).drop_duplicate_positions()

    return le_sj


def le_to_sj(gr: pr.PyRanges,
             introns_gr: pr.PyRanges,
             le_type: Literal['distal', 'contained', 'annotated'],
             id_col: Optional[str] = "transcript_id",
             region_num_col: Optional[str] = "intron_number") -> tuple[pr.PyRanges, set]:
    '''_summary_

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    introns_gr : pr.PyRanges
        _description_
    le_type : Literal[&#39;distal&#39;, &#39;contained&#39;, &#39;annotated&#39;]
        _description_

    Returns
    -------
    tuple[pr.PyRanges, set]
        _description_
    '''
    assert le_type in ['distal', 'contained', 'annotated']

    if le_type == "distal":
        # check introns_gr has strand-aware, transcript-wise intron number (to extact terminal introns)
        assert region_num_col in introns_gr.columns
        assert id_col in introns_gr.columns
    

    le_ids = set(gr.Name)

    # Determine SJs based on type of last exon
    if le_type == "annotated":
        le_sj = annot_le_to_sj(gr, introns_gr)

    elif le_type == "contained":
        le_sj = intron_le_to_sj(gr, introns_gr)

    elif le_type == "distal":
        le_sj = annot_le_to_sj(gr, introns_gr, id_col, region_num_col)

    # finally generate a cleaned BED-ready file of novel, annotated intron contained SJs
    le_sj_clean = le_sj[["Name"]]
    le_sj_clean.Score = "."
    
    # double check if any IDs/regions have been dropped
    dropped_le_ids = le_ids.difference(set(le_sj_clean.Name))

    return le_sj_clean, dropped_le_ids



def main(le_bed_path,
         ref_gtf_path,
         output_prefix):
    
    ref_gtf = pr.read_gtf(ref_gtf_path)
    le = pr.read_bed(le_bed_path)

    # ref gtf has lots of attributes don't need, subset to minimal cols for analysis
    ref_gtf = ref_gtf[["Feature", "gene_id", "transcript_id", "gene_name"]]

    # get IDs from BED file of cryptic LEs
    le_ids = set(le.Name)

    # Identify different groups of last exons
    # 1. Distal to known gene coordinates (i.e. doesn't overlap known gene)
    # 2. Overlaps known gene, completely contained within annotated intron
    # 3. Overlaps known gene, is annotated last exon

    # Extract gene entries from GTF
    ref_genes = ref_gtf.subset(lambda df: df.Feature == "gene")

    # extract introns (SJs) for each transcript
    ref_introns = ref_gtf.features.introns(by="transcript")

    # assign a strand-aware order of introns within transcripts (1 = first, n = last)
    ref_introns = add_region_number(ref_introns)

    # store a df noting which transcript IDs share the same intron coordinates (useful if collapse later and want to track/retain)
    # intron2tx = ref_introns.as_df()[["Chromosome", "Start", "End", "Strand", "transcript_id"]]

    # get IDs that don't overlap with known genes
    le_distal = le.overlap(ref_genes, strandedness="same", invert=True)
    le_distal_ids = set(le_distal.Name)
    print(f"Putative number of LEs distal/downstream to annotated gene - {len(le_distal_ids)}")

    le_n_distal = le.subset(lambda df: ~df.Name.isin(le_distal_ids))
    
    # get LEs completely contained within introns
    le_cont = le_n_distal.overlap(ref_introns, strandedness="same", how="containment")
    le_cont_ids = set(le_cont.Name)

    print(f"Putative number of novel ALEs completely contained within annotated introns - {len(le_cont_ids)}")

    # remaining IDs are annotated
    le_annot_ids = le_ids.difference(le_distal_ids.union(le_cont_ids))
    le_annot = le_n_distal.subset(lambda df: df.Name.isin(le_annot_ids))
    print(f"Putative number of annotated ALEs- {len(le_annot_ids)}")

    # missing IDs due to categorisation
    le_ids_cat = le_annot_ids.union(le_cont_ids).union(le_distal)
    le_ids_missing = le_ids.difference(le_ids_cat)
    print(f"Number of missing IDs due to not fitting any category - {len(le_ids_missing)}")

    # generate SJs
    le_distal_sj, le_distal_ids_missing = le_to_sj(le_distal, ref_introns, le_type="distal")
    le_cont_sj, le_cont_ids_missing = le_to_sj(le_cont, ref_introns, le_type="contained")
    le_annot_sj, le_annot_ids_missing = le_to_sj(le_cont, ref_introns, le_type="annotated")

    print(f"Number of gene-distal IDs that lack a SJ - {len(le_distal_ids_missing)}")
    print(f"Number of intron-contained IDs that lack a SJ - {len(le_cont_ids_missing)}")
    print(f"Number of annotated IDs that lack a SJ - {len(le_annot_ids_missing)}")

    # generate df of missing ids and their reason
    missing_dict = {"not_categorised": le_ids_missing,
                    "gene_distal": le_distal_ids_missing,
                    "intron_contained": le_cont_ids_missing,
                    "annotated": le_annot_ids_missing}
    
    # TODO: must be a more straightforward way to do this...
    missing_df = (pd.DataFrame.from_dict(missing_dict, orient="index") # index = keys, values = columns
                  .melt(value_name="region_name", ignore_index=False) # pivot to longer format
                  .dropna()
                  .drop(columns="variable")
                  .reset_index()
                  .rename(columns={"index": "le_category"})
                  )

    # combine into single bed
    sj_bed = pr.concat([le_distal_sj, le_cont_sj, le_annot_sj]).sort()

    # output to file
    print("Writing output files...")
    sj_bed.to_bed(output_prefix + ".le.unmodified.junctions.bed")
    sj_bed.assign("End", lambda df: df.End.add(1)).to_bed(output_prefix + ".le.end_shift.junctions.bed")
    missing_df.to_csv(output_prefix + ".le.missing.tsv", sep="\t", index=False, header=True)


if __name__ == '__main__':

    descrpn = """usage: python last_exons_to_sj.py LAST_EXON_BED REFERENCE_GTF OUTPUT_PREFIX [-h/--help]

Positional arguments:
LAST_EXON_BED - BED file containing last exon coordinates produced by PAPA for which to infer splice junctions
REFERENCE_GTF - path to GTF file with reference/annotated transcripts
OUTPUT_PREFIX - prefix for output files. BED file of splice junctions suffixed with '.le.unmodified.junctions.bed', '.le.end_shift.junctions.bed' for SJs to use with bedops_parse_star_junctions pipeline (see details). TSV file suffixed with '.le.missing.tsv' storing input intervals for which splice junctions could not be inferred/extracted.
-h/--help - print usage message
    
Details:

what is 'end_shift.junctions.bed'? STAR's junction.tab fiels represents SJs in 1-based manner, and not necessarily representing only the intron coordinates (according to AL). To make SJs compatible with the pipeline I have to add 1 to the end coordinate. 

"""

    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
        print(descrpn)
        sys.exit()

    main(sys.argv[1], sys.argv[2], sys.argv[3])




   


