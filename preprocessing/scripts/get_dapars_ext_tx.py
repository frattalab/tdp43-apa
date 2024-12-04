#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import argparse
import sys
from typing import Literal
import logging


def _df_swap_coord(df: pd.DataFrame, change: Literal["Start", "End"], replace_col: str):
    """Swap one end of a coordinate range (i.e. either start or end) with a joined value

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
    """ """"""

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
    """_summary_

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    """

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


def get_terminal_regions(
    gr: pr.PyRanges,
    feature_col="Feature",
    feature_key="exon",
    id_col="transcript_id",
    region_number_col="exon_number",
    number_type: Literal["stranded", "unstranded"] = "stranded",
    which_region: str = "last",
    filter_single=False,
):
    """Return the first/last interval in group of intervals

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
    number_type : str, optional
        _description_, by default ["stranded", "unstranded"]
    which_region : str, optional
        _description_, by default "last"
    filter_single : bool, optional
        _description_, by default False

    Returns
    -------
    _type_
        _description_
    """

    assert number_type in ["stranded", "unstranded"]
    assert which_region in ["first", "last"]
    assert region_number_col in gr.columns.tolist()
    assert feature_col in gr.columns.tolist()
    assert id_col in gr.columns.tolist()

    # Make sure only 'exon' features are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [
        feature_key
    ], "only {} entries should be present in gr".format(feature_key)

    # Make sure region_number_col is int
    # assert gr.as_df()[region_number_col].dtype

    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    gr = gr.apply(
        lambda df: df.sort_values(by=[id_col, region_number_col], ascending=True)
    )

    # Filter out single-exon transcripts
    if filter_single:
        print("Filtering for multi-exon transcripts...")
        print("Before: {}".format(len(set(gr.as_df()[id_col].tolist()))))

        # Setting to 'False' marks all duplicates as True (so keep these)
        gr = gr.subset(lambda df: df.duplicated(subset=[id_col], keep=False), nb_cpu=1)

        print("After: {}".format(len(set(gr.as_df()[id_col].tolist()))))

    if number_type == "stranded":
        # source = None means that 1 = first region of group regardless of strand
        # Pick last region entry by max region number for each transcript (id_col)
        # Pick first region entry by min region number for each transcript (id_col)

        # keep="last" sets last in ID to 'False' and all others true (negate to keep last only)
        # keep="first" sets first in ID to 'False'

        out_gr = gr.subset(
            lambda df: ~(df.duplicated(subset=[id_col], keep=which_region)), nb_cpu=1
        )

    else:
        # Numbering doesn't respect strand
        # Need to flip selecting first/last in group depending on strand
        # minus strand - pick min if Minus strand, max if plus strand

        if which_region == "first":
            # + strand - pick first in group, - strand - pick last in group

            out_gr = gr.subset(
                lambda df:
                # 1. plus strand & first in group/ID
                (df["Strand"] == "+") & ~(df.duplicated(subset=[id_col], keep="first"))
                |
                # 2. minus strand & last in group/ID
                (df["Strand"] == "-") & ~(df.duplicated(subset=[id_col], keep="last")),
                nb_cpu=1,
            )

        elif which_region == "last":
            # + strand - pick last in group/ID
            # - strand - pick first in group/ID
            out_gr = gr.subset(
                lambda df:
                # 1. plus strand & last in group/ID
                (df["Strand"] == "+") & ~(df.duplicated(subset=[id_col], keep="last"))
                |
                # 2. minus strand & first in group/ID
                (df["Strand"] == "-") & ~(df.duplicated(subset=[id_col], keep="first")),
                nb_cpu=1,
            )

    return out_gr.sort()


def get_internal_regions(
    gr,
    feature_col="Feature",
    feature_key="exon",
    id_col="transcript_id",
    region_number_col="exon_number",
):
    """
    Return gr of internal exons for each transcript_id
    In process, exon_number_col will be converted to type 'int'
    """

    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [
        feature_key
    ], "only {} entries should be present in gr".format(feature_key)

    # Pull out exons, convert exon_number to int
    exons_gr = gr.assign(
        region_number_col,
        lambda df: df[region_number_col].astype(float).astype("Int64"),
        nb_cpu=1,
    )

    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    exons_gr = exons_gr.apply(
        lambda df: df.sort_values(by=[id_col, region_number_col], ascending=True),
        nb_cpu=1,
    )

    # Filter out 1st + last exons for each ID
    # first exons for each transcript (.ne(1))
    # keep="last" sets last dup value to 'False' & all others True
    # This will filter out last exons

    out_gr = exons_gr.subset(
        lambda df: (df[region_number_col].ne(1).astype(bool))
        & (df.duplicated(subset=["transcript_id"], keep="last")),
        nb_cpu=1,
    )

    return out_gr.sort()


def main(ref_gtf: str, ext_bed: str, output_prefix: str):
    """Extend annotated last exons based on provided extension coordinates

    Parameters
    ----------
    ref_gtf : str
        Path to input (RefSeq) GTF file containing annotated transcript models to extend at 3'end
    ext_bed : str
        Path to BED file containing last exon coordinates of 3'UTR extensions
    output_prefix : str
        Prefix for output files
    """

    # Read the GTF file using pyranges
    logging.info("Reading the GTF file...")
    gtf = pr.read_gtf(ref_gtf)

    # Read the BED file using pyranges
    logging.info("Reading the BED file...")
    bed = pr.read_bed(ext_bed)

    # Workflow:
    # 1. Separate regions to update from those staying as in annotation
    # 1a. Extract last exons & 3'UTRs from GTF into a combined pyranges object.
    # 1b. Also create a separate object excluding these events
    # 2. left join the annotated LEs and 3'UTRs with bed
    # 3. Update the 3'end coordinates of annotated les to the joined coordinates
    # 4. Output to file

    # Step 1a: Extract last exons & 3'UTRs from GTF into a combined pyranges object
    logging.info("Extracting last exons and 3'UTRs from GTF...")
    exons = gtf[gtf.Feature == "exon"].apply(lambda df: df.astype({"exon_number": int}))
    last_exons = get_terminal_regions(exons, which_region="last")
    three_utrs = gtf[gtf.Feature == "3UTR"]

    three_ends = pr.concat([last_exons, three_utrs])
    logging.info("Separating out other intervals which are not being updated...")

    # Step 1b: Create object lacking any region found in three_ends
    # extract all non 3UTR/exon feature entries
    other_feats = gtf[~gtf.Feature.isin({"3UTR", "exon"})]
    nonlast_exons = pr.concat(
        [get_terminal_regions(exons, which_region="first"), get_internal_regions(exons)]
    )

    other_gtf = pr.concat([other_feats, nonlast_exons])

    # Step 2. Left join (pr.join) annotated LEs and 3'UTRs with bed
    logging.info(
        "Updating 3'end coordinates of last exons/3'UTRs to provided 3'UTR extension ends..."
    )
    bed_columns = bed.columns.tolist()
    joined = three_ends.join(bed, strandedness="same")

    # Step 3. Update 3'end coordinates of joined using _df_update_3p
    updated_joined = joined.apply(lambda df: _df_update_3p(df))

    # Drop any joined columns
    logging.info(
        "Return to original configuration (non-updated regions, column names/attributes)..."
    )
    only_bed_columns = [col for col in bed_columns if col not in three_ends.columns]
    updated_joined = updated_joined.drop(like="_b$").drop(only_bed_columns)

    # combine updated entries with remaining other entries
    updated_gtf = pr.concat([other_gtf, updated_joined]).sort()

    logging.info(
        f"Output original 3'ends - {output_prefix + '.original_3p_ends.gtf'} - and updated complete GTF - {output_prefix + '.updated_3p_ends.gtf'} - to file..."
    )
    three_ends.to_gtf(output_prefix + ".original_3p_ends.gtf")
    updated_gtf.to_gtf(output_prefix + ".updated_3p_ends.gtf")

    logging.info("Script execution completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extend annotated last exons based on provided extension coordinates"
    )
    parser.add_argument(
        "-g",
        "--gtf",
        type=str,
        help="Path to input (RefSeq) GTF file containing annotated transcript models to extend at 3'end",
    )
    parser.add_argument(
        "-b",
        "--bed",
        type=str,
        help="Path to BED file containing last exon coordinates of 3'UTR extensions",
    )
    parser.add_argument("-o", "--output", type=str, help="Prefix for output GTF files")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Print the parsed command line argument values
    logging.info("Input args: %r", args)

    main(args.gtf, args.bed, args.output)
