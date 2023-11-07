import pyranges as pr
import pandas as pd
from Bio.Seq import Seq


def get_terminal_regions(gr: pr.PyRanges,
                         feature_col: str = "Feature",
                         feature_key: str = "exon",
                         id_col: str = "transcript_id",
                         region_number_col: str = "exon_number",
                         which_region: str = "last",
                         filter_single: bool = False,
                         ):
    '''
    Return gr of last exons for each transcript_id
    In process, region_number_col will be converted to type 'int'

    i.e. for minus-strand - largest exon_number for transcript corresponds to FIRST EXON, not last
    Annotated (i.e. Ensembl) reported exon_numbers DO RESPECT STRAND (i.e. max always = last exon)

    if Do respect strand, put source = None (default)
    if Don't respect strand, put source = "stringtie" (i.e. plus strand = max, minus strand = min)

    which_region can be set to 'last' or 'first'
    '''

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


    # 1 = first region of group regardless of strand
    # Pick last region entry by max region number for each transcript (id_col)
    # Pick first region entry by min region number for each transcript (id_col)

    # keep="last" sets last in ID to 'False' and all others true (negate to keep last only)
    # keep="first" sets first in ID to 'False'

    out_gr = mod_gr.subset(lambda df: ~(df.duplicated(subset=[id_col], keep=which_region)),
                                                          nb_cpu=1
                              )


    return out_gr



def get_internal_regions(gr: pr.PyRanges ,
                         feature_col: str = "Feature",
                         feature_key: str = "exon",
                         id_col: str ="transcript_id",
                         region_number_col: str = "exon_number",
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


def translate_exon(gr: pr.PyRanges, fasta: str = "data/GRCh38.primary_assembly.genome.fa", drop_seqs: bool = True, add_stop: bool = True) -> pr.PyRanges:
    '''Generate a peptide sequence for a continuous genomic interval (e.g. exon)

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    fasta : str, optional
        path to FASTA file of genome sequence, by default "data/GRCh38.primary_assembly.genome.fa"
    drop_seqs : bool, optional
        whether to drop columns containing extracted genomic sequence, by default True

    Returns
    -------
    pr.PyRanges
        input PyRanges object with added 'peptide_seq' column containing str of peptide sequences
    '''

    assert "Frame" in gr.columns
    assert "dna_seq" not in gr.columns

    if add_stop:
        stop_str = "*"
    else:
        stop_str = ""

    # extract genomic sequences for intervals, add as column
    seqs_gr = pr.get_sequence(gr, fasta)
    gr.dna_seq = seqs_gr

    # convert dna_seq to a Seq object to translate (BioPython)
    # remove nucleotides from start if necessary to get complete codon
    gr = gr.assign("seq_dna_seq",
                    lambda df: df.apply(lambda row: Seq(row["dna_seq"][int(row["Frame"]):]),
                                        axis=1)
                                        )
    
    # generate peptide sequence
    gr = gr.assign("peptide_seq",
                   lambda df: df["seq_dna_seq"].apply(lambda x: str(x.translate(to_stop=True)) + stop_str)
                   )

    if drop_seqs:
        return gr.drop(["dna_seq", "seq_dna_seq"])
    
    else:
        return gr
    
    
def translate_tx(gr: pr.PyRanges, grpby: str, fasta: str = "data/GRCh38.primary_assembly.genome.fa", drop_seqs: bool = True, add_stop: bool = False) -> pd.DataFrame:
    '''Generate a peptide sequence for a combination of intervals (e.g. transcripts)

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    grpby : str
        column containing group identifiers
    fasta : str, optional
        path to FASTA file of genome sequence, by default "data/GRCh38.primary_assembly.genome.fa"
    drop_seqs : bool, optional
        whether to drop columns containing extracted genomic sequence, by default True
    add_stop : bool, optional
        whether to add '*' symbol to end of inferred peptide sequence to indicate stop codon, by default False

    Returns
    -------
    pd.DataFrame
        _description_
    '''
    
    assert "Frame" in gr.columns
    assert "exon_number" in gr.columns
    assert gr.exon_number.dtype == "int64"

    if add_stop:
        stop_str = "*"
    else:
        stop_str = ""
    
    # get the frame of the first interval of each group
    tx2frame = gr.as_df()[[grpby, "exon_number", "Frame"]].astype({"Frame": int}).sort_values(by=[grpby, "exon_number"]).drop_duplicates(subset=grpby, keep="first")

    # get a df of transcript sequence for each interval
    tr_seqs = pr.get_transcript_sequence(gr, group_by=grpby, path=fasta)

    # merge back in frame for each group
    tr_seqs = tr_seqs.merge(tx2frame, on=grpby, how="left")

    # convert sequences to Seq object, trimming from start of sequence depenedent on frame
    tr_seqs["seq"] = tr_seqs.apply(lambda row: Seq(row["Sequence"][row["Frame"]:]), axis=1)

    # generate peptide sequence
    tr_seqs["peptide_seq"] = tr_seqs["seq"].apply(lambda x: str(x.translate(to_stop=True)) + stop_str)
    
    if drop_seqs:
        return tr_seqs.drop(columns=["Sequence", "seq"])
    
    else:  
        return tr_seqs
    

# function to match two peptide sequences & find position where they stop matching (assuming they start at the same position)
def longest_matching_substring(str1: str, str2: str) -> int:
    '''Return the index of the final position in longest exactly matched substring from the beginning of two strings

    this is intended for strings that begin with identical values. it only checks the longest matching substring from the first position in each string. if there is one mismatch, the function terminates
    

    Parameters
    ----------
    str1 : str
        first string of pair wish to find the longest matching substring from the beginning of the string
    str2 : str
        second string of pair wish to find the longest matching substring from the beginning of the string

    Returns
    -------
    int
        index of the final position (in either string) of the longest matched substring between the two. If you wish to slice a string to retain/exclude the longest match, remember to add one to returned value
    '''
    longest_substring = 0 
    for i, (char1, char2) in enumerate(zip(str1, str2)):
        if char1 == char2:
            longest_substring = i
        else:
            break

    return longest_substring 

# string1 = "abcdefg"
# string2 = "abcxyz"

# print("Testing with strings:", string1, string2)
# print("Original function result:", longest_matching_substring(string1, string2))
# print("Original function result extracting matched string1:", string1[:longest_matching_substring(string1, string2) + 1])
# print("Original function result extracting unique part of string1:", string1[longest_matching_substring(string1, string2) + 1:])