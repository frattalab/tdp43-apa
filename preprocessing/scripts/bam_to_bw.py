#!/usr/bin/env python3

import pyranges as pr
import pysam
import pyrle
from collections import namedtuple
from typing import Tuple, Dict, List, Literal
import argparse
import sys

# object to store different outputs of bam_to_strand_alignments
Alignments = namedtuple('Alignments', 
                        field_names=["gr", "fragments_plus", "fragments_minus", "refskip_idxs", "orphans_plus", "orphans_minus"],
                        defaults=[None]*6, 
                        )

# Object to store region coordinates to provide to pysam.AlignmentFile.fetch()
# Note: strand is not passed to pysam.AlignmentFile.fetch(), but can be used when trying to match/filter alignments for a given strand
# None is default for not specifying a region/contig
Region = namedtuple("Region",
                    field_names=["contig", "start", "stop", "strand"],
                    defaults=[None]*4)

def get_chromsizes_dict(path: str) -> dict:
    '''Read in a chromsizes txt file to a dictionary of contig names and length

    Parameters
    ----------
    path : str
        _description_

    Returns
    -------
    dict
        _description_
    '''

    with open(path) as infile:
        chromsizes = {line.split('\t')[0]: line.split('\t')[1].strip() for line in infile}

    return chromsizes


def cigar_to_gap_idxs(lst: List[int]) -> List[int]:
    '''Get indexes of reference skips with respect to all gaps in an alignment

    Parameters
    ----------
    lst : List[int]
        operations returned by pysam.AlignedSegment.cigartuples e.g. [op for op,_ in pysam.AlignedSegment.cigartuples]

    Returns
    -------
    List[int]
        indexes of ref-skips in alignment with respect to all ref-skips/deletions (= gap in alignment with respect to reference)
    '''
    
    # pysam docs - N (ref-skip) = 3, D (del with respect to reference) = 2

    # Use a list comprehension to find the indexes where the value is 3 or 2
    indexes = [i for i, value in enumerate(lst) if value == 3 or value == 2]

    # grab the indexes of 3/ref-skip with respect to all gaps in the read alignment (3/2)
    filtered_indexes = [i for i, lst_idx in enumerate(indexes) if lst[lst_idx] == 3]

    return filtered_indexes


def check_aln_strand(read: pysam.AlignedSegment,
                     strandedness: str = Literal["--rf", "--fr"]):
    '''Determined transcribed strand of origin for a read given strandedness of library prep

    Parameters
    ----------
    read : pysam.AlignedSegment
        _description_
    strandedness : str, optional
        _description_, by default Literal["--rf", "--fr"]

    Returns
    -------
    _type_
        _description_
    '''
    
    assert strandedness in ["--rf", "--fr"]

    # determine alignment strand
    if strandedness == "--rf":
        # reverse stranded - if read 1 then maps to strand of origin, otherwise opposite strand
        if read.is_reverse:
            # maps to minus strand - if read1 then opposite strand of origin, otherwise same strand
            if read.is_read1:
                return 1
            else:
                return -1
    
        # read maps to plus strand - if read1 then opposite strand of origin, otherwise same strand
        elif read.is_read1:
            return -1
    
        elif read.is_read2:
            return 1
        
    else:
        # forward stranded - if read1 then maps to strand of origin, otherwise opposite
        if read.is_read1:
            return -1 if read.is_reverse else 1
        
        else:
            return 1 if read.is_reverse else -1


def bam_to_strand_alignments(bam_path: str,
                             region: Region,
                             strandedness: str,
                             as_pyranges: bool = True,
                             return_refskips: bool = True,
                             return_orphans: bool = True) -> Alignments:
    '''_summary_

    Parameters
    ----------
    bam_path : str
        _description_
    region : Region
        _description_
    strandedness : str
        _description_
    as_pyranges : bool, optional
        _description_, by default True
    return_refskips : bool, optional
        _description_, by default True
    return_orphans : bool, optional
        _description_, by default True

    Returns
    -------
    Alignments
        _description_
    '''
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    # Create a dictionary to cache read objects (per reference name/chrom)
    # {<read_name>: <pysam.AlignedSegment>}
    read_cache_plus = {}
    read_cache_minus = {}

    # Store aligned regions of read pairs
    fragments_plus = {"Chromosome": [],
                      "Start": [],
                      "End": [],
                      "Strand": [],
                      "read_name": [],
                      "read_name_num": []}

    fragments_minus = {"Chromosome": [],
                       "Start": [],
                       "End": [],
                       "Strand": [],
                       "read_name": [],
                       "read_name_num": []}

    # initialise a dict to store indexes of refskips for reads that have gapped alignments {read_name_num: [<refskip_idxs>]}
    refskips = {}

    # track the number of successfully written read pairs (per reference name/chrom)
    read_count_plus = 0
    read_count_minus = 0

    for alignment in bamfile.fetch(contig=region.contig, start=region.start, stop=region.stop):

        # read must be paired in alignment, and a primary alignment
        if alignment.is_proper_pair and not alignment.is_secondary:

            # determine putative strand of origin
            putative_strand = check_aln_strand(alignment, strandedness)

            if putative_strand == 1:
                # plus strand - check if mate is already in cache
                if alignment.query_name in read_cache_plus:

                    # Pop the mate from the cache
                    mate = read_cache_plus.pop(alignment.query_name)
                    
                    # define read_name_num - additional id of <read_name>_<r1/r2> (to allow read-specific processing)
                    aln_name_suff = "_r1" if alignment.is_read1 else "_r2"
                    mate_name_suff = "_r1" if mate.is_read1 else "_r2"

                    aln_rname_num = alignment.query_name + aln_name_suff
                    mate_rname_num = mate.query_name + mate_name_suff

                    # get indexes of refskips with respect to all gaps in alignment
                    if return_refskips:
                        aln_refskips = cigar_to_gap_idxs([op for op,_ in alignment.cigartuples])
                        mate_refskips = cigar_to_gap_idxs([op for op,_ in mate.cigartuples])

                        # if alignment doesn't have any refskips then don't populate dict
                        if len(aln_refskips) > 0:
                            refskips[aln_rname_num] = aln_refskips
                        
                        if len(mate_refskips) > 0:
                            refskips[mate_rname_num] = mate_refskips
                    
                    # extract aligned coords of read and its mate to a dictionary
                    read_posns = alignment.get_blocks()
                    mate_posns = mate.get_blocks()

                    num_blocks = len(read_posns) + len(mate_posns)


                    # Store start and end coords of aligned segments
                    fragments_plus["Start"].extend(tup[0] for tup in read_posns)
                    fragments_plus["End"].extend(tup[1] for tup in read_posns)

                    # add rname_num for alignment positions
                    fragments_plus["read_name_num"].extend([aln_rname_num] * len(read_posns))

                    # now add mate coords to dict
                    fragments_plus["Start"].extend(tup[0] for tup in mate_posns)
                    fragments_plus["End"].extend(tup[1] for tup in mate_posns)

                    # add rname_num for mate positions
                    fragments_plus["read_name_num"].extend([mate_rname_num] * len(mate_posns))

                    # extend Chromosome, Strand & query name lists
                    fragments_plus["Chromosome"].extend([region.contig] * num_blocks)
                    fragments_plus["Strand"].extend(["+"] * num_blocks)
                    fragments_plus["read_name"].extend([alignment.query_name] * num_blocks)

                    read_count_plus += 1

                else:
                    # Add the alignment to the cache
                    read_cache_plus[alignment.query_name] = alignment

            else:
                # minus strand - check if mate already in cache
                if alignment.query_name in read_cache_minus:

                    # Pop the mate from the cache
                    mate = read_cache_minus.pop(alignment.query_name)

                    # define read_name_num - additional id of <read_name>_<r1/r2> (to allow read-specific processing)
                    aln_name_suff = "_r1" if alignment.is_read1 else "_r2"
                    mate_name_suff = "_r1" if mate.is_read1 else "_r2"

                    aln_rname_num = alignment.query_name + aln_name_suff
                    mate_rname_num = mate.query_name + mate_name_suff

                    # get indexes of refskips with respect to all gaps in alignment
                    if return_refskips:
                        aln_refskips = cigar_to_gap_idxs([op for op,_ in alignment.cigartuples])
                        mate_refskips = cigar_to_gap_idxs([op for op,_ in mate.cigartuples])

                        # if alignment doesn't have any refskips then don't populate dict
                        if len(aln_refskips) > 0:
                            refskips[aln_rname_num] = aln_refskips
                        
                        if len(mate_refskips) > 0:
                            refskips[mate_rname_num] = mate_refskips

                    # extract aligned coords of read and its mate to a dictionary
                    read_posns = alignment.get_blocks()
                    mate_posns = mate.get_blocks()

                    num_blocks = len(read_posns) + len(mate_posns)

                    # Store start and end coords of aligned segments
                    fragments_minus["Start"].extend(tup[0] for tup in read_posns)
                    fragments_minus["End"].extend(tup[1] for tup in read_posns)

                    # add rname_num for alignment positions
                    fragments_minus["read_name_num"].extend([aln_rname_num] * len(read_posns))

                    # now add mate coords to dict
                    fragments_minus["Start"].extend(tup[0] for tup in mate_posns)
                    fragments_minus["End"].extend(tup[1] for tup in mate_posns)
                    
                    # add rname_num for mate positions
                    fragments_minus["read_name_num"].extend([mate_rname_num] * len(mate_posns))

                    # extend Chromosome and Strand lists
                    fragments_minus["Chromosome"].extend([region.contig] * num_blocks)
                    fragments_minus["Strand"].extend(["-"] * num_blocks)
                    fragments_minus["read_name"].extend([alignment.query_name] * num_blocks)

                    read_count_minus += 1

                else:
                    # Add the alignment to the cache
                    read_cache_minus[alignment.query_name] = alignment


    bamfile.close()

    # construct output dict to be unpacked into Alignments namedtuple for output
    out_dict = {}
    
    # convert to pyranges or return dicts
    if as_pyranges:
        out_dict["gr"] = pr.concat([pr.from_dict(fragments_plus, int64=True), pr.from_dict(fragments_minus, int64=True)])
    else:
        out_dict["fragments_plus"] = fragments_plus
        out_dict["fragments_minus"] = fragments_minus

    if return_orphans:
        out_dict["orphans_plus"] = read_cache_plus
        out_dict["orphans_minus"] = read_cache_minus

    if return_refskips:
        out_dict["refskip_idxs"] = refskips

    return Alignments(**out_dict)


def gr_to_rle(gr: pr.PyRanges, merge: bool = True, merge_by: str = "read_name", by_strand: bool = True, by_rpm: bool = False) -> pyrle.RleDict:
    '''_summary_

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    merge : bool, optional
        _description_, by default True
    merge_by : str, optional
        _description_, by default "read_name"
    by_strand : bool, optional
        _description_, by default True
    by_rpm : bool, optional
        _description_, by default False

    Returns
    -------
    pyrle.RleDict
        _description_
    '''

    if merge:
        assert merge_by in gr.columns

        # collapse coverage of intervals of a group (e.g. read pairs) to prevent double-counting
        gr = gr.merge(strand=by_strand, by=merge_by)

    return gr.to_rle(strand=by_strand, rpm=by_rpm)


def rle_to_bigwig(rle, chrom_sizes: dict, out_prefix, by_strand: bool = True, plus_strand_suffix: str = ".plus.bw", minus_strand_suffix: str = ".minus.bw", by_rpm: bool = False, value_col: str = "Score"):
    '''_summary_

    Parameters
    ----------
    rle : _type_
        _description_
    chrom_sizes : _type_
        _description_
    out_prefix : _type_
        _description_
    by_strand : bool, optional
        _description_, by default True
    plus_strand_suffix : str, optional
        _description_, by default ".plus.bw"
    minus_strand_suffix : str, optional
        _description_, by default ".minus.bw"
    '''

    gr = rle.to_ranges()

    if by_strand:
        out_plus = out_prefix + plus_strand_suffix
        out_minus = out_prefix + minus_strand_suffix

        gr["+"].to_bigwig(out_plus, chrom_sizes, rpm=by_rpm, value_col=value_col)
        gr["-"].to_bigwig(out_minus, chrom_sizes, rpm=by_rpm, value_col=value_col)

    else:
        out = out_prefix + ".bw"
        gr.to_bigwig(out, chrom_sizes, rpm=by_rpm, value_col=value_col)


def main(bam_path: str,
         bed_path: str,
         chromsizes_path: str,
         strandedness: str,
         output_prefix: str):


    # read in chromsizes to dict of {chrom: length}
    chromsizes = get_chromsizes_dict(chromsizes_path)
    # print(chromsizes)
    
    # read in regions_bed
    bed = pr.read_bed(bed_path)

    # merge regions by strand to prevent same reads being read in twice
    bed = bed.merge(strand=True)

    # Extract aligned segments for properly mapped read pairs
    # First generate 'Region' objects for each interval in bed file
    regions_dict = bed.apply(lambda df: df.apply(lambda row: Region(contig=row["Chromosome"],
                                                     start=row["Start"], stop=row["End"], strand=row["Strand"]),
                                                     axis=1
                                                     ).tolist(),
                                                    as_pyranges=False)
    
    # print(regions_dict)

    # combine the values (regions for each chrom-strand pair) into a single iterable of Region objects
    regions = []
    for region_list in regions_dict.values():
        regions.extend(region_list)

    num_regions = len(regions)
    print(f"Total number of regions - {num_regions}")
    
    # print(regions)
    
    # Now extract stranded alignments for each specified region (as Alignment objects)
    print("Extracting coverage over specified regions...")
    alns_list = []
    for i, region in enumerate(regions):
        if (i + 1) % 100 == 0:
            print(f"Region number - {i + 1} / {num_regions}")

        alns_list.append(bam_to_strand_alignments(bam_path, region, "--" + strandedness, return_refskips=False))

    # for idx, aln in enumerate(alns_list):
    #     print(f"Number of aligned segments for {regions[idx]}: {len(aln.gr)}")
    
    # print(alns_list)

    # Now generate a combined PyRanges object of aligned blocks for all regions
    # TODO: option to rescue coverage from orphan reads if necessary
    print("Calculating region-wide per base coverage...")
    alns = pr.concat([aln.gr for aln in alns_list])

    print(f"Number of aligned segments for all reads - {len(alns)}")

    # convert aligned blocks to genome-wide coverage, merging covereage by read pairs
    cov = gr_to_rle(alns, merge=True)

    # print(cov)
    
    # output to bigwig
    print("Outputting coverage per strand to bigwig...")
    rle_to_bigwig(cov, chromsizes, output_prefix)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Extract BigWig files containing stranded RNA-seq coverage for specified intervals from a BAM file")

    parser.add_argument('-b', '--bam', dest='bam_path', required=True, help="Path to BAM file")
    parser.add_argument('-r', '--regions-bed', dest='regions_bed', required=True, help="Regions for coverage extraction in BED format")
    parser.add_argument('-c', '--chromsizes', dest='chromsizes_path', required=True, help="Path to chromosome sizes file")
    parser.add_argument('-o', '--output-prefix', dest='output_prefix', required=True, help="Prefix for output file names")
    parser.add_argument(
        '-s', '--strandedness',
        choices=['rf', 'fr'],
        dest='strandedness',
        required=True,
        help="Orientation of library in 'StringTie' convention. rf or fr"
    )


    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()


    print("BAM file:", args.bam_path)
    print("Regions BED file:", args.regions_bed)
    print("Chromosome sizes file:", args.chromsizes_path)
    print("Output prefix:", args.output_prefix)
    print("Strandedness:", args.strandedness)

    main(args.bam_path,  args.regions_bed,  args.chromsizes_path, args.output_prefix)
