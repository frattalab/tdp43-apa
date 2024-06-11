#!/usr/bin/env python3

import pyranges as pr
import pysam
import pyrle
from collections import namedtuple
from typing import Tuple, Dict, List, Literal
from functools import reduce
from operator import add
import argparse
import sys

# object to store different outputs of bam_to_strand_alignments
Alignments = namedtuple('Alignments', 
                        field_names=["fragments", "orphans", "num_reads"],
                        defaults=[None]*3, 
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


def update_aligned_region(fragments: Dict[str, list], alignment: pysam.AlignedSegment, strand: Literal["+", "-"]) -> None:
    '''Update dictionary with minimal region information for a read alignment

    Parameters
    ----------
    fragments : Dict[str, list]
        _description_
    alignment : pysam.AlignedSegment
        _description_
    strand : Literal[
        _description_

    Raises
    ------
    ValueError
        _description_
    '''
    
    if strand not in ['+', '-']:
        raise ValueError("strand must be either '+' or '-'")

    read_posns = alignment.get_blocks()

    num_blocks = len(read_posns)

    # define read_name_num - additional id of <read_name>_<r1/r2> (to allow read-specific processing)
    aln_name_suff = "_r1" if alignment.is_read1 else "_r2"
    aln_rname_num = alignment.query_name + aln_name_suff

    # Store start and end coords of aligned segments
    fragments["Start"].extend(tup[0] for tup in read_posns)
    fragments["End"].extend(tup[1] for tup in read_posns)

    # add rname_num for alignment positions
    fragments["read_name_num"].extend([aln_rname_num] * len(read_posns))

    # extend Chromosome and Strand lists
    fragments["Chromosome"].extend([alignment.reference_name] * num_blocks)
    fragments["Strand"].extend([strand] * num_blocks)
    fragments["read_name"].extend([alignment.query_name] * num_blocks)


def bam_to_strand_alignments(bam_path: str,
                             region: Region,
                             strandedness: str,
                             as_pyranges: bool = False,
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
        Whether to return paired alignments + orphans as pr.PyRanges (True) or pyrle.Rle object (False), by default True
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
    read_cache = {}

    # Store aligned regions of read pairs (both mates overlapping with interval)
    fragments = {"Chromosome": [],
                      "Start": [],
                      "End": [],
                      "Strand": [],
                      "read_name": []}

    # track the number of successfully written read pairs (per reference name/chrom)
    # read_count_plus = 0


    for alignment in bamfile.fetch(contig=region.contig, start=region.start, stop=region.stop):

        # read must be paired in alignment, and a primary alignment
        if alignment.is_proper_pair and not alignment.is_secondary:

            # determine putative strand of origin
            putative_strand = check_aln_strand(alignment, strandedness)

            if (putative_strand == 1 and region.strand == "+") or (putative_strand == -1 and region.strand == "-"):
                # check if mate is already in cache
                if alignment.query_name in read_cache:
                    # Pop the mate from the cache
                    mate = read_cache.pop(alignment.query_name)
                    
                    # extract aligned coords of read and its mate to a dictionary
                    read_posns = alignment.get_blocks()
                    mate_posns = mate.get_blocks()
                    num_blocks = len(read_posns) + len(mate_posns)

                    # Store start and end coords of current read + mate
                    fragments["Start"].extend(tup[0] for tup in read_posns)
                    fragments["End"].extend(tup[1] for tup in read_posns)

                    fragments["Start"].extend(tup[0] for tup in mate_posns)
                    fragments["End"].extend(tup[1] for tup in mate_posns)

                    # extend Chromosome, Strand & query name lists
                    fragments["Chromosome"].extend([region.contig] * num_blocks)
                    fragments["Strand"].extend([region.strand] * num_blocks)
                    fragments["read_name"].extend([alignment.query_name] * num_blocks)

                else:
                    # Add the alignment to the cache
                    read_cache[alignment.query_name] = alignment

    # finished iterating
    bamfile.close()

    # construct output dict to be unpacked into Alignments namedtuple for output
    out_dict = {}

    # convert to pyranges object
    fragments = pr.from_dict(fragments, int64=True)
    if len(fragments) == 0:
        num_frags = 0
    else:
        num_frags = fragments.read_name.nunique()

    if not as_pyranges:
        # convert to PyRle object
        fragments = gr_to_rle(fragments, merge=True, merge_by="read_name", by_strand=True)

    out_dict["fragments"] = fragments

    if return_orphans:
        # extract info from orphans
        orphans = {"Chromosome": [],
                      "Start": [],
                      "End": [],
                      "Strand": []
                      }
        
        # to add to final total of reads
        num_orphans = len(read_cache)
        
        for alignment in read_cache.values():
            # extract aligned coords of read and its mate to a dictionary
            read_posns = alignment.get_blocks()
            num_blocks = len(read_posns)

            # Store start and end coords of current read + mate
            orphans["Start"].extend(tup[0] for tup in read_posns)
            orphans["End"].extend(tup[1] for tup in read_posns)

            # extend Chromosome, Strand lists (read name irrelevant here as no overlap)
            orphans["Chromosome"].extend([region.contig] * num_blocks)
            orphans["Strand"].extend([region.strand] * num_blocks)

        # convert to pyranges object
        orphans = pr.from_dict(orphans, int64=True)

        if not as_pyranges:
            # convert to PyRle object
            orphans = gr_to_rle(orphans, merge=False, by_strand=True)
        
        out_dict["orphans"] = orphans
        out_dict["num_reads"] = num_frags + num_orphans
    
    else:
        out_dict["num_reads"] = num_frags


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

    if gr.empty:
        return pyrle.RleDict()

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
         output_prefix: str,
         keep_orphans: bool = True):


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
    regions = [region_list for chrom_strand in regions_dict.values() for region_list in chrom_strand]

    num_regions = len(regions)
    print(f"Total number of regions - {num_regions}")
    
    # print(regions)
    
    # Now extract stranded alignments for each specified region (as Alignment objects)
    print("Extracting coverage over specified regions...")
    alns_list = []
    for i, region in enumerate(regions):
        if (i + 1) % 100 == 0:
            print(f"Region number - {i + 1} / {num_regions}")

        alns_list.append(bam_to_strand_alignments(bam_path, region, "--" + strandedness, as_pyranges=False, return_orphans=keep_orphans))


    # Now generate a combined PyRanges object of aligned blocks for all regions
    # TODO: option to rescue coverage from orphan reads if necessary
    print("Done!")
    print(f"Number of extracted reads across all regions - {sum(aln.num_reads for aln in alns_list)}")
    
    print("Obtaining per-base coverage across all provided regions from overlapping pairs...")
    # combine individual pyrle objects for each region into a single PyRle object
    fragments_all = reduce(add, (aln.fragments for aln in alns_list))
    
    if keep_orphans:
        # repeat for orphans
        print("Adding per-base coverage of regions from orphan paired-end reads...")
        orphans_all = reduce(add, (aln.orphans for aln in alns_list))
        
        # Combine coverage
        fragments_all = fragments_all + orphans_all

    print("Outputting coverage per strand to bigwig...")
    rle_to_bigwig(fragments_all, chromsizes, output_prefix)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Extract BigWig files containing stranded RNA-seq coverage of paired-end reads for specified intervals from a BAM file")

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
    parser.add_argument("--keep-orphans", action="store_true", help="Whether to keep reads where mate does not align within specified region (recommended for short intervals)")


    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    print("Input parameters:")
    print("BAM file:", args.bam_path)
    print("Regions BED file:", args.regions_bed)
    print("Chromosome sizes file:", args.chromsizes_path)
    print("Output prefix:", args.output_prefix)
    print("Strandedness:", args.strandedness)
    print("Keep orphans (i.e. properly paired, but only 1 aligning to region):", args.keep_orphans)
    print("Running analysis...")

    main(args.bam_path, args.regions_bed, args.chromsizes_path, args.strandedness, args.output_prefix, args.keep_orphans)
