# %%
#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
from typing import Union, List, Optional, Tuple, Dict, Set
from pathlib import Path
from datetime import datetime
import os
import math
import multiprocessing

def nearest_threshold_wrapper(gr: pr.PyRanges,
                              gr2: pr.PyRanges,
                              thresholds: Union[int, List[int]],
                              id_col: str = "le_id",
                              nearest_only: bool = True,
                              return_grs: bool = False,
                              nearest_kwargs: dict = {"k": 10,
                                                      "strandedness": "same",
                                                      "overlap": True,
                                                      "how": None},
                              ) -> Tuple[Dict[str, set], Optional[Dict[str, pr.PyRanges]]]:
    '''Get events with nearest overlap passing a maximum distance threshold(s)

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    gr2 : pr.PyRanges
        _description_
    thresholds : Union[int, List[int]]
        _description_
    id_col : str, optional
        _description_, by default "le_id"
    nearest_only : bool, optional
        _description_, by default True
    return_grs : bool, optional
        _description_, by default False
    nearest_kwargs : _type_, optional
        _description_, by default {"k": 10, "strandedness": "same", "overlap": True, "how": None}

    Returns
    -------
    Tuple[Dict[str, set], Optional[Dict[str, pr.PyRanges]]]
        return {threshold: set(id_col values)} and optionally {threshold: pr.PyRanges} (if return_grs=True)
    '''
    

    if nearest_only:
        if "k" in nearest_kwargs.keys():
            kwargs = {key: v for key, v in nearest_kwargs.items()if key != "k"}
        gr_nr = gr.nearest(gr2, **kwargs)

    else:
        assert "k" in nearest_kwargs.keys()
        gr_nr = gr.k_nearest(gr2, **nearest_kwargs)


    # convert to absolute values (don't care if up/downstream rn)
    gr_nr = gr_nr.assign("DistanceAbs", lambda df: df.Distance.abs())
    
    # if nearest_only:
    #     # # for each le_id, select the smallest distance
    #     gr_nr = gr_nr.apply(lambda df: df.sort_values(by=[id_col, "DistanceAbs"]).drop_duplicates(subset=[id_col])).sort()

    if isinstance(thresholds, int):
        # subset to those pasing distance threshold
        pass_grs = {str(thresholds): gr_nr.subset(lambda df: df.DistanceAbs.le(thresholds))}

    else:
        pass_grs = {str(threshold): gr_nr.subset(lambda df: df.DistanceAbs.le(threshold)) for threshold in thresholds}

    # return set of ids at each threshold
    pass_ids = {threshold: set(gr.as_df()[id_col]) for threshold, gr in pass_grs.items()}

    if return_grs:
        return pass_ids, pass_grs
    else:
        return pass_ids


def create_threshold_binary_dataframe(
    ids: Set[str], 
    nr_threshold_ids: Dict[str, Set[str]], 
    id_column: str = 'le_id', 
    column_prefix: str = 'distancethresh_'
) -> pd.DataFrame:
    '''Create df summarising the results of nearest_threshold_wrapper() into a pass/fail df for a set of ids at each provided cutoff

    Parameters
    ----------
    ids : Set[str]
        set of IDs of interest
    nr_threshold_ids : Dict[str, Set[str]]
        output of nearest_threshold_wrapper()
    id_column : str, optional
        _description_, by default 'le_id'
    column_prefix : str, optional
        _description_, by default 'distancethresh_'

    Returns
    -------
    pd.DataFrame
        _description_
    '''
    # Convert list of all ids to DataFrame
    le_ids_df = pd.DataFrame({id_column: list(ids)})
    
    # Add columns for each key in the dictionary, prefixed with 'distance_'
    for key, ids in nr_threshold_ids.items():
        col_name = f'{column_prefix}{key}'
        le_ids_df[col_name] = le_ids_df[id_column].apply(lambda x: 1 if x in ids else 0)
    
    return le_ids_df


def create_threshold_count_dataframe(
    le_ids: Set[str], 
    nr_threshold_grs: Dict[str, pr.PyRanges], 
    id_column: str = 'le_id', 
    nr_id_column: str = 'Name',
    column_prefix: str = 'distancethresh_'
) -> pd.DataFrame:
    '''For a given set of target IDs, count the number of unique intervals that fall within a (series of) distance threshold

    Parameters
    ----------
    le_ids : Set[str]
        set of IDs of interest
    nr_threshold_grs : Dict[str, pr.PyRanges]
        output of nearest_threshold_wrapper(return_grs=True) - grs containing intervals passing distance threshold
    id_column : str, optional
        _description_, by default 'le_id'
    column_prefix : str, optional
        _description_, by default 'distance_'
    nr_id_column : str, optional
        identifier column for 'joined nearest intervals' (i.e. gr2 from nearest_threshold_wrapper), by default 'Name'

    Returns
    -------
    pd.DataFrame
        _description_
    '''

    # Convert list of all ids to DataFrame
    le_ids_df = pd.DataFrame({id_column: list(le_ids)})
    
    # Add columns for each key in the dictionary, prefixed with 'distance_'
    for key, gr in nr_threshold_grs.items():
        col_name = f'{column_prefix}{key}'
        
        # Group by the specified id column and count unique 'nr_id_column' values
        gr_df = gr.as_df().groupby(id_column)[nr_id_column].nunique().reset_index()
        
        # Merge the counts with the main DataFrame (missing = none within distance, so fill with 0)
        le_ids_df = le_ids_df.merge(gr_df, on=id_column, how='left').fillna(0)
        le_ids_df = le_ids_df.rename(columns={nr_id_column: col_name})
    
    return le_ids_df.astype({f'{column_prefix}{k}': int for k in nr_threshold_grs.keys()})


# Split the file list into evenly sized batches
def split_file_list(file_list, num_batches):
    batch_size = math.ceil(len(file_list) / num_batches)
    return [file_list[i:i + batch_size] for i in range(0, len(file_list), batch_size)]


# run nearest_threshold_wrapper over list of text files containing IDs to evaluate (intended for matched annotated ids)
def annot_ids_path_nearest_threshold_wrapper(file_names: list, annot_ids_dir: str,
                                             gr: pr.PyRanges,
                                             gr2: pr.PyRanges,
                                             distance_thresholds: Union[int, List[int]],
                                             progress_counter):
    
    
    results = {}
    
    for i, txt_file in enumerate(file_names, 1):
        # Strip file extension to get identifier
        file_identifier = Path(txt_file).stem
        
        # Read identifiers into a set
        with open(os.path.join(annot_ids_dir, txt_file), 'r') as f:
            identifiers = set(line.strip() for line in f)
        
        # Filter papa_pas for rows corresponding to these identifiers
        gr_filt = gr[papa_pas_annot.le_id.isin(identifiers)]
        
        # double check all IDs found
        if not set(gr_filt.le_id) == identifiers:
            print(f"Not all identifiers found in papa_pas for file {file_identifier}")
            continue
        
        # Call nearest_threshold_wrapper
        wrapper_dict = nearest_threshold_wrapper(
            gr_filt, 
            gr2, 
            distance_thresholds, 
            nearest_only=True, 
            return_grs=False
        )
        
        # Call create_threshold_binary_dataframe
        threshold_df = create_threshold_binary_dataframe(identifiers, wrapper_dict)
        
        # Store results
        results[file_identifier] = threshold_df

        # Update progress counter and print progress every n files files
        progress_counter.value += 1
        if progress_counter.value % 25 == 0:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Progress: {progress_counter.value} files processed", end="\n")
    
    return results


# Combine results from all processes into a single dictionary
def combine_results(all_results):
    combined = {}
    for result in all_results:
        combined.update(result)
    return combined


# Main function to execute parallel processing
def parallel_annot_ids_nearest_threshold_wrapper(file_list: List[str],
                                                 annot_ids_dir: str,
                                                 gr: pr.PyRanges,
                                                 gr2: pr.PyRanges,
                                                 distance_thresholds: Union[int, List[int]],
                                                 num_processes: int = 1):
    
    # Split the file list into batches
    file_batches = split_file_list(file_list, num_processes)

    # Create a manager to handle shared resources like the progress counter
    with multiprocessing.Manager() as manager:
        # Create a shared counter for tracking progress
        progress_counter = manager.Value('i', 0)  # 'i' for integer, initialized to 0
        print(f"Total files (sets of annotated IDs) to process - {len(file_list)}")
        
        # Use multiprocessing pool
        with multiprocessing.Pool(processes=num_processes) as pool:
            # Prepare the arguments for the worker
            pool_args = [(file_batch, annot_ids_dir, gr, gr2, distance_thresholds, progress_counter)
                         for file_batch in file_batches]
            
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Beginning processing of files")
            # Use `map` to process the batches
            batch_results = pool.starmap(annot_ids_path_nearest_threshold_wrapper, pool_args)

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Processing complete: Combining outputs across batches")

        combined_results = {}
        for result in batch_results:
            combined_results.update(result)
        
        return combined_results
    
    # # Create a pool of processes
    # with multiprocessing.Pool(processes=num_processes) as pool:
    #     # Process each batch in parallel and gather results
    #     all_results = pool.starmap(annot_ids_path_nearest_threshold_wrapper,
    #                                [(batch, annot_ids_dir, gr, gr2, distance_thresholds) 
    #                                 for batch in file_batches]
    #                                 )
    
    # # Combine the results from all processes into a single dictionary
    # final_results = combine_results(all_results)
    # return final_results

# %%
# last exon annotations (to extract pas)
papa_pas = pr.read_gtf("../data/novel_ref_combined.last_exons.gtf")[["le_id", "ref_gene_name", "ref_gene_id"]].three_end()
papa_pas = papa_pas.drop_duplicate_positions(strand=True)
outdir = "../processed/curation/cryptic_annot_comparison"
pas_counts = pd.read_csv(os.path.join(outdir, "2024-09-03_le_id_pas_counts.tsv"), sep="\t")
annot_ids_dir = "../processed/curation/cryptic_annot_comparison/ids/"
distance_thresholds = [0,10,25,50,100,200,500]
# number of processes for parallel processing of annotated events
nproc = 8
papa_pas

# %%
# PATRs - pooling across all KD samples from all experiments
patr_pas_all = pr.read_bed("../data/bulk_polya_reads/tdp_ko_collection/pas_clusters/condition__TDP43KD/two_class_simple/polya_clusters.bed")
patr_pas_all

# %%
# extract cryptic PAS
cryp_le_ids = set(pas_counts.loc[pas_counts.cryptic_status.eq(1), "le_id"])
papa_pas_cryp = papa_pas.subset(lambda df: df.le_id.isin(cryp_le_ids))
papa_pas_cryp

# %%
# identify PAS support for cryptics at range of distance thresholds
# # get ids & grs where cryptics have pas junction support from any dataset across range of thresholds
cryp_patr_all_nr_ids, cryp_patr_all_nr_grs = nearest_threshold_wrapper(papa_pas_cryp, patr_pas_all, distance_thresholds, nearest_only=False, return_grs=True)
cryp_patr_all_nr_df = create_threshold_binary_dataframe(cryp_le_ids, cryp_patr_all_nr_ids)
cryp_patr_all_nr_df

# %%
# calculate sum / num cryptic ids supported at each threshold (column)
cryp_patr_all_nr_df.drop(columns="le_id").sum(axis=0)


# %%

    
# List all .txt files in the directory
txt_files = [f for f in os.listdir(annot_ids_dir) if f.endswith('.txt')]
total_files = len(txt_files)

# first, collect all annotated IDs plan to extract, then subset the papa PAS GTF
# (saves lugging around a huge pyranges object, possibly speeds up the subsetting)
all_annot_le_ids = set()
for file in txt_files:
    with open(os.path.join(annot_ids_dir, file), 'r') as f:
        for line in f:
            all_annot_le_ids.add(line.rstrip())

print(f"Total unique annotated le_ids to analyse - {len(all_annot_le_ids)}")
print(f"Number of intervals in PAPA PAS GTF - {len(papa_pas)}")
papa_pas_annot = papa_pas.subset(lambda df: df.le_id.isin(all_annot_le_ids))

# %%
# Calculate for all iterations of sampled, covariate-matched (expression, num PAS) annotated events
results = parallel_annot_ids_nearest_threshold_wrapper(txt_files, annot_ids_dir, papa_pas_annot, patr_pas_all, distance_thresholds, num_processes=nproc)


# %%
# output binary dataframes at each threshold for cryptics and combined annotated iterations

# combine annotated dataframes across iterations
comb_annot_all_nr_df = pd.concat({k: v for k,v in results.items()}, names=['iteration']).reset_index(level='iteration')
comb_annot_all_nr_df

# %%
# 227 * 1000 is going to create a 227k row CSV/TSV (quite big)
# Main interest is in number of events supported at a given threshold
# Precaluclate counts of supported events at each threshold for each iteration
distancethresh_cols = [col for col in comb_annot_all_nr_df.columns if col.startswith('distancethresh')]
comb_annot_all_nr_df_counts = comb_annot_all_nr_df.groupby('iteration')[distancethresh_cols].sum().reset_index()
comb_annot_all_nr_df_counts

# %%
# repeat for cryptic events
cryp_patr_all_nr_df_counts = cryp_patr_all_nr_df.drop(columns="le_id").sum(axis=0)
cryp_patr_all_nr_df_counts
# %%
# output counts to plain TSV
comb_annot_all_nr_df_counts.to_csv(os.path.join(outdir, "2024-09-13_annotated_iterations.patr_all.nearest_threshold.counts.tsv"), sep="\t", header=True, index=False)
(cryp_patr_all_nr_df_counts
 # convert to df
 .reset_index(name="n").rename(columns={"index": "distancethresh"})
 .to_csv(os.path.join(outdir, "2024-09-13_cryptics.patr_all.nearest_threshold.counts.tsv"), sep="\t", header=True, index=False)
 )

# output binary matrices of diff IDs to TSV (and GZIPPED TSV for annotated)
comb_annot_all_nr_df.to_csv(os.path.join(outdir,
                                         "2024-09-13_annotated_iterations.patr_all.nearest_threshold.matrix.tsv.gz"), sep="\t", header=True, index=False)
cryp_patr_all_nr_df.to_csv(os.path.join(outdir,
                                         "2024-09-13_cryptics.patr_all.nearest_threshold.matrix.tsv"), sep="\t", header=True, index=False)


