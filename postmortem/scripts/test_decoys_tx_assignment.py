import pandas as pd
import argparse
import sys

'''
Quick and dirty test script to compare the transcript_id to le_id assignment between the original PAPA implementation and the decoy-augmented implementation (add_decoys_to_gtf.py).
Interested in at least one match, because decoys script implements collapsing of redundant intervals to a single representative ID
Lack of matches indicate that le_ids have been shifted/are misaligned
Script raises an exception if any failure encountered (completely missing le_ids and/or misannotated le_ids)
Script can be tested using output of test_data_decoys_tx_assignment.py (LE5 dropped, LE2,4,6 have no matching transcript_ids)
'''

def compare_tsv_files(orig_file: str, new_file: str):
    # Read the TSV files
    tx2le_orig = pd.read_csv(orig_file, sep='\t')
    tx2le_new = pd.read_csv(new_file, sep='\t')

    # Perform a left join on 'le_id'
    merged = pd.merge(tx2le_orig, tx2le_new, on='le_id', how='left', suffixes=('_orig', '_new'))

    # Check for dropped le_ids (will have NaN in transcript_id_new)
    dropped_le_ids = merged[merged['transcript_id_new'].isna()]['le_id']

    # For non-dropped le_ids, check if there's at least one matching transcript_id
    def check_transcript_match(group):
        return (group['transcript_id_orig'] == group['transcript_id_new']).any()

    # Check for >= 1 matching tx_ids 
    # returns series with le_id as index, values = True/False if matching/no matching transcript ids
    transcript_match = merged[~merged['transcript_id_new'].isna()].groupby('le_id').apply(check_transcript_match)

    # Identify le_ids with no matching transcript_ids
    unmatched_transcripts = transcript_match[~transcript_match].index

    # Print results
    print("Results of the comparison:")

    dropped_ids = False
    print(f"1. Number of le_ids in original file: {len(tx2le_orig['le_id'].unique())}")
    print(f"2. Number of dropped le_ids: {len(dropped_le_ids)}")
    if len(dropped_le_ids) > 0:
        dropped_ids = True
        print("Dropped le_ids:")
        print(dropped_le_ids.tolist())

    unmatched_tx = False
    print(f"3. Number of le_ids with no matching transcript_ids: {len(unmatched_transcripts)}")
    if len(unmatched_transcripts) > 0:
        unmatched_tx = True
        print("le_ids with no matching transcript_ids:")
        print(unmatched_transcripts.tolist())

    print("\nSummary:")
    if len(dropped_le_ids) == 0 and len(unmatched_transcripts) == 0:
        print("All le_ids are present in both files and have at least one matching transcript_id.")
    else:
        print("There are discrepancies between the files. Please review the details above.")
        raise Exception(f"Following errors were encountered (True/False): dropped le_ids = {dropped_ids} , no matching transcript_ids for same le_id = {unmatched_tx}")


def main(args):

    compare_tsv_files(args.orig_tx2le, args.new_tx2le)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check for consistent assignment of transcript_ids to le_ids between two tx2le tables.")
    parser.add_argument("orig_tx2le", help="Path to the original 'tx2le' TSV file")
    parser.add_argument("new_tx2le", help="Path to the new 'tx2le' TSV file")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    main(args)



# # Read the TSV files
# tx2le_orig = pd.read_csv('tx2le_orig.tsv', sep='\t')
# tx2le_new = pd.read_csv('tx2le_new.tsv', sep='\t')

# # Perform a left join on 'le_id'
# merged = pd.merge(tx2le_orig, tx2le_new, on='le_id', how='left', suffixes=('_orig', '_new'))

# # Check for dropped le_ids (will have NaN in transcript_id_new)
# dropped_le_ids = merged[merged['transcript_id_new'].isna()]['le_id']

# # For non-dropped le_ids, check if there's at least one matching transcript_id
# def check_transcript_match(group):
#     return (group['transcript_id_orig'] == group['transcript_id_new']).any()


# transcript_match = merged[~merged['transcript_id_new'].isna()].groupby('le_id').apply(check_transcript_match)

# # Identify le_ids with no matching transcript_ids
# unmatched_transcripts = transcript_match[~transcript_match].index

# # Print results
# print("Results of the comparison:")
# print(f"1. Number of le_ids in original file: {len(tx2le_orig['le_id'].unique())}")
# print(f"2. Number of dropped le_ids: {len(dropped_le_ids)}")



# print("\nSummary:")
# if len(dropped_le_ids) == 0 and len(unmatched_transcripts) == 0:
#     print("All le_ids are present in both files and have at least one matching transcript_id.")
# else:
#     print("There are discrepancies between the files. Please review the details above.")
    
