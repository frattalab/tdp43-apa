import pandas as pd

# Create test data for tx2le_orig
tx2le_orig_data = {
    'le_id': ['LE1', 'LE2', 'LE3', 'LE4', 'LE5', 'LE6'],
    'transcript_id': ['T1', 'T2', 'T3', 'T4', 'T5', 'T6']
}
tx2le_orig = pd.DataFrame(tx2le_orig_data)

# Create test data for tx2le_new
tx2le_new_data = {
    'le_id': ['LE1', 'LE2', 'LE3', 'LE4', 'LE6', 'LE7'],
    'transcript_id': ['T1', 'T2_new', 'T3', 'T4_new', 'T6_new', 'T7']
}
tx2le_new = pd.DataFrame(tx2le_new_data)

# Save test data to TSV files
tx2le_orig.to_csv('tx2le_orig.tsv', sep='\t', index=False)
tx2le_new.to_csv('tx2le_new.tsv', sep='\t', index=False)

print("Test data created and saved to tx2le_orig.tsv and tx2le_new.tsv")

# Display the test data
print("\ntx2le_orig:")
print(tx2le_orig)
print("\ntx2le_new:")
print(tx2le_new)