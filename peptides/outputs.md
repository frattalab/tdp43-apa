## 2023-11-08_i3cortical_cryptic_bleedthrough_fullpeptides.tsv

TSV file containing peptide sequences of bleedthrough-last exon terminating transcripts.

Strategy:

My pipeline defines bleedthrough last exons as:

- Exactly matches an annotated internal exon at its 5'end (i.e. incoming splice junction)
- Last exon 3'end is downstream of the annotated exon's 3'end (i.e. it extends the known exon and terminates/truncates the transcript).

These assumptions allow me to simplify extracting peptide sequences:

- Match the bleedthrough last exon with its respective matching annotated coding exon, replacing the annotated exon 3'end coordinate with the bleedthrough's 3'end coordinate
- Extract all upstream annotated coding exons for the given annotated transcript
- Generate a 'truncated transcript', combining the upstream annotated coding exons with the bleedthrough last exon.
- Translate the truncated transcript (& the full annotated transcript as a sanity check)

Assorted notes:

- 'annotated coding exon' just means annotated CDS sequences from a gencode GTF file (I used v40)
- Stop codons are included in the sequences and marked by *

Column descriptions:

- cryptic_transcript_id - <ensembl_transcript_id>;<last_exon_isoform_id> - just a unique identifier for each bleedthrough last exon + annotated transcript pair
- gene_name - gene symbol
- peptide_seq_cryptic - full peptide sequence of the last exon inserted into the matching annotated transcript
- peptide_seq_full - full peptide sequence of the full-length matching annotated transcript
- peptide_seq_cryptic_uniq - peptide sequence uniquely encoded by the bleedthrough last exon. In some cases the bleedthrough begins with a stop codon (these are marked with a *)
- longest_match_slice - index (zero-based) of the first position along the sequence where peptide_seq_full & peptide_seq_cryptic no longer completely align. I used this to extract peptide_seq_cryptic_uniq and retained here for convenience/sanity check.
