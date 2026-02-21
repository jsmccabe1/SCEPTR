#!/usr/bin/env python3

import pandas as pd
import sys
import os

def extract_uniprot_ids(blast_file, output_file):
    """
    Extract UniProt IDs from DIAMOND BLAST results.
    
    Parameters:
    blast_file (str): Path to the BLAST results file
    output_file (str): Path to the output mapping file
    """
    
    print(f"Reading BLAST results from {blast_file}")
    
    # Read the BLAST results
    try:
        blast_df = pd.read_csv(blast_file, sep='\t', header=None, 
                              names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                    'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                    'evalue', 'bitscore'])
        print(f"Read {len(blast_df)} BLAST hits")
    except Exception as e:
        print(f"Error reading BLAST results: {e}")
        # Create empty placeholder file and exit
        with open(output_file, 'w') as f:
            f.write("Sequence_ID\tUniProt_ID\tPercent_Identity\tE_value\n")
        return
    
    # Extract sequence IDs and UniProt IDs (assuming UniProt IDs are in the 'sseqid' column)
    mapping_df = blast_df[['qseqid', 'sseqid', 'pident', 'evalue']].copy()
    
    # Clean UniProt IDs - extract accession from format like 'sp|P12345|ABC1_HUMAN'
    def extract_uniprot_id(id_str):
        if '|' in id_str:
            parts = id_str.split('|')
            if len(parts) >= 2:
                return parts[1]
        return id_str
    
    mapping_df['sseqid'] = mapping_df['sseqid'].apply(extract_uniprot_id)
    
    # Rename columns
    mapping_df.columns = ['Sequence_ID', 'UniProt_ID', 'Percent_Identity', 'E_value']
    
    # Save to TSV file
    mapping_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Processed {len(mapping_df)} UniProt mappings")
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <blast_results.tsv> <output_mapping.tsv>")
        sys.exit(1)
    
    blast_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Check if input file exists
    if not os.path.exists(blast_file):
        print(f"Error: Input file {blast_file} does not exist!")
        sys.exit(1)
    
    extract_uniprot_ids(blast_file, output_file)
