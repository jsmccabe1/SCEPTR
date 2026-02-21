#!/usr/bin/env python3

import requests
import pandas as pd
import time
import sys
import os
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

# Create a cache directory for storing downloaded UniProt data
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "uniprot_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

def query_uniprot(uniprot_id):
    """
    Query the UniProt API with caching for information about a specific UniProt ID.
    """
    # Check if we have a cached result
    cache_file = os.path.join(CACHE_DIR, f"{uniprot_id}.json")
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading cache for {uniprot_id}: {e}")
    
    # Make API request if not in cache
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    max_retries = 3
    retry_delay = 1  # seconds
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            if response.status_code == 200:
                result = response.json()
                # Cache the result
                try:
                    with open(cache_file, 'w') as f:
                        json.dump(result, f)
                except Exception as e:
                    print(f"Error caching result for {uniprot_id}: {e}")
                return result
            elif response.status_code == 429:  # Too Many Requests
                wait_time = retry_delay * (2 ** attempt)
                print(f"Rate limit hit, waiting {wait_time} seconds...")
                time.sleep(wait_time)
                continue
            else:
                print(f"Failed to fetch data for UniProt ID: {uniprot_id} (Status: {response.status_code})")
                return None
        except Exception as e:
            print(f"Error querying UniProt for ID {uniprot_id}: {e}")
            time.sleep(retry_delay)
    
    print(f"Failed to fetch data for UniProt ID: {uniprot_id} after {max_retries} attempts")
    return None

def extract_go_terms(data, uniprot_id):
    """
    Extract GO terms from UniProt API response.
    """
    go_terms = {
        'GO_BP': [],  # Biological Process
        'GO_CC': [],  # Cellular Component
        'GO_MF': []   # Molecular Function
    }
    
    try:
        if data and 'uniProtKBCrossReferences' in data:
            for reference in data['uniProtKBCrossReferences']:
                if reference['database'] == 'GO' and 'properties' in reference and len(reference['properties']) > 0:
                    for prop in reference['properties']:
                        if 'value' in prop:
                            go_term = prop['value']
                            if go_term.startswith('P:'):
                                go_terms['GO_BP'].append(go_term)
                            elif go_term.startswith('C:'):
                                go_terms['GO_CC'].append(go_term)
                            elif go_term.startswith('F:'):
                                go_terms['GO_MF'].append(go_term)
    except Exception as e:
        print(f"Error extracting GO terms for {uniprot_id}: {e}")
    
    # Join multiple GO terms into a single string separated by "; "
    go_terms['GO_BP'] = '; '.join(go_terms['GO_BP']) if go_terms['GO_BP'] else None
    go_terms['GO_CC'] = '; '.join(go_terms['GO_CC']) if go_terms['GO_CC'] else None
    go_terms['GO_MF'] = '; '.join(go_terms['GO_MF']) if go_terms['GO_MF'] else None
    
    return go_terms

def process_uniprot_id(index, uniprot_id):
    """Process a single UniProt ID and return its GO terms"""
    if pd.notna(uniprot_id):
        print(f"Processing GO terms for UniProt ID: {uniprot_id}")
        result = query_uniprot(uniprot_id)
        if result:
            return index, uniprot_id, extract_go_terms(result, uniprot_id)
    return index, uniprot_id, {'GO_BP': None, 'GO_CC': None, 'GO_MF': None}

def annotate_go_terms(annotations_file, output_file, max_workers=5):
    """
    Annotate proteins with GO terms from UniProt using parallel processing.
    """
    print(f"Reading protein annotations file: {annotations_file}")
    
    # Load the protein annotations file
    try:
        data = pd.read_csv(annotations_file, sep='\t')
        print(f"Loaded protein annotations with {len(data)} entries")
    except Exception as e:
        print(f"Error reading the protein annotations file: {e}")
        # Create empty output file
        with open(output_file, 'w') as f:
            f.write("Sequence_ID\tUniProt_ID\tPercent_Identity\tE_value\tProtein_Name\tFunction\tDomains\tGO_BP\tGO_CC\tGO_MF\n")
        return
    
    # Prepare columns to store GO annotations
    data['GO_BP'] = None
    data['GO_CC'] = None
    data['GO_MF'] = None
    
    # Process in chunks to avoid memory issues
    chunk_size = 100
    total_rows = len(data)
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for i in range(0, total_rows, chunk_size):
            chunk = data.iloc[i:min(i+chunk_size, total_rows)]
            print(f"Processing chunk {i//chunk_size + 1}/{(total_rows-1)//chunk_size + 1}")
            
            # Submit all jobs in this chunk
            future_to_uniprot = {
                executor.submit(process_uniprot_id, index, row['UniProt_ID']): (index, row['UniProt_ID'])
                for index, row in chunk.iterrows() if pd.notna(row['UniProt_ID'])
            }
            
            # Process completed jobs
            for future in as_completed(future_to_uniprot):
                index, uniprot_id, go_terms = future.result()
                data.at[index, 'GO_BP'] = go_terms['GO_BP']
                data.at[index, 'GO_CC'] = go_terms['GO_CC']
                data.at[index, 'GO_MF'] = go_terms['GO_MF']
            
            # Save intermediate results every chunk
            data.to_csv(output_file + ".in_progress", sep='\t', index=False)
    
    # Save final annotated data
    data.to_csv(output_file, sep='\t', index=False)
    print(f"GO term annotation complete. Final data saved to '{output_file}'")

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <protein_annotations.tsv> <output_annotations.tsv>")
        sys.exit(1)
    
    annotations_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Check if input file exists
    if not os.path.exists(annotations_file):
        print(f"Error: Input file {annotations_file} does not exist!")
        sys.exit(1)
    
    # Check if output file already in progress
    if os.path.exists(output_file + ".in_progress"):
        print(f"Found in-progress output file. Attempting to resume...")
        try:
            # Use the in-progress file as the input
            annotate_go_terms(output_file + ".in_progress", output_file)
        except Exception as e:
            print(f"Failed to resume from in-progress file: {e}")
            annotate_go_terms(annotations_file, output_file)
    else:
        annotate_go_terms(annotations_file, output_file)
