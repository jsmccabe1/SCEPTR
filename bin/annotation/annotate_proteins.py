#!/usr/bin/env python3

import requests
import pandas as pd
import time
import sys
import os

def query_uniprot(uniprot_id):
    """
    Query the UniProt API for information about a specific UniProt ID.
    
    Parameters:
    uniprot_id (str): The UniProt ID to query
    
    Returns:
    dict: The JSON response from the UniProt API, or None if the query failed
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    max_retries = 3
    retry_delay = 1  # seconds
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            if response.status_code == 200:
                return response.json()
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

def extract_annotations(data):
    """
    Extract protein name, function, and domain information from UniProt API response.
    
    Parameters:
    data (dict): The JSON response from the UniProt API
    
    Returns:
    dict: Dictionary containing the extracted annotations
    """
    annotations = {
        'Protein_Name': None,
        'Function': None,
        'Domains': None
    }
    
    try:
        if data and 'proteinDescription' in data and 'recommendedName' in data['proteinDescription']:
            annotations['Protein_Name'] = data['proteinDescription']['recommendedName']['fullName']['value']
        
        if data and 'comments' in data:
            for comment in data['comments']:
                if comment['commentType'] == 'FUNCTION' and 'texts' in comment and len(comment['texts']) > 0:
                    annotations['Function'] = comment['texts'][0]['value']
                    break
        
        if data and 'features' in data:
            domains = []
            for feature in data['features']:
                if feature['type'] == 'DOMAIN' and 'description' in feature:
                    domains.append(feature['description'])
            annotations['Domains'] = '; '.join(domains) if domains else None
    except Exception as e:
        print(f"Error extracting annotations: {e}")
    
    return annotations

def annotate_proteins(mapping_file, output_file):
    """
    Annotate proteins with information from UniProt.
    
    Parameters:
    mapping_file (str): Path to the TSV file containing Sequence_ID to UniProt_ID mappings
    output_file (str): Path to the output TSV file for saving the annotated data
    """
    print(f"Reading mapping file: {mapping_file}")
    
    # Load the UniProt mapping file
    try:
        data = pd.read_csv(mapping_file, sep='\t')
        print(f"Loaded mapping file with {len(data)} entries")
    except Exception as e:
        print(f"Error reading the mapping file: {e}")
        # Create empty output file
        with open(output_file, 'w') as f:
            f.write("Sequence_ID\tUniProt_ID\tPercent_Identity\tE_value\tProtein_Name\tFunction\tDomains\n")
        return
    
    # Prepare columns for annotations
    data['Protein_Name'] = None
    data['Function'] = None
    data['Domains'] = None
    
    # Process in chunks to avoid overloading the API
    chunk_size = 100
    total_rows = len(data)
    
    for i in range(0, total_rows, chunk_size):
        chunk = data.iloc[i:min(i+chunk_size, total_rows)]
        print(f"Processing chunk {i//chunk_size + 1}/{(total_rows-1)//chunk_size + 1}")
        
        for index, row in chunk.iterrows():
            uniprot_id = row['UniProt_ID']
            if pd.notna(uniprot_id):
                print(f"Processing UniProt ID: {uniprot_id}")
                result = query_uniprot(uniprot_id)
                if result:
                    annotations = extract_annotations(result)
                    data.at[index, 'Protein_Name'] = annotations['Protein_Name']
                    data.at[index, 'Function'] = annotations['Function']
                    data.at[index, 'Domains'] = annotations['Domains']
                # Add a small delay to avoid overwhelming the API
                time.sleep(0.1)
    
    # Save annotated data
    data.to_csv(output_file, sep='\t', index=False)
    print(f"Protein annotation complete. Data saved to '{output_file}'")

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <uniprot_mapping.tsv> <output_annotations.tsv>")
        sys.exit(1)
    
    mapping_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Check if input file exists
    if not os.path.exists(mapping_file):
        print(f"Error: Input file {mapping_file} does not exist!")
        sys.exit(1)
    
    annotate_proteins(mapping_file, output_file)
