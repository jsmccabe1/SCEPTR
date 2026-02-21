#!/usr/bin/env python3
"""
SCEPTR: Prepare GO enrichment data from annotation TSV.

Parses the annotation file's GO term columns into a long-format CSV
suitable for topGO enrichment analysis. Column names are standardised
to match what tpm_go_enrichment_analysis.R expects.

Output columns: sequence_id, UniProt_ID, GO_ID, GO_Term, Category
"""

import sys
import pandas as pd
import argparse


def parse_go_column(value):
    """Parse a semicolon-delimited GO term string into (go_id, go_term) pairs.
    
    Expected format: "term description [GO:XXXXXXX]; another term [GO:YYYYYYY]"
    """
    if pd.isna(value) or not str(value).strip() or str(value).strip() == 'nan':
        return []
    
    results = []
    for entry in str(value).split(';'):
        entry = entry.strip()
        if not entry or entry == 'nan':
            continue
        
        if '[GO:' in entry:
            parts = entry.split('[GO:')
            if len(parts) == 2:
                go_term = parts[0].strip()
                go_id = 'GO:' + parts[1].rstrip(']').strip()
                results.append((go_id, go_term))
            else:
                results.append(('Unknown', entry))
        else:
            results.append(('Unknown', entry))
    
    return results


def prepare_go_enrichment(input_file, output_file, analysis_file):
    print(f"Preparing GO enrichment data from: {input_file}")
    
    df = pd.read_csv(input_file, sep='\t')
    print(f"Loaded {len(df)} annotations")
    print(f"Columns: {df.columns.tolist()}")
    
    if len(df) == 0:
        print("Warning: No annotation data found")
        pd.DataFrame(columns=['sequence_id', 'UniProt_ID', 'GO_ID', 'GO_Term', 'Category']).to_csv(output_file, index=False)
        pd.DataFrame(columns=['Category', 'Total_Terms', 'Sequences_With_Terms']).to_csv(analysis_file, index=False)
        return
    
    # Map annotation columns to GO categories
    go_column_map = {
        'GO_Biological_Process': 'BP',
        'GO_Cellular_Component': 'CC',
        'GO_Molecular_Function': 'MF',
    }
    
    go_data = []
    
    for _, row in df.iterrows():
        seq_id = row.get('sequence_id', '')
        uniprot_id = row.get('uniprot_id', '')
        
        for col_name, category in go_column_map.items():
            if col_name not in df.columns:
                continue
            
            for go_id, go_term in parse_go_column(row.get(col_name)):
                go_data.append({
                    'sequence_id': seq_id,
                    'UniProt_ID': uniprot_id,
                    'GO_ID': go_id,
                    'GO_Term': go_term,
                    'Category': category,
                })
    
    if not go_data:
        print("No GO terms found in annotations")
        pd.DataFrame(columns=['sequence_id', 'UniProt_ID', 'GO_ID', 'GO_Term', 'Category']).to_csv(output_file, index=False)
        pd.DataFrame(columns=['Category', 'Total_Terms', 'Sequences_With_Terms']).to_csv(analysis_file, index=False)
        return
    
    go_df = pd.DataFrame(go_data)
    
    # Remove entries with unknown GO IDs
    known_go = go_df[go_df['GO_ID'] != 'Unknown']
    unknown_count = len(go_df) - len(known_go)
    if unknown_count > 0:
        print(f"Removed {unknown_count} entries with unknown GO IDs")
    go_df = known_go
    
    go_df.to_csv(output_file, index=False)
    
    # Analysis summary
    analysis_rows = []
    for category in ['BP', 'CC', 'MF']:
        cat_data = go_df[go_df['Category'] == category]
        analysis_rows.append({
            'Category': category,
            'Total_Terms': cat_data['GO_ID'].nunique(),
            'Sequences_With_Terms': cat_data['sequence_id'].nunique(),
            'Total_Associations': len(cat_data),
        })
    
    pd.DataFrame(analysis_rows).to_csv(analysis_file, index=False)
    
    print(f"GO enrichment data prepared:")
    print(f"  Total GO term associations: {len(go_df)}")
    print(f"  Unique sequences with GO terms: {go_df['sequence_id'].nunique()}")
    for cat in ['BP', 'CC', 'MF']:
        n = len(go_df[go_df['Category'] == cat])
        print(f"  {cat} associations: {n}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prepare GO enrichment data')
    parser.add_argument('input_file', help='Input annotations file')
    parser.add_argument('--output', required=True, help='Output GO terms file')
    parser.add_argument('--analysis', required=True, help='Output analysis file')
    
    args = parser.parse_args()
    prepare_go_enrichment(args.input_file, args.output, args.analysis)
