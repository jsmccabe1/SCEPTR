#!/usr/bin/env python3
"""
SCEPTR: Merge Salmon quantification with functional annotations.

Joins TPM values from Salmon quant.sf with the annotation TSV, matching
on sequence ID (with fallback to base transcript ID without .p suffix).
Also deduplicates by UniProt_ID, keeping the highest-TPM representative
per protein to avoid inflating gene counts in downstream enrichment.

Output columns: original annotation columns + TPM
"""

import sys
import pandas as pd
import argparse
import os
import glob
import re


def merge_annotations_expression(annotations_file, quant_dir, output_file):
    print(f"Merging annotations from: {annotations_file}")
    print(f"With expression data from: {quant_dir}")
    
    # ---- Read annotations ----
    annotations_df = pd.read_csv(annotations_file, sep='\t')
    print(f"Loaded {len(annotations_df)} annotations")
    
    # ---- Find quant.sf ----
    if os.path.isfile(quant_dir) and quant_dir.endswith('quant.sf'):
        quant_file = quant_dir
    elif os.path.isdir(quant_dir):
        candidates = glob.glob(os.path.join(quant_dir, "quant.sf"))
        if not candidates:
            candidates = glob.glob(os.path.join(quant_dir, "**", "quant.sf"), recursive=True)
        if not candidates:
            print(f"ERROR: Could not find quant.sf in {quant_dir}")
            sys.exit(1)
        quant_file = candidates[0]
    else:
        # Try current directory
        candidates = glob.glob("**/quant.sf", recursive=True)
        if not candidates:
            print(f"ERROR: Could not find quant.sf")
            sys.exit(1)
        quant_file = candidates[0]
    
    print(f"Using quantification file: {quant_file}")
    expression_df = pd.read_csv(quant_file, sep='\t')
    print(f"Loaded expression data for {len(expression_df)} transcripts")
    
    # ---- Merge: try exact match on sequence_id first ----
    merged_df = pd.merge(
        annotations_df,
        expression_df[['Name', 'TPM']],
        left_on='sequence_id',
        right_on='Name',
        how='left'
    )
    
    exact_matched = merged_df['TPM'].notna().sum()
    print(f"After exact match: {exact_matched} sequences matched")
    
    # ---- Fallback: strip .p suffix for unmatched ----
    unmatched_mask = merged_df['TPM'].isna()
    n_unmatched = unmatched_mask.sum()
    
    if n_unmatched > 0:
        print(f"Attempting base-ID matching for {n_unmatched} unmatched sequences...")
        
        # Build base-ID lookup from expression data
        expr_base = expression_df[['Name', 'TPM']].copy()
        expr_base['base_id'] = expr_base['Name'].apply(
            lambda x: re.sub(r'\.p\d+$', '', str(x))
        )
        
        # For unmatched annotations, compute base ID and merge by key
        unmatched_idx = merged_df.index[unmatched_mask]
        base_ids = merged_df.loc[unmatched_idx, 'sequence_id'].apply(
            lambda x: re.sub(r'\.p\d+$', '', str(x))
        )
        
        # Merge on base ID (by key, not position)
        base_lookup = expr_base.drop_duplicates(subset='base_id', keep='first').set_index('base_id')['TPM']
        merged_df.loc[unmatched_idx, 'TPM'] = base_ids.map(base_lookup).values
        
        total_matched = merged_df['TPM'].notna().sum()
        print(f"After base-ID matching: {total_matched} total sequences matched")
    
    # Drop helper column
    if 'Name' in merged_df.columns:
        merged_df = merged_df.drop('Name', axis=1)
    
    # Fill missing TPM with 0
    merged_df['TPM'] = merged_df['TPM'].fillna(0)
    
    # Move TPM to second column
    cols = merged_df.columns.tolist()
    if 'TPM' in cols:
        cols.remove('TPM')
        cols.insert(1, 'TPM')
        merged_df = merged_df[cols]
    
    # ---- Deduplicate by UniProt_ID: keep highest TPM per protein ----
    # This is critical for downstream enrichment - multiple TransDecoder
    # ORFs can map to the same UniProt entry, inflating gene counts
    if 'uniprot_id' in merged_df.columns:
        before_dedup = len(merged_df)
        # Sort by TPM descending so first occurrence is highest
        merged_df = merged_df.sort_values('TPM', ascending=False)
        merged_df = merged_df.drop_duplicates(subset='uniprot_id', keep='first')
        after_dedup = len(merged_df)
        if before_dedup > after_dedup:
            print(f"Deduplicated: {before_dedup} → {after_dedup} "
                  f"(removed {before_dedup - after_dedup} duplicate UniProt mappings, "
                  f"kept highest TPM per protein)")
        # Re-sort by TPM for tier analysis
        merged_df = merged_df.sort_values('TPM', ascending=False).reset_index(drop=True)
    
    merged_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nMerge completed:")
    print(f"  Total sequences: {len(merged_df)}")
    print(f"  Sequences with TPM > 0: {(merged_df['TPM'] > 0).sum()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge annotations with expression data')
    parser.add_argument('--annotations', required=True, help='Annotations file')
    parser.add_argument('--quant', required=True, help='Quantification directory or file')
    parser.add_argument('--output', required=True, help='Output merged file')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.annotations):
        print(f"ERROR: Annotations file not found: {args.annotations}")
        sys.exit(1)
    
    merge_annotations_expression(args.annotations, args.quant, args.output)
