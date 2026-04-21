#!/usr/bin/env python3
"""
Augment SCEPTR's integrated_annotations_expression.tsv with InterProScan
(Pfam) hits and the Pfam->GO mappings InterProScan emits via InterPro2GO.

This is the merge step that turns InterProScan output into something
SCEPTR's downstream categorisation can consume.

Inputs
------
  --annotations  integrated_annotations_expression.tsv produced by
                 MergeExpressionAnnotations (TPM + UniProt-derived
                 annotations).
  --interproscan InterProScan TSV (15 columns, the standard format).
  --output       output path; the augmented annotation TSV.

Behaviour
---------
  - For each protein already in the input annotation, append InterProScan
    Pfam/InterPro descriptions to its protein_name field, and append any
    InterProScan-derived GO IDs to its GO_Biological_Process field. The
    GO IDs are written in the same `term name [GO:0001234]` format that
    SCEPTR's downstream parser already expects, so no other code needs
    to change.
  - For proteins that have InterProScan hits but were missing from the
    input annotation (i.e. they did not produce a UniProt BLAST hit),
    add new rows. TPM is filled in from the same quant.sf if --quant is
    provided; otherwise TPM is left as NaN.
  - The output column schema is identical to the input. Downstream code
    cannot tell that augmentation happened, except that more rows have
    annotations.

Usage
-----
  python3 augment_annotations.py \\
      --annotations integrated_annotations_expression.tsv \\
      --interproscan dtrenchii.fasta.tsv \\
      --output integrated_annotations_expression.tsv

The output path may be the same as the input (in-place augmentation).
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# Standard InterProScan TSV column names (15 columns).
IPRSCAN_COLUMNS = [
    'protein', 'md5', 'length', 'analysis', 'sig_acc', 'sig_desc',
    'start', 'stop', 'evalue', 'status', 'date',
    'ipr_acc', 'ipr_desc', 'go_terms', 'pathways',
]

# These are the SCEPTR annotation columns we may need to write into.
ANNOT_COLUMNS = [
    'sequence_id', 'TPM', 'uniprot_id', 'protein_name', 'gene_names',
    'organism', 'length',
    'GO_Biological_Process', 'GO_Cellular_Component', 'GO_Molecular_Function',
]


def load_iprscan(path):
    """Read an InterProScan TSV and return {protein_id: {pfam_ids,
    ipr_ids, go_ids, descriptions}}.

    An empty or header-only TSV returns an empty dict. This happens
    when InterProScan runs successfully but finds no Pfam hits, a
    legitimate outcome for highly divergent proteomes that should not
    crash the pipeline.
    """
    if Path(path).stat().st_size == 0:
        return {}

    try:
        df = pd.read_csv(path, sep='\t', names=IPRSCAN_COLUMNS, header=None,
                         na_values=['-'], dtype=str)
    except pd.errors.EmptyDataError:
        return {}

    if df.empty:
        return {}

    by_protein = defaultdict(lambda: {
        'pfam_ids': set(), 'ipr_ids': set(),
        'go_ids': set(), 'descriptions': set(),
    })
    for _, row in df.iterrows():
        pid = row['protein']
        if pd.notna(row['sig_acc']):
            by_protein[pid]['pfam_ids'].add(str(row['sig_acc']))
        if pd.notna(row['sig_desc']):
            by_protein[pid]['descriptions'].add(str(row['sig_desc']).strip())
        if pd.notna(row['ipr_acc']):
            by_protein[pid]['ipr_ids'].add(str(row['ipr_acc']))
        if pd.notna(row['ipr_desc']):
            by_protein[pid]['descriptions'].add(str(row['ipr_desc']).strip())
        if pd.notna(row['go_terms']):
            for go_id in re.findall(r'GO:\d{7}', str(row['go_terms'])):
                by_protein[pid]['go_ids'].add(go_id)
    return dict(by_protein)


def format_iprscan_go_text(go_ids):
    """Format InterProScan GO IDs into the same form SCEPTR's parser
    expects in GO_Biological_Process columns: 'name [GO:nnnnnnn]'.
    Since we don't have GO term names without loading the OBO, we use
    a placeholder name; SCEPTR's categorise_genes only matches GO IDs
    via regex, so the name is informational only.
    """
    return '; '.join(f'InterProScan term [{g}]' for g in sorted(go_ids))


def append_field(existing, addition, sep='; '):
    """Append addition to existing field, handling NaN/empty cases."""
    s = '' if (existing is None or pd.isna(existing)) else str(existing).strip()
    a = '' if not addition else str(addition).strip()
    if not s:
        return a
    if not a:
        return s
    if a in s:
        return s
    return s + sep + a


def augment(annotations_df, iprscan):
    """Return a new DataFrame containing augmented existing rows + new
    rows for InterProScan-only proteins.
    """
    df = annotations_df.copy()

    for col in ANNOT_COLUMNS:
        if col not in df.columns:
            df[col] = '' if col != 'TPM' else np.nan

    existing_ids = set(df['sequence_id'].astype(str))

    n_augmented = 0
    n_added = 0

    # Augment existing rows in place
    for idx, row in df.iterrows():
        pid = str(row['sequence_id'])
        if pid not in iprscan:
            continue
        n_augmented += 1
        info = iprscan[pid]
        if info['descriptions']:
            extra_desc = '; '.join(sorted(info['descriptions']))
            df.at[idx, 'protein_name'] = append_field(
                row.get('protein_name', ''), extra_desc, sep=' | ')
        if info['go_ids']:
            extra_go = format_iprscan_go_text(info['go_ids'])
            df.at[idx, 'GO_Biological_Process'] = append_field(
                row.get('GO_Biological_Process', ''), extra_go, sep='; ')

    # Add new rows for proteins not in the existing annotation but
    # found by InterProScan.
    new_rows = []
    for pid, info in iprscan.items():
        if pid in existing_ids:
            continue
        descriptions = '; '.join(sorted(info['descriptions']))
        go_text = format_iprscan_go_text(info['go_ids'])
        new_rows.append({
            'sequence_id': pid,
            'TPM': np.nan,
            'uniprot_id': '',
            'protein_name': descriptions,
            'gene_names': '',
            'organism': 'InterProScan-only',
            'length': '',
            'GO_Biological_Process': go_text,
            'GO_Cellular_Component': '',
            'GO_Molecular_Function': '',
        })
        n_added += 1

    if new_rows:
        df_new = pd.DataFrame(new_rows)
        # Ensure column order matches existing df
        df_new = df_new.reindex(columns=df.columns)
        df = pd.concat([df, df_new], ignore_index=True)

    return df, n_augmented, n_added


def fill_tpm_for_new_rows(df, quant_path):
    """Optional: look up TPM for InterProScan-only rows from a Salmon
    quant.sf so they participate in the expression ranking.
    """
    if not quant_path or not Path(quant_path).exists():
        return df
    quant = pd.read_csv(quant_path, sep='\t')
    tpm_map = dict(zip(quant['Name'].astype(str), quant['TPM']))

    n_filled = 0
    for idx, row in df.iterrows():
        if pd.isna(row.get('TPM')):
            pid = str(row['sequence_id'])
            # Try exact match
            if pid in tpm_map:
                df.at[idx, 'TPM'] = tpm_map[pid]
                n_filled += 1
                continue
            # Try base transcript ID (strip TransDecoder .pN suffix)
            base = re.sub(r'\.p\d+$', '', pid)
            if base in tpm_map:
                df.at[idx, 'TPM'] = tpm_map[base]
                n_filled += 1
    print(f"  Filled TPM for {n_filled} InterProScan-only rows")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--annotations', required=True,
                    help='integrated_annotations_expression.tsv from '
                         'MergeExpressionAnnotations')
    ap.add_argument('--interproscan', required=True,
                    help='InterProScan TSV output (15-column standard format)')
    ap.add_argument('--output', required=True,
                    help='Output path. May be the same as --annotations.')
    ap.add_argument('--quant', default=None,
                    help='Optional Salmon quant.sf to fill TPM for '
                         'InterProScan-only proteins')
    args = ap.parse_args()

    print(f"Loading annotation: {args.annotations}")
    annotations_df = pd.read_csv(args.annotations, sep='\t')
    print(f"  {len(annotations_df)} rows")

    print(f"Loading InterProScan TSV: {args.interproscan}")
    if not Path(args.interproscan).exists():
        print(f"ERROR: InterProScan TSV not found: {args.interproscan}",
              file=sys.stderr)
        sys.exit(1)
    iprscan = load_iprscan(args.interproscan)
    print(f"  {len(iprscan)} proteins with Pfam hits")

    print("Augmenting...")
    augmented, n_aug, n_add = augment(annotations_df, iprscan)
    print(f"  augmented {n_aug} existing rows")
    print(f"  added {n_add} new rows from InterProScan-only proteins")
    print(f"  total rows: {len(augmented)} (was {len(annotations_df)})")

    if args.quant:
        augmented = fill_tpm_for_new_rows(augmented, args.quant)

    print(f"Writing: {args.output}")
    augmented.to_csv(args.output, sep='\t', index=False)
    print("Done.")


if __name__ == '__main__':
    main()
