#!/usr/bin/env python3
"""
CLI for generating a unified SCEPTR interactive report.

Takes functional + cellular report data JSONs and produces one HTML
report with tabbed navigation across all analyses.

Usage:
    python generate_report_cli.py integrated_results.tsv \
        --functional functional/sceptr_BP_MF_report_data.json \
        --cellular cellular/sceptr_CC_report_data.json \
        --output sceptr_report.html

Author: James McCabe
Module: SCEPTR ExPlot
"""

import os
import sys
import argparse
import logging
import pandas as pd

_script_dir = os.path.dirname(os.path.abspath(__file__))
_module_dir = os.path.abspath(os.path.join(_script_dir, '..'))
sys.path.insert(0, _module_dir)

from reporting import interactive_report

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.report_cli')


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='SCEPTR: Generate Unified Interactive Report')
    parser.add_argument('integrated_results',
                        help='Path to integrated results TSV')
    parser.add_argument('--functional', required=True,
                        help='Path to functional report data JSON')
    parser.add_argument('--cellular', default=None,
                        help='Path to cellular report data JSON (optional)')
    parser.add_argument('--output', default='sceptr_report.html',
                        help='Output HTML path')
    return parser.parse_args()


def main():
    args = parse_arguments()

    logger.info("Loading expression data for gene names...")
    try:
        df = pd.read_csv(args.integrated_results, sep='\t')
    except Exception:
        df = pd.read_csv(args.integrated_results, sep=',')

    if 'TPM' not in df.columns:
        candidates = [c for c in df.columns if c.lower() in ['tpm', 'expression', 'fpkm']]
        if candidates:
            df['TPM'] = df[candidates[0]]
    df_sorted = df.sort_values(by='TPM', ascending=False)
    total_genes = len(df)

    logger.info(f"Loading functional data: {args.functional}")
    func_data = interactive_report.load_report_data(args.functional)

    if args.cellular and os.path.exists(args.cellular):
        logger.info(f"Loading cellular data: {args.cellular}")
        cell_data = interactive_report.load_report_data(args.cellular)

        interactive_report.generate_combined_report(
            func_data, cell_data, args.output, total_genes,
            df_sorted=df_sorted)
    else:
        logger.info("No cellular data - generating functional-only report")
        interactive_report.generate_interactive_report(
            func_data['results'], func_data['all_results'],
            os.path.splitext(args.output)[0], total_genes,
            cont_results=func_data.get('cont_results'),
            chart_type='BP_MF',
            report_title='Functional Profiling',
            description='biological process and molecular function',
            df_sorted=df_sorted)

    logger.info("Report generation complete!")


if __name__ == "__main__":
    main()
