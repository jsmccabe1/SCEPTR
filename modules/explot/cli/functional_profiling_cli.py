#!/usr/bin/env python3
"""
CLI for SCEPTR ExPlot functional profiling.

Processes gene annotation/expression data to categorise genes into functional
categories, calculate enrichment statistics, and generate visualisations.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import os
import sys
import json
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

# Resolve module path relative to this script
_script_dir = os.path.dirname(os.path.abspath(__file__))
_module_dir = os.path.abspath(os.path.join(_script_dir, '..'))
sys.path.insert(0, _module_dir)

import go_utils
import categorisation
import enrichment
from visualisation import radar_charts, bar_charts
from reporting import html_report

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.functional_cli')


def parse_arguments():
    parser = argparse.ArgumentParser(description='SCEPTR ExPlot: Functional Profiling')
    parser.add_argument('integrated_results', help='Path to integrated results TSV')
    parser.add_argument('--custom_categories', default=None,
                        help='Path to custom functional categories JSON')
    parser.add_argument('--prefix', default='sceptr', help='Output file prefix')
    parser.add_argument('--multi_category', action='store_true', default=True,
                        help='Allow genes in multiple categories')
    parser.add_argument('--min_percentage', type=float, default=5.0)
    parser.add_argument('--tiers', default='50,100,250,500',
                        help='Expression tiers (comma-separated)')
    parser.add_argument('--p_threshold', type=float, default=0.05)
    parser.add_argument('--correction_method', default='fdr_bh')
    parser.add_argument('--expand_go', action='store_true', default=False,
                        help='Expand keywords with GO hierarchy (off by default)')
    parser.add_argument('--go_expansion_depth', type=int, default=2,
                        help='Max descendant depth for GO expansion (default: 2)')
    parser.add_argument('--outdir', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--debug', action='store_true')
    return parser.parse_args()


def load_default_categories():
    """Load default functional categories from JSON file."""
    possible_paths = [
        os.path.join(_script_dir, '..', 'categories', 'functional_categories.json'),
        os.path.join(_script_dir, '..', 'categories', 'functional',
                     'general_functional_categories.json'),
        '/SCEPTR/modules/explot/categories/functional_categories.json',
    ]
    for path in possible_paths:
        if os.path.exists(path):
            try:
                with open(path) as f:
                    return json.load(f)
            except Exception as e:
                logger.error(f"Error loading {path}: {e}")

    logger.warning("Using hardcoded fallback functional categories")
    return {
        "Central Metabolism": ["glycolysis", "fatty acid", "metabolism", "enzymatic",
                               "oxidative phosphorylation", "dehydrogenase"],
        "Translation & Ribosome Biogenesis": ["ribosomal", "translation", "rRNA processing",
                                                "elongation factor"],
        "Protein Folding & Quality Control": ["chaperone", "protein folding", "ubiquitin",
                                                "proteasome"],
        "Signal Transduction": ["kinase", "phosphorylation", "signaling", "GTPase"],
        "DNA Replication & Repair": ["DNA repair", "DNA replication", "helicase"],
        "Transcriptional Regulation": ["transcription factor", "RNA polymerase", "chromatin"],
        "Cell Division & Proliferation": ["mitosis", "meiosis", "cell cycle", "cell division"],
        "Transport & Membrane": ["vesicle transport", "membrane transport", "ion channel"],
        "RNA Processing": ["RNA processing", "mRNA splicing", "spliceosome"],
        "Cytoskeletal Organization": ["cytoskeleton", "actin", "tubulin", "microtubule"],
        "Uncharacterized": ["hypothetical protein", "unknown function", "uncharacterized"],
    }


def main():
    args = parse_arguments()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info("SCEPTR ExPlot functional profiling")
    logger.info(f"  Input: {args.integrated_results}")
    logger.info(f"  Prefix: {args.prefix}")
    logger.info(f"  Multi-category: {args.multi_category}")
    logger.info(f"  GO expansion: {args.expand_go}")

    # Set up output directory structure
    func_dir = os.path.join(args.outdir, 'functional')
    fig_dir = os.path.join(func_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    # Prefixes for data files vs figure files
    data_prefix = os.path.join(func_dir, args.prefix)
    fig_prefix = os.path.join(fig_dir, args.prefix)

    # Parse tiers
    try:
        tiers = [int(t) for t in args.tiers.split(',')]
    except ValueError:
        logger.error(f"Invalid tiers: {args.tiers}")
        tiers = [50, 100, 250, 500]

    # Load categories
    if args.custom_categories:
        try:
            with open(args.custom_categories) as f:
                functional_categories = json.load(f)
            logger.info(f"Loaded {len(functional_categories)} custom categories")
        except Exception as e:
            logger.error(f"Error loading custom categories: {e}")
            functional_categories = load_default_categories()
    else:
        functional_categories = load_default_categories()
        logger.info(f"Using default categories ({len(functional_categories)} categories)")

    # Detect format and extract keywords + anchor GO IDs + core keywords
    cat_format = categorisation.detect_category_format(functional_categories)
    keyword_map, anchor_map, core_keyword_sets = categorisation.normalize_categories(functional_categories)
    logger.info(f"Category format: {cat_format}")
    core_count = sum(1 for s in core_keyword_sets.values() if len(s) > 0)
    if core_count > 0:
        logger.info(f"Core keywords defined for {core_count}/{len(core_keyword_sets)} categories")

    # Expand anchor GO IDs if present (always-on for v2 format)
    expanded_go_sets = None
    has_anchors = any(len(ids) > 0 for ids in anchor_map.values())
    go_dag = None

    if has_anchors:
        logger.info("Expanding anchor GO IDs via GO hierarchy...")
        try:
            go_dag = go_utils.load_go_dag()
            expanded_go_sets = go_utils.expand_anchor_go_ids(anchor_map, go_dag)
        except Exception as e:
            logger.warning(f"Could not expand anchor GO IDs: {e}")
            logger.warning("Falling back to keyword-only assignment")

    # Optional keyword-based GO expansion (existing behaviour, off by default)
    if args.expand_go:
        logger.info("Expanding keywords with GO terms...")
        if go_dag is None:
            go_dag = go_utils.load_go_dag()
        expanded_categories = go_utils.expand_category_keywords(
            keyword_map, go_dag,
            max_descendant_depth=args.go_expansion_depth)
        logger.info("Keyword GO expansion complete")

        # Save expanded categories for transparency
        try:
            with open(f"{data_prefix}_expanded_categories.json", 'w') as f:
                json.dump(expanded_categories, f, indent=2)
        except Exception as e:
            logger.warning(f"Could not save expanded categories: {e}")
    else:
        expanded_categories = functional_categories
        if not has_anchors:
            logger.info("GO expansion disabled (use --expand_go to enable)")

    # Load expression data
    logger.info(f"Loading expression data...")
    try:
        df = pd.read_csv(args.integrated_results, sep='\t')
    except Exception:
        try:
            df = pd.read_csv(args.integrated_results, sep=',')
        except Exception as e:
            logger.error(f"Cannot read input file: {e}")
            sys.exit(1)

    logger.info(f"Loaded {len(df)} genes")

    if 'TPM' not in df.columns:
        candidates = [c for c in df.columns if c.lower() in ['tpm', 'expression', 'fpkm']]
        if candidates:
            df['TPM'] = df[candidates[0]]
            logger.info(f"Using '{candidates[0]}' as expression column")
        else:
            logger.error("No TPM/expression column found")
            sys.exit(1)

    df_sorted = df.sort_values(by='TPM', ascending=False)

    gene_subsets = {f"top_{t}": df_sorted.head(t) for t in tiers}

    # Categorise all genes
    logger.info("Categorising all genes...")
    all_results = categorisation.categorise_genes(
        df_sorted, expanded_categories, allow_multiple=args.multi_category,
        expanded_go_sets=expanded_go_sets,
        core_keyword_sets=core_keyword_sets)

    # Save assignment method statistics if dual-method was used
    if all_results.get('method_summary'):
        try:
            with open(f"{data_prefix}_BP_MF_assignment_methods.json", 'w') as f:
                json.dump(all_results['method_summary'], f, indent=2)
            logger.info(f"Saved assignment method statistics")
        except Exception as e:
            logger.warning(f"Could not save assignment methods: {e}")

    # Save core specificity statistics if core keywords are defined
    if all_results.get('core_summary'):
        has_any_core = any(
            cs['core'] > 0 for cs in all_results['core_summary'].values())
        if has_any_core:
            try:
                with open(f"{data_prefix}_BP_MF_core_specificity.json", 'w') as f:
                    json.dump(all_results['core_summary'], f, indent=2)
                logger.info(f"Saved core specificity statistics")
            except Exception as e:
                logger.warning(f"Could not save core specificity: {e}")

    # Process each tier
    results = {}
    for tier_name, subset_df in gene_subsets.items():
        logger.info(f"Processing {tier_name}...")
        tier_results = categorisation.categorise_genes(
            subset_df, expanded_categories, allow_multiple=args.multi_category,
            expanded_go_sets=expanded_go_sets,
            core_keyword_sets=core_keyword_sets)

        enrichment_stats = enrichment.calculate_enrichment(
            tier_results['category_counts'],
            len(subset_df), len(df_sorted),
            all_results['category_counts'],
            gene_overlap=tier_results['multi_category_genes'] if args.multi_category else None)

        enrichment_stats = enrichment.adjust_p_values(
            enrichment_stats, method=args.correction_method)

        # Compute core specificity percentage per category for this tier
        tier_core_pct = {}
        tier_core_summary = tier_results.get('core_summary', {})
        for cat, cs in tier_core_summary.items():
            total = cs['core'] + cs['extended']
            tier_core_pct[cat] = round(cs['core'] / total * 100, 1) if total > 0 else 0.0

        results[tier_name] = {
            'counts': tier_results['category_counts'],
            'categorised': tier_results['total_categorised'],
            'uncategorised': tier_results['uncategorised'],
            'enrichment': enrichment_stats,
            'genes_in_category': tier_results['genes_in_category'],
            'multi_category_genes': tier_results['multi_category_genes'],
            'core_pct': tier_core_pct,
            'core_summary': tier_core_summary,
        }

    # Log significant categories
    for tier_name, tier_data in results.items():
        sig = enrichment.get_significant_categories(
            tier_data['enrichment'], p_threshold=args.p_threshold, use_adjusted_p=True)
        logger.info(f"{tier_name}: {len(sig)} significant categories")

    # Generate visualisations
    tier_names = [f"top_{t}" for t in tiers]
    cat_names = list(keyword_map.keys())

    try:
        radar_charts.create_radar_chart(
            results, tier_names, cat_names, fig_prefix, "BP_MF",
            min_percentage=args.min_percentage)

        bar_charts.create_bar_charts(
            results, tier_names, cat_names, fig_prefix, "BP_MF")

        # Create enrichment bar charts for each tier
        for tier_name in tier_names:
            bar_charts.create_enrichment_bar_chart(
                results, tier_name, fig_prefix, "BP_MF")

        # Create grouped bar chart across tiers
        bar_charts.create_grouped_bar_chart(
            results, tier_names, cat_names, fig_prefix, "BP_MF")

        # Create multi-tier enrichment bar chart
        bar_charts.create_multi_tier_enrichment_bar_chart(
            results, tier_names, fig_prefix, "BP_MF")

    except Exception as e:
        logger.error(f"Visualisation error: {e}")
        import traceback
        traceback.print_exc()

    # Save enrichment results TSV
    enrichment.save_enrichment_tsv(results, tier_names, data_prefix, "BP_MF")

    # Generate HTML report
    try:
        html_report.generate_functional_report(
            results, all_results, data_prefix, len(df),
            figures_dir=fig_dir)
    except Exception as e:
        logger.error(f"HTML report error: {e}")
        # Fallback: save JSON
        for tier_data in results.values():
            for cat, s in tier_data['enrichment'].items():
                for k, v in s.items():
                    if isinstance(v, (np.integer, np.floating)):
                        tier_data['enrichment'][cat][k] = float(v)
        with open(f"{data_prefix}_functional_results.json", 'w') as f:
            json.dump(results, f, indent=2, default=str)

    logger.info("Functional profiling complete!")


if __name__ == "__main__":
    main()
