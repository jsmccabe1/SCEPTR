#!/usr/bin/env python3
"""
CLI for SCEPTR ExPlot cellular component profiling.

Processes gene annotation/expression data to categorise genes into cellular
component categories, calculate enrichment, and generate visualisations.

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

_script_dir = os.path.dirname(os.path.abspath(__file__))
_module_dir = os.path.abspath(os.path.join(_script_dir, '..'))
sys.path.insert(0, _module_dir)

import go_utils
import categorisation
import enrichment
import continuous_enrichment
from visualisation import bar_charts, continuous_curves
from reporting import interactive_report

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.cellular_cli')


def parse_arguments():
    parser = argparse.ArgumentParser(description='SCEPTR: Cellular Component Enrichment Profiling')
    parser.add_argument('integrated_results', help='Path to integrated results TSV')
    parser.add_argument('--custom_categories', default=None,
                        help='Path to custom cellular categories JSON')
    parser.add_argument('--prefix', default='sceptr', help='Output file prefix')
    parser.add_argument('--multi_category', action='store_true', default=True)
    parser.add_argument('--min_percentage', type=float, default=5.0)
    parser.add_argument('--tiers', default='50,100,250,500')
    parser.add_argument('--p_threshold', type=float, default=0.05)
    parser.add_argument('--correction_method', default='fdr_bh')
    parser.add_argument('--expand_go', action='store_true', default=False,
                        help='Expand keywords with GO hierarchy (off by default)')
    parser.add_argument('--go_expansion_depth', type=int, default=2)
    parser.add_argument('--outdir', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--continuous', action='store_true', default=True,
                        help='Compute continuous enrichment functions (default: on)')
    parser.add_argument('--no_continuous', action='store_true', default=False,
                        help='Disable continuous enrichment computation')
    parser.add_argument('--continuous_step', type=int, default=5,
                        help='Step size for continuous enrichment (default: 5)')
    parser.add_argument('--continuous_k_min', type=int, default=10,
                        help='Minimum tier size for continuous enrichment (default: 10)')
    parser.add_argument('--continuous_k_max', type=int, default=None,
                        help='Maximum tier size for continuous enrichment (default: N/2)')
    parser.add_argument('--profile_permutations', type=int, default=1000,
                        help='Permutations for global profile test (default: 1000)')
    parser.add_argument('--debug', action='store_true')
    return parser.parse_args()


def load_default_categories():
    possible_paths = [
        os.path.join(_script_dir, '..', 'categories', 'cellular_categories.json'),
        os.path.join(_script_dir, '..', 'categories', 'cellular',
                     'general_cellular_categories.json'),
        '/SCEPTR/modules/explot/categories/cellular_categories.json',
    ]
    for path in possible_paths:
        if os.path.exists(path):
            try:
                with open(path) as f:
                    return json.load(f)
            except Exception as e:
                logger.error(f"Error loading {path}: {e}")

    logger.warning("Using hardcoded fallback cellular categories")
    return {
        "Nuclear Compartments": ["nucleus", "nuclear membrane", "nuclear pore", "chromatin",
                                  "nucleolus"],
        "Mitochondrial Compartments": ["mitochondrion", "mitochondrial matrix",
                                        "mitochondrial membrane"],
        "Cytoskeletal Structures": ["cytoskeleton", "microtubule", "actin filament",
                                     "centrosome"],
        "Membrane Compartments": ["plasma membrane", "cell membrane", "endosome", "lysosome"],
        "Secretory Pathway": ["endoplasmic reticulum", "golgi apparatus", "vesicle"],
        "Cytoplasmic Inclusions": ["cytoplasm", "ribonucleoprotein complex", "proteasome"],
        "Extracellular Compartments": ["extracellular space", "extracellular matrix"],
        "Unlocalized Components": ["unlocalized", "unknown localization",
                                    "hypothetical protein"],
    }


def main():
    args = parse_arguments()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info("SCEPTR cellular component enrichment profiling")
    logger.info(f"  Input: {args.integrated_results}")
    logger.info(f"  Prefix: {args.prefix}")
    logger.info(f"  GO expansion: {args.expand_go}")

    # Set up output directory structure
    cell_dir = os.path.join(args.outdir, 'cellular')
    fig_dir = os.path.join(cell_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    # Prefixes for data files vs figure files
    data_prefix = os.path.join(cell_dir, args.prefix)
    fig_prefix = os.path.join(fig_dir, args.prefix)

    try:
        tiers = [int(t) for t in args.tiers.split(',')]
    except ValueError:
        tiers = [50, 100, 250, 500]

    # Load categories
    if args.custom_categories:
        try:
            with open(args.custom_categories) as f:
                cellular_categories = json.load(f)
            logger.info(f"Loaded {len(cellular_categories)} custom categories")
        except Exception as e:
            logger.error(f"Error loading custom categories: {e}")
            cellular_categories = load_default_categories()
    else:
        cellular_categories = load_default_categories()
        logger.info(f"Using default categories ({len(cellular_categories)} categories)")

    # Detect format and extract keywords + anchor GO IDs + core keywords
    cat_format = categorisation.detect_category_format(cellular_categories)
    keyword_map, anchor_map, core_keyword_sets = categorisation.normalize_categories(cellular_categories)
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
            go_dag = go_utils.load_go_dag('CC')
            expanded_go_sets = go_utils.expand_anchor_go_ids(anchor_map, go_dag)
        except Exception as e:
            logger.warning(f"Could not expand anchor GO IDs: {e}")
            logger.warning("Falling back to keyword-only assignment")

    # Optional keyword-based GO expansion (existing behaviour, off by default)
    if args.expand_go:
        logger.info("Expanding keywords with CC GO terms...")
        if go_dag is None:
            go_dag = go_utils.load_go_dag('CC')
        expanded_categories = go_utils.expand_category_keywords(
            keyword_map, go_dag,
            max_descendant_depth=args.go_expansion_depth)

        try:
            with open(f"{data_prefix}_expanded_cc_categories.json", 'w') as f:
                json.dump(expanded_categories, f, indent=2)
        except Exception:
            pass
    else:
        expanded_categories = cellular_categories
        if not has_anchors:
            logger.info("GO expansion disabled (use --expand_go to enable)")

    # Load data
    try:
        df = pd.read_csv(args.integrated_results, sep='\t')
    except Exception:
        try:
            df = pd.read_csv(args.integrated_results, sep=',')
        except Exception as e:
            logger.error(f"Cannot read input: {e}")
            sys.exit(1)

    logger.info(f"Loaded {len(df)} genes")

    if 'TPM' not in df.columns:
        candidates = [c for c in df.columns if c.lower() in ['tpm', 'expression', 'fpkm']]
        if candidates:
            df['TPM'] = df[candidates[0]]
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
            with open(f"{data_prefix}_CC_assignment_methods.json", 'w') as f:
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
                with open(f"{data_prefix}_CC_core_specificity.json", 'w') as f:
                    json.dump(all_results['core_summary'], f, indent=2)
                logger.info(f"Saved core specificity statistics")
            except Exception as e:
                logger.warning(f"Could not save core specificity: {e}")

    # Process tiers
    results = {}
    for tier_name, subset_df in gene_subsets.items():
        logger.info(f"Processing {tier_name}...")
        tier_results = categorisation.categorise_genes(
            subset_df, expanded_categories, allow_multiple=args.multi_category,
            expanded_go_sets=expanded_go_sets,
            core_keyword_sets=core_keyword_sets)

        enrichment_stats = enrichment.calculate_enrichment(
            tier_results['category_counts'], len(subset_df), len(df_sorted),
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

    for tier_name, tier_data in results.items():
        sig = enrichment.get_significant_categories(
            tier_data['enrichment'], p_threshold=args.p_threshold, use_adjusted_p=True)
        logger.info(f"{tier_name}: {len(sig)} significant categories")

    # Continuous enrichment analysis
    cont_results = None
    if args.continuous and not args.no_continuous:
        logger.info("Computing continuous enrichment functions (CC)...")

        membership, cont_cat_names, cont_bg_counts, _ = \
            continuous_enrichment.build_membership_matrix(
                df_sorted, expanded_categories, categorisation.categorise_genes,
                expanded_go_sets=expanded_go_sets,
                core_keyword_sets=core_keyword_sets)

        N = len(df_sorted)
        k_max = args.continuous_k_max

        k_values, enrichment_matrix = continuous_enrichment.compute_continuous_enrichment(
            membership, cont_cat_names, cont_bg_counts, N,
            k_min=args.continuous_k_min, k_max=k_max, step=args.continuous_step)
        logger.info(f"  Computed enrichment at {len(k_values)} evaluation points")

        _, dkl_values = continuous_enrichment.compute_continuous_dkl(
            membership, cont_cat_names, cont_bg_counts, N,
            k_min=args.continuous_k_min, k_max=k_max, step=args.continuous_step)

        shape_stats = continuous_enrichment.classify_profile_shapes(
            k_values, enrichment_matrix, cont_cat_names)

        logger.info(f"Running global profile test ({args.profile_permutations} permutations)...")
        _, profile_stats = continuous_enrichment.permutation_global_test(
            membership, cont_cat_names, cont_bg_counts, N,
            k_min=args.continuous_k_min, k_max=k_max, step=args.continuous_step,
            n_permutations=args.profile_permutations)

        n_sig = sum(1 for s in profile_stats.values() if s['supremum_p'] < 0.05)
        logger.info(f"  {n_sig}/{len(cont_cat_names)} categories with significant profiles")

        continuous_enrichment.save_continuous_enrichment_tsv(
            k_values, enrichment_matrix, cont_cat_names, data_prefix, "CC")

        try:
            continuous_curves.create_continuous_enrichment_plot(
                k_values, enrichment_matrix, cont_cat_names, profile_stats,
                fig_prefix, "CC", default_tiers=tiers)
            continuous_curves.create_continuous_dkl_plot(
                k_values, dkl_values, fig_prefix, "CC", default_tiers=tiers)
        except Exception as e:
            logger.error(f"Continuous visualisation error: {e}")
            import traceback
            traceback.print_exc()

        cont_results = {
            'k_values': k_values,
            'enrichment_matrix': enrichment_matrix,
            'dkl_values': dkl_values,
            'profile_stats': profile_stats,
            'shape_stats': shape_stats,
            'cat_names': cont_cat_names,
        }

    # Generate static publication figures (trimmed set)
    tier_names = [f"top_{t}" for t in tiers]
    cat_names = list(keyword_map.keys())

    try:
        bar_charts.create_multi_tier_enrichment_bar_chart(
            results, tier_names, fig_prefix, "CC")
    except Exception as e:
        logger.error(f"Visualisation error: {e}")
        import traceback
        traceback.print_exc()

    # Save enrichment results TSV
    enrichment.save_enrichment_tsv(results, tier_names, data_prefix, "CC")

    # Generate interactive HTML report
    try:
        interactive_report.generate_interactive_report(
            results, all_results, data_prefix, len(df),
            cont_results=cont_results,
            chart_type='CC',
            report_title='Cellular Component Profiling',
            description='cellular component localisation',
            figures_dir=fig_dir,
            df_sorted=df_sorted)
    except Exception as e:
        logger.error(f"Interactive report error: {e}")
        import traceback
        traceback.print_exc()

    # Save report data for combined report generation
    try:
        interactive_report.save_report_data(
            results, all_results, cont_results,
            f"{data_prefix}_CC_report_data.json")
    except Exception as e:
        logger.warning(f"Could not save report data: {e}")

    logger.info("Cellular component profiling complete!")


if __name__ == "__main__":
    main()
