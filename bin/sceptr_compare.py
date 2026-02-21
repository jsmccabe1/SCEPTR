#!/usr/bin/env python3
"""
SCEPTR cross-sample comparison module.

Compares functional enrichment profiles between two conditions using a
gene-label permutation test, and computes concordance metrics. Produces
a differential enrichment TSV, concordance TSV, publication-quality figures,
and a self-contained HTML comparison report.

This is the production implementation of the permutation comparison
described in the SCEPTR manuscript (Section 3.5, Supplementary Table S5).

Usage:
    python bin/sceptr_compare.py \
        --condition_a results_mock/integrated_data/integrated_annotations_expression.tsv \
        --condition_b results_infected/integrated_data/integrated_annotations_expression.tsv \
        --label_a "Mock" --label_b "Infected" \
        --category_set vertebrate_host \
        --tiers 50,100,250,500 \
        --n_permutations 10000 \
        --output_dir results_comparison/ \
        --output_prefix sceptr_comparison

Author: James McCabe
Module: SCEPTR
"""

import os
import sys
import json
import argparse
import logging
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

# Resolve module paths
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..'))
_explot_dir = os.path.join(_project_root, 'modules', 'explot')
sys.path.insert(0, _explot_dir)

import categorisation
import enrichment as enrichment_mod
from visualisation import comparison_charts
from reporting import comparison_report

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


DEFAULT_TIERS = [50, 100, 250, 500]
DEFAULT_N_PERMS = 10000


# ---------------------------------------------------------------------------
# Category loading
# ---------------------------------------------------------------------------

def load_categories(category_set: str, custom_path: str = None) -> dict:
    """Load functional category definitions."""
    if custom_path:
        with open(custom_path) as f:
            return json.load(f)

    cat_dir = os.path.join(_explot_dir, 'categories', 'functional')
    path = os.path.join(cat_dir, f'{category_set}_functional_categories.json')
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)

    # Fallback: general
    fallback = os.path.join(cat_dir, 'general_functional_categories.json')
    if os.path.exists(fallback):
        logger.warning(f"Category set '{category_set}' not found, using general")
        with open(fallback) as f:
            return json.load(f)

    raise FileNotFoundError(
        f"No category file found for '{category_set}' in {cat_dir}")


def try_expand_anchors(categories: dict) -> dict:
    """Expand anchor GO IDs through GO hierarchy if available."""
    _, anchor_map, _ = categorisation.normalize_categories(categories)
    has_anchors = any(len(ids) > 0 for ids in anchor_map.values())
    if not has_anchors:
        return None
    try:
        import go_utils
        go_dag = go_utils.load_go_dag()
        return go_utils.expand_anchor_go_ids(anchor_map, go_dag)
    except Exception as e:
        logger.warning(f"Could not expand GO IDs: {e}")
        return None


# ---------------------------------------------------------------------------
# Pre-computed category assignment matrix (from permutation_comparison.py)
# ---------------------------------------------------------------------------

def precompute_category_assignments(df, categories, expanded_go_sets):
    """Pre-compute per-gene category membership as a boolean matrix.

    Returns:
        cat_matrix: np.ndarray (n_genes, n_categories), dtype bool
        cat_names: list of category names (excluding Uncharacterised)
        bg_counts: dict of {category: count} in full dataset
    """
    res = categorisation.categorise_genes(
        df, categories, allow_multiple=True,
        expanded_go_sets=expanded_go_sets, verbose=False)

    cat_names = sorted([c for c in categories.keys() if c != 'Uncharacterised'])
    n_genes = len(df)
    cat_matrix = np.zeros((n_genes, len(cat_names)), dtype=bool)

    gene_categories = res.get('gene_categories', {})
    cat_name_to_idx = {c: i for i, c in enumerate(cat_names)}

    for gid, cats in gene_categories.items():
        try:
            row_idx = int(gid)
        except (ValueError, TypeError):
            continue
        if 0 <= row_idx < n_genes:
            for cat in cats:
                idx = cat_name_to_idx.get(cat)
                if idx is not None:
                    cat_matrix[row_idx, idx] = True

    return cat_matrix, cat_names, res['category_counts']


def compute_fold_changes_from_matrix(cat_matrix, tier_indices, n_total, cat_names):
    """Compute fold changes from pre-computed category matrix."""
    n_tier = len(tier_indices)
    bg_counts = cat_matrix.sum(axis=0)
    tier_counts = cat_matrix[tier_indices].sum(axis=0)

    fc_dict = {}
    for i, cat in enumerate(cat_names):
        if bg_counts[i] == 0:
            fc_dict[cat] = 0.0
        else:
            expected = (bg_counts[i] / n_total) * n_tier
            fc_dict[cat] = tier_counts[i] / expected if expected > 0 else 0.0
    return fc_dict


# ---------------------------------------------------------------------------
# Per-sample tiered enrichment (with CIs and significance)
# ---------------------------------------------------------------------------

def compute_sample_enrichment(df, categories, expanded_go_sets, tiers):
    """Compute full enrichment statistics for a single sample.

    Returns:
        enrichment_results: {tier: {category: {fold_change, p_value, ...}}}
        all_results: categorisation results for full dataset
    """
    _, _, core_keyword_sets = categorisation.normalize_categories(categories)

    df_sorted = df.sort_values('TPM', ascending=False).reset_index(drop=True)

    all_results = categorisation.categorise_genes(
        df_sorted, categories, allow_multiple=True,
        expanded_go_sets=expanded_go_sets,
        core_keyword_sets=core_keyword_sets,
        verbose=False)

    results = {}
    for tier in tiers:
        tier_name = f"top_{tier}"
        subset_df = df_sorted.head(tier)
        tier_results = categorisation.categorise_genes(
            subset_df, categories, allow_multiple=True,
            expanded_go_sets=expanded_go_sets,
            core_keyword_sets=core_keyword_sets,
            verbose=False)

        enrichment_stats = enrichment_mod.calculate_enrichment(
            tier_results['category_counts'],
            len(subset_df), len(df_sorted),
            all_results['category_counts'],
            gene_overlap=tier_results['multi_category_genes'])

        enrichment_stats = enrichment_mod.adjust_p_values(enrichment_stats)

        # Core specificity
        tier_core_pct = {}
        for cat, cs in tier_results.get('core_summary', {}).items():
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
        }

    return results, all_results


# ---------------------------------------------------------------------------
# Permutation test
# ---------------------------------------------------------------------------

def run_permutation_test(df_a, df_b, categories, expanded_go_sets,
                         tiers, n_permutations, seed=42):
    """Run gene-label permutation test comparing two conditions.

    Returns:
        results_df: DataFrame with per-category, per-tier results
        cat_names: list of category names
    """
    df_a = df_a.sort_values('TPM', ascending=False).reset_index(drop=True)
    df_b = df_b.sort_values('TPM', ascending=False).reset_index(drop=True)

    n_a = len(df_a)
    n_b = len(df_b)

    logger.info(f"Condition A: {n_a} genes")
    logger.info(f"Condition B: {n_b} genes")

    # Pre-compute category assignments
    logger.info("Pre-computing category assignments...")
    cat_matrix_a, cat_names, _ = precompute_category_assignments(
        df_a, categories, expanded_go_sets)
    cat_matrix_b, _, _ = precompute_category_assignments(
        df_b, categories, expanded_go_sets)

    results = []
    rng = np.random.default_rng(seed)

    for tier in tiers:
        logger.info(f"Tier {tier}: Running {n_permutations} permutations...")

        # Observed fold changes
        tier_idx_a = np.arange(min(tier, n_a))
        tier_idx_b = np.arange(min(tier, n_b))

        obs_fc_a = compute_fold_changes_from_matrix(
            cat_matrix_a, tier_idx_a, n_a, cat_names)
        obs_fc_b = compute_fold_changes_from_matrix(
            cat_matrix_b, tier_idx_b, n_b, cat_names)

        obs_diff = {cat: obs_fc_b[cat] - obs_fc_a[cat] for cat in cat_names}

        # Pool category matrices
        pooled_matrix = np.vstack([cat_matrix_a, cat_matrix_b])
        n_pooled = n_a + n_b
        tpm_pooled = np.concatenate([df_a['TPM'].values, df_b['TPM'].values])

        # Permutation loop
        perm_diffs = {cat: np.zeros(n_permutations) for cat in cat_names}

        for p in range(n_permutations):
            perm_idx = rng.permutation(n_pooled)
            group_a_idx = perm_idx[:n_a]
            group_b_idx = perm_idx[n_a:]

            tpm_a = tpm_pooled[group_a_idx]
            tpm_b = tpm_pooled[group_b_idx]

            tier_a = group_a_idx[np.argsort(tpm_a)[::-1][:min(tier, len(group_a_idx))]]
            tier_b = group_b_idx[np.argsort(tpm_b)[::-1][:min(tier, len(group_b_idx))]]

            fc_a = compute_fold_changes_from_matrix(
                pooled_matrix, tier_a, n_a, cat_names)
            fc_b = compute_fold_changes_from_matrix(
                pooled_matrix, tier_b, n_b, cat_names)

            for cat in cat_names:
                perm_diffs[cat][p] = fc_b[cat] - fc_a[cat]

            if (p + 1) % 2000 == 0:
                logger.info(f"  Permutation {p + 1}/{n_permutations}")

        # Compute p-values (two-sided)
        for cat in cat_names:
            obs_abs = abs(obs_diff[cat])
            perm_abs = np.abs(perm_diffs[cat])
            p_value = (np.sum(perm_abs >= obs_abs) + 1) / (n_permutations + 1)

            results.append({
                'Category': cat,
                'Tier': tier,
                'FC_A': round(obs_fc_a[cat], 4),
                'FC_B': round(obs_fc_b[cat], 4),
                'FC_Diff': round(obs_diff[cat], 4),
                'Perm_P_Value': round(p_value, 6),
            })

    results_df = pd.DataFrame(results)

    # Apply per-tier BH correction
    try:
        from statsmodels.stats.multitest import multipletests

        results_df['Perm_FDR'] = 1.0
        for tier in tiers:
            mask = results_df['Tier'] == tier
            pvals = results_df.loc[mask, 'Perm_P_Value'].values
            rejected, adj_pvals, _, _ = multipletests(pvals, method='fdr_bh')
            results_df.loc[mask, 'Perm_FDR'] = adj_pvals
    except Exception as e:
        logger.warning(f"BH correction failed: {e}")
        results_df['Perm_FDR'] = results_df['Perm_P_Value']

    # Add effect size and direction columns
    results_df['Significant'] = results_df['Perm_FDR'] < 0.05
    results_df['Direction'] = results_df['FC_Diff'].apply(
        lambda x: 'B' if x > 0.1 else ('A' if x < -0.1 else 'NS'))

    return results_df, cat_names


# ---------------------------------------------------------------------------
# Concordance metrics
# ---------------------------------------------------------------------------

def spearman_ci(rho, n, alpha=0.05):
    """Fisher z-transform confidence interval for Spearman rho."""
    if n < 4:
        return rho, -1.0, 1.0
    z = np.arctanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    from scipy.stats import norm
    z_crit = norm.ppf(1 - alpha / 2)
    ci_lo = np.tanh(z - z_crit * se)
    ci_hi = np.tanh(z + z_crit * se)
    return rho, float(ci_lo), float(ci_hi)


def compute_concordance(results_a, results_b, df_a, df_b, tiers, cat_names):
    """Compute concordance metrics between two samples.

    Returns:
        concordance: list of dicts with per-tier metrics
    """
    from scipy.stats import spearmanr

    concordance = []

    for tier in tiers:
        tier_name = f"top_{tier}"

        # FC profile correlation
        fc_a = []
        fc_b = []
        for cat in cat_names:
            ea = results_a.get(tier_name, {}).get('enrichment', {}).get(cat, {})
            eb = results_b.get(tier_name, {}).get('enrichment', {}).get(cat, {})
            fc_a.append(ea.get('fold_change', 0))
            fc_b.append(eb.get('fold_change', 0))

        if len(fc_a) >= 3:
            rho, p = spearmanr(fc_a, fc_b)
            rho_val, ci_lo, ci_hi = spearman_ci(rho, len(fc_a))
        else:
            rho_val, ci_lo, ci_hi = 0.0, -1.0, 1.0

        # Jaccard similarity of top-k gene sets
        df_a_sorted = df_a.sort_values('TPM', ascending=False)
        df_b_sorted = df_b.sort_values('TPM', ascending=False)

        # Use gene identifiers for overlap
        id_col = None
        for col in ['gene_names', 'uniprot_id', 'protein_name', 'sequence_id']:
            if col in df_a.columns and col in df_b.columns:
                id_col = col
                break

        if id_col:
            top_a = set(df_a_sorted.head(tier)[id_col].dropna().astype(str))
            top_b = set(df_b_sorted.head(tier)[id_col].dropna().astype(str))
            intersection = top_a & top_b
            union = top_a | top_b
            jaccard = len(intersection) / len(union) if len(union) > 0 else 0.0
            n_shared = len(intersection)
        else:
            jaccard = 0.0
            n_shared = 0

        concordance.append({
            'Tier': tier,
            'Spearman_Rho': round(rho_val, 4),
            'Spearman_CI_Lower': round(ci_lo, 4),
            'Spearman_CI_Upper': round(ci_hi, 4),
            'Jaccard_Similarity': round(jaccard, 4),
            'N_Shared_Genes': n_shared,
        })

    return concordance


# ---------------------------------------------------------------------------
# CI columns for differential enrichment
# ---------------------------------------------------------------------------

def add_ci_columns(diff_df, results_a, results_b, tiers):
    """Add confidence interval columns from per-sample enrichment."""
    ci_lower_a = []
    ci_upper_a = []
    ci_lower_b = []
    ci_upper_b = []

    for _, row in diff_df.iterrows():
        cat = row['Category']
        tier_name = f"top_{int(row['Tier'])}"

        ea = results_a.get(tier_name, {}).get('enrichment', {}).get(cat, {})
        eb = results_b.get(tier_name, {}).get('enrichment', {}).get(cat, {})

        ci_lower_a.append(ea.get('ci_lower', 0))
        ci_upper_a.append(ea.get('ci_upper', 0))
        ci_lower_b.append(eb.get('ci_lower', 0))
        ci_upper_b.append(eb.get('ci_upper', 0))

    diff_df['CI_Lower_A'] = ci_lower_a
    diff_df['CI_Upper_A'] = ci_upper_a
    diff_df['CI_Lower_B'] = ci_lower_b
    diff_df['CI_Upper_B'] = ci_upper_b

    return diff_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='SCEPTR cross-sample comparison',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python bin/sceptr_compare.py \\
        --condition_a results_mock/integrated_data/integrated_annotations_expression.tsv \\
        --condition_b results_infected/integrated_data/integrated_annotations_expression.tsv \\
        --label_a "Mock" --label_b "Infected" \\
        --category_set vertebrate_host \\
        --output_dir results_comparison/
""")

    parser.add_argument('--condition_a', required=True,
                        help='Path to condition A integrated_annotations_expression.tsv')
    parser.add_argument('--condition_b', required=True,
                        help='Path to condition B integrated_annotations_expression.tsv')
    parser.add_argument('--label_a', default='Condition_A',
                        help='Label for condition A (default: Condition_A)')
    parser.add_argument('--label_b', default='Condition_B',
                        help='Label for condition B (default: Condition_B)')
    parser.add_argument('--category_set', default='general',
                        help='Category set name (default: general)')
    parser.add_argument('--functional_categories', default=None,
                        help='Path to custom functional categories JSON')
    parser.add_argument('--cellular_categories', default=None,
                        help='Path to custom cellular categories JSON (reserved for future use)')
    parser.add_argument('--tiers', default='50,100,250,500',
                        help='Comma-separated tier sizes (default: 50,100,250,500)')
    parser.add_argument('--n_permutations', type=int, default=DEFAULT_N_PERMS,
                        help=f'Number of permutations (default: {DEFAULT_N_PERMS})')
    parser.add_argument('--output_dir', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--output_prefix', default='sceptr_comparison',
                        help='Prefix for output files')
    parser.add_argument('--expand_go', action='store_true', default=False,
                        help='Expand keywords with GO hierarchy')
    parser.add_argument('--go_expansion_depth', type=int, default=2,
                        help='Max GO expansion depth (default: 2)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')

    args = parser.parse_args()

    # Parse tiers
    tiers = [int(t.strip()) for t in args.tiers.split(',')]

    # Set up output directories
    os.makedirs(args.output_dir, exist_ok=True)
    fig_dir = os.path.join(args.output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    data_prefix = os.path.join(args.output_dir, args.output_prefix)
    fig_prefix = os.path.join(fig_dir, args.output_prefix)

    # Load categories
    logger.info(f"Loading categories: {args.category_set}")
    categories = load_categories(args.category_set, args.functional_categories)
    logger.info(f"  {len(categories)} categories loaded")

    # Expand GO anchors
    expanded_go_sets = try_expand_anchors(categories)

    # Load datasets
    logger.info(f"Loading condition A: {args.condition_a}")
    df_a = pd.read_csv(args.condition_a, sep='\t')
    logger.info(f"  {len(df_a)} genes")

    logger.info(f"Loading condition B: {args.condition_b}")
    df_b = pd.read_csv(args.condition_b, sep='\t')
    logger.info(f"  {len(df_b)} genes")

    # Validate TPM column
    for label, df in [('A', df_a), ('B', df_b)]:
        if 'TPM' not in df.columns:
            candidates = [c for c in df.columns if c.lower() in ['tpm', 'expression', 'fpkm']]
            if candidates:
                df['TPM'] = df[candidates[0]]
            else:
                logger.error(f"No TPM column found in condition {label}")
                sys.exit(1)

    # Shared gene universe
    id_col = None
    for col in ['uniprot_id', 'gene_names', 'protein_name', 'sequence_id']:
        if col in df_a.columns and col in df_b.columns:
            id_col = col
            break

    if id_col:
        ids_a = set(df_a[id_col].dropna().astype(str))
        ids_b = set(df_b[id_col].dropna().astype(str))
        n_shared = len(ids_a & ids_b)
        logger.info(f"Shared genes ({id_col}): {n_shared}")
        if n_shared < 100:
            logger.warning(f"Only {n_shared} shared genes - comparison may be unreliable")
    else:
        n_shared = min(len(df_a), len(df_b))
        logger.warning("No common ID column found - cannot compute gene overlap")

    # -----------------------------------------------------------------------
    # Step 1: Per-sample enrichment
    # -----------------------------------------------------------------------
    logger.info("Computing per-sample enrichment for condition A...")
    results_a, all_results_a = compute_sample_enrichment(
        df_a, categories, expanded_go_sets, tiers)

    logger.info("Computing per-sample enrichment for condition B...")
    results_b, all_results_b = compute_sample_enrichment(
        df_b, categories, expanded_go_sets, tiers)

    # -----------------------------------------------------------------------
    # Step 2: Permutation test
    # -----------------------------------------------------------------------
    logger.info(f"Running permutation test ({args.n_permutations:,} permutations)...")
    diff_df, cat_names = run_permutation_test(
        df_a, df_b, categories, expanded_go_sets,
        tiers=tiers, n_permutations=args.n_permutations, seed=args.seed)

    # Add CI columns
    diff_df = add_ci_columns(diff_df, results_a, results_b, tiers)

    # -----------------------------------------------------------------------
    # Step 3: Concordance metrics
    # -----------------------------------------------------------------------
    logger.info("Computing concordance metrics...")
    concordance = compute_concordance(
        results_a, results_b, df_a, df_b, tiers, cat_names)

    # -----------------------------------------------------------------------
    # Step 4: Save TSVs
    # -----------------------------------------------------------------------
    diff_tsv = f"{data_prefix}_differential_enrichment.tsv"
    diff_df.to_csv(diff_tsv, sep='\t', index=False)
    logger.info(f"Saved differential enrichment: {diff_tsv}")

    conc_df = pd.DataFrame(concordance)
    conc_tsv = f"{data_prefix}_concordance.tsv"
    conc_df.to_csv(conc_tsv, sep='\t', index=False)
    logger.info(f"Saved concordance metrics: {conc_tsv}")

    # -----------------------------------------------------------------------
    # Step 5: Generate figures
    # -----------------------------------------------------------------------
    logger.info("Generating figures...")

    try:
        comparison_charts.create_radar_overlay(
            results_a, results_b, cat_names,
            args.label_a, args.label_b, fig_prefix)
    except Exception as e:
        logger.warning(f"Radar overlay failed: {e}")

    try:
        comparison_charts.create_differential_heatmap(
            diff_df, cat_names, tiers,
            args.label_a, args.label_b, fig_prefix)
    except Exception as e:
        logger.warning(f"Differential heatmap failed: {e}")

    try:
        comparison_charts.create_fc_barplot(
            diff_df, cat_names, tiers,
            args.label_a, args.label_b, fig_prefix)
    except Exception as e:
        logger.warning(f"FC barplot failed: {e}")

    try:
        comparison_charts.create_volcano_plot(
            diff_df, args.label_a, args.label_b, fig_prefix)
    except Exception as e:
        logger.warning(f"Volcano plot failed: {e}")

    try:
        comparison_charts.create_gradient_overlay(
            diff_df, cat_names, tiers,
            args.label_a, args.label_b, fig_prefix)
    except Exception as e:
        logger.warning(f"Gradient overlay failed: {e}")

    # -----------------------------------------------------------------------
    # Step 6: Generate HTML report
    # -----------------------------------------------------------------------
    logger.info("Generating HTML comparison report...")
    try:
        comparison_report.generate_comparison_report(
            diff_df=diff_df,
            concordance_data=concordance,
            label_a=args.label_a,
            label_b=args.label_b,
            n_genes_a=len(df_a),
            n_genes_b=len(df_b),
            n_shared=n_shared,
            n_permutations=args.n_permutations,
            tiers=tiers,
            output_prefix=data_prefix,
            figures_dir=fig_dir,
        )
    except Exception as e:
        logger.error(f"HTML report generation failed: {e}")
        import traceback
        traceback.print_exc()

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print(f"SCEPTR Cross-Sample Comparison: {args.label_a} vs {args.label_b}")
    print(f"{'=' * 70}")
    print(f"Condition A: {len(df_a):,} genes ({args.label_a})")
    print(f"Condition B: {len(df_b):,} genes ({args.label_b})")
    print(f"Shared genes: {n_shared:,}")
    print(f"Permutations: {args.n_permutations:,}")
    print()

    for tier in tiers:
        tier_data = diff_df[diff_df['Tier'] == tier]
        n_sig = tier_data['Significant'].sum()
        print(f"  Tier {tier}: {n_sig} significant categories (BH FDR < 0.05)")

    print()
    print("Concordance:")
    for entry in concordance:
        print(f"  Tier {entry['Tier']}: Spearman rho = {entry['Spearman_Rho']:.3f}, "
              f"Jaccard = {entry['Jaccard_Similarity']:.3f}")

    print(f"\nOutputs:")
    print(f"  Differential enrichment: {diff_tsv}")
    print(f"  Concordance metrics: {conc_tsv}")
    print(f"  HTML report: {data_prefix}_comparison_report.html")
    print(f"  Figures: {fig_dir}/")
    print(f"{'=' * 70}\n")


if __name__ == '__main__':
    main()
