#!/usr/bin/env python3
"""
Statistical enrichment analysis for SCEPTR ExPlot.

Calculates enrichment of functional/cellular categories in expression tiers
using Fisher's exact test with BH FDR correction and exact binomial confidence
intervals for fold changes.

Note on multi-category assignment:
    When genes can belong to multiple categories, the gene counts per category
    may sum to more than the total number of genes. Fisher's exact test assumes
    independent draws, which is violated. We note this in the results and
    recommend interpreting multi-category enrichment as approximate.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import logging
import math
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
from scipy import stats

try:
    from statsmodels.stats.multitest import multipletests
except ImportError:
    logging.warning("statsmodels not found. FDR correction unavailable.")
    def multipletests(pvals, method='fdr_bh', **kwargs):
        return [False] * len(pvals), list(pvals), None, None

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.enrichment')


# ---------------------------------------------------------------------------
# Core enrichment
# ---------------------------------------------------------------------------

def calculate_enrichment(
    category_counts: Dict[str, int],
    subset_size: int,
    total_genes: int,
    all_counts: Dict[str, int],
    gene_overlap: Optional[int] = None,
    method: str = 'fisher_exact',
) -> Dict[str, Dict[str, Any]]:
    """
    Calculate enrichment statistics for each category.

    Args:
        category_counts: Counts in the subset (tier)
        subset_size: Number of genes in the tier
        total_genes: Total genes in the dataset
        all_counts: Counts in the full dataset
        gene_overlap: Multi-category genes (for reporting only)
        method: 'fisher_exact' or 'hypergeometric'

    Returns:
        {category: {fold_change, p_value, significant, ci_lower, ci_upper, ...}}
    """
    enrichment = {}

    for category, count in category_counts.items():
        bg_count = all_counts.get(category, 0)

        if count == 0 or bg_count == 0:
            enrichment[category] = {
                'count': count,
                'expected': 0,
                'fold_change': 0,
                'p_value': 1.0,
                'significant': False,
                'ci_lower': 0,
                'ci_upper': 0,
            }
            continue

        expected = (bg_count / total_genes) * subset_size
        fold_change = count / expected if expected > 0 else float('inf')

        # Exact binomial CI for fold change
        fc, ci_lower, ci_upper = calculate_fold_change_ci(
            count, subset_size, bg_count, total_genes)

        # P-value
        if method == 'fisher_exact':
            p_value = calculate_fisher_exact_p(count, subset_size, bg_count, total_genes)
        elif method == 'hypergeometric':
            p_value = calculate_hypergeometric_p(count, subset_size, bg_count, total_genes)
        else:
            raise ValueError(f"Unknown method: {method}")

        enrichment[category] = {
            'count': count,
            'expected': round(expected, 2),
            'fold_change': round(fold_change, 2),
            'p_value': p_value,
            'significant': p_value < 0.05,
            'ci_lower': round(ci_lower, 2),
            'ci_upper': round(ci_upper, 2),
        }

    # Flag if multi-category mode was used
    if gene_overlap is not None and gene_overlap > 0:
        for cat in enrichment:
            enrichment[cat]['multi_category_note'] = (
                f"Multi-category mode: {gene_overlap} genes in >1 category. "
                f"Enrichment statistics are approximate.")

    return enrichment


# ---------------------------------------------------------------------------
# Statistical tests
# ---------------------------------------------------------------------------

def calculate_fisher_exact_p(count, subset_size, background_count, total_genes):
    """Fisher's exact test for enrichment (one-sided, greater)."""
    # Clamp to prevent negative contingency cells
    a = min(count, subset_size, background_count)
    b = subset_size - a
    c = background_count - a
    d = total_genes - subset_size - c

    # Ensure non-negative
    if d < 0:
        # Adjust: this happens when subset + background > total (multi-category)
        logger.debug(f"Contingency adjustment: count={count}, subset={subset_size}, "
                     f"bg={background_count}, total={total_genes}")
        d = 0
        c = total_genes - subset_size

    table = np.array([[a, b], [c, d]])

    try:
        _, p = stats.fisher_exact(table, alternative='greater')
        return p
    except Exception as e:
        logger.error(f"Fisher's exact test error: {e}")
        return 1.0


def calculate_hypergeometric_p(count, subset_size, background_count, total_genes):
    """Hypergeometric test for enrichment."""
    try:
        return stats.hypergeom.sf(count - 1, total_genes, background_count, subset_size)
    except Exception as e:
        logger.error(f"Hypergeometric test error: {e}")
        return 1.0


def calculate_fold_change_ci(count, subset_size, background_count, total_genes,
                              confidence=0.95):
    """
    Fold change with exact binomial (Clopper-Pearson) confidence interval.

    Returns (fold_change, ci_lower, ci_upper).
    """
    if background_count == 0 or subset_size == 0:
        return 0.0, 0.0, float('inf')

    p_bg = background_count / total_genes
    expected = p_bg * subset_size
    fold_change = count / expected if expected > 0 else 0.0

    if count == 0:
        return 0.0, 0.0, float('inf')

    # Use binom_test CI via scipy
    n_trial = min(subset_size, total_genes)
    if count > n_trial:
        count = n_trial  # clamp for multi-category edge case

    try:
        result = stats.binomtest(count, n_trial, p_bg)
        ci = result.proportion_ci(confidence_level=confidence, method='exact')
        ci_lower = (ci.low * n_trial) / expected if expected > 0 else 0.0
        ci_upper = (ci.high * n_trial) / expected if expected > 0 else float('inf')
        return fold_change, ci_lower, ci_upper
    except Exception:
        # Fallback for older scipy without binomtest
        try:
            ci_result = stats.binom_test(count, n_trial, p_bg)
            # Can't get CI from binom_test, use Wilson approximation
            p_hat = count / n_trial
            z = stats.norm.ppf(1 - (1 - confidence) / 2)
            denom = 1 + z**2 / n_trial
            centre = (p_hat + z**2 / (2 * n_trial)) / denom
            spread = z * math.sqrt(p_hat * (1 - p_hat) / n_trial + z**2 / (4 * n_trial**2)) / denom
            p_lower = max(0, centre - spread)
            p_upper = min(1, centre + spread)
            ci_lower = (p_lower * n_trial) / expected if expected > 0 else 0.0
            ci_upper = (p_upper * n_trial) / expected if expected > 0 else float('inf')
            return fold_change, ci_lower, ci_upper
        except Exception:
            return fold_change, 0.0, float('inf')


# ---------------------------------------------------------------------------
# Multiple testing correction
# ---------------------------------------------------------------------------

def adjust_p_values(enrichment_results: Dict[str, Dict[str, Any]],
                    method: str = 'fdr_bh') -> Dict[str, Dict[str, Any]]:
    """Apply multiple testing correction to enrichment p-values."""
    categories = list(enrichment_results.keys())
    p_values = [enrichment_results[cat]['p_value'] for cat in categories]

    try:
        rejected, adj_pvals, _, _ = multipletests(p_values, method=method)
        for cat, adj_p, rej in zip(categories, adj_pvals, rejected):
            enrichment_results[cat]['adjusted_p_value'] = adj_p
            enrichment_results[cat]['significant'] = bool(rej)
    except Exception as e:
        logger.error(f"Multiple testing correction failed: {e}")
        for cat in categories:
            enrichment_results[cat]['adjusted_p_value'] = enrichment_results[cat]['p_value']

    return enrichment_results


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def get_significant_categories(
    enrichment_results: Dict[str, Dict[str, Any]],
    p_threshold: float = 0.05,
    fold_change_threshold: float = 1.5,
    use_adjusted_p: bool = True,
) -> List[str]:
    """Get categories with significant enrichment."""
    sig = []
    for cat, s in enrichment_results.items():
        p = s.get('adjusted_p_value', s['p_value']) if use_adjusted_p else s['p_value']
        if p < p_threshold and s['fold_change'] >= fold_change_threshold:
            sig.append(cat)
    return sig


def rank_categories_by_enrichment(
    enrichment_results: Dict[str, Dict[str, Any]],
    use_adjusted_p: bool = True,
) -> List[Tuple[str, float, float]]:
    """Rank categories by fold_change × -log10(p-value)."""
    ranked = []
    for cat, s in enrichment_results.items():
        fc = s['fold_change']
        p = s.get('adjusted_p_value', s['p_value']) if use_adjusted_p else s['p_value']
        p = max(p, 1e-300)
        score = fc * min(-math.log10(p), 300)
        ranked.append((cat, fc, p, score))
    ranked.sort(key=lambda x: x[3], reverse=True)
    return [(cat, fc, p) for cat, fc, p, _ in ranked]


def save_enrichment_tsv(
    results: Dict[str, Dict[str, Any]],
    tier_names: List[str],
    prefix: str,
    analysis_type: str,
) -> str:
    """
    Save enrichment results as a TSV file.

    Columns: Tier, Category, Count, Expected, Fold_Change, CI_Lower, CI_Upper,
             P_Value, Adjusted_P_Value, Significant, Core_Pct

    Args:
        results: Per-tier results dict (keys like 'top_50', values with
                 'enrichment' and optionally 'core_pct' sub-dicts)
        tier_names: Ordered tier names
        prefix: Output file prefix
        analysis_type: 'BP_MF' or 'CC'

    Returns:
        Path to the written TSV file.
    """
    out_path = f"{prefix}_{analysis_type}_enrichment_results.tsv"

    # Check if any tier has core_pct data
    has_core = any(
        results.get(tn, {}).get('core_pct', {})
        for tn in tier_names
    )

    header = ['Tier', 'Category', 'Count', 'Expected', 'Fold_Change',
              'CI_Lower', 'CI_Upper', 'P_Value', 'Adjusted_P_Value',
              'Significant']
    if has_core:
        header.append('Core_Pct')

    rows = []
    for tn in tier_names:
        tier_data = results.get(tn, {})
        enrichment_stats = tier_data.get('enrichment', {})
        core_pct = tier_data.get('core_pct', {})

        for cat, s in enrichment_stats.items():
            row = [
                tn,
                cat,
                s.get('count', 0),
                s.get('expected', 0),
                s.get('fold_change', 0),
                s.get('ci_lower', 0),
                s.get('ci_upper', 0),
                s.get('p_value', 1.0),
                s.get('adjusted_p_value', s.get('p_value', 1.0)),
                s.get('significant', False),
            ]
            if has_core:
                row.append(core_pct.get(cat, 0.0))
            rows.append(row)

    try:
        with open(out_path, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for row in rows:
                f.write('\t'.join(str(v) for v in row) + '\n')
        logger.info(f"Saved enrichment TSV: {out_path}")
    except Exception as e:
        logger.error(f"Could not save enrichment TSV: {e}")

    return out_path
