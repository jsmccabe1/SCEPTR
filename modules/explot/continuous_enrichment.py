#!/usr/bin/env python3
"""
Continuous enrichment function computation for SCEPTR.

Computes E_C(k) at fine resolution across the full expression gradient,
producing continuous enrichment curves rather than discrete tier summaries.
Includes permutation-based global profile tests and continuous D_KL.

The enrichment function E_C(k) = (|T_k ∩ C| / k) / (|C| / N) is computed
at every k from k_min to k_max using cumulative sums for O(N×C) efficiency.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional

logger = logging.getLogger('sceptr.explot.continuous')


# ---------------------------------------------------------------------------
# Membership matrix construction
# ---------------------------------------------------------------------------

def build_membership_matrix(df_sorted, categories, categorise_fn,
                            expanded_go_sets=None, core_keyword_sets=None):
    """Build per-gene boolean category membership matrix.

    Args:
        df_sorted: DataFrame sorted by TPM descending, index reset
        categories: Category definitions dict
        categorise_fn: categorisation.categorise_genes function
        expanded_go_sets: Expanded GO ID sets (optional)
        core_keyword_sets: Core keyword sets (optional)

    Returns:
        membership: (N, C) boolean array
        cat_names: list of category names (excluding Uncharacterised)
        bg_counts: dict {category: count in full dataset}
    """
    N = len(df_sorted)

    all_results = categorise_fn(
        df_sorted, categories, allow_multiple=True,
        expanded_go_sets=expanded_go_sets,
        core_keyword_sets=core_keyword_sets,
        verbose=False)

    bg_counts = all_results['category_counts']
    gene_categories = all_results.get('gene_categories', {})

    cat_names = sorted([c for c in categories.keys() if c != 'Uncharacterised'])
    cat_to_idx = {c: i for i, c in enumerate(cat_names)}
    n_cats = len(cat_names)

    membership = np.zeros((N, n_cats), dtype=bool)
    for gid, cats in gene_categories.items():
        try:
            row = int(gid)
        except (ValueError, TypeError):
            continue
        if 0 <= row < N:
            for cat in cats:
                idx = cat_to_idx.get(cat)
                if idx is not None:
                    membership[row, idx] = True

    return membership, cat_names, bg_counts, all_results


# ---------------------------------------------------------------------------
# Continuous enrichment computation
# ---------------------------------------------------------------------------

def compute_continuous_enrichment(membership, cat_names, bg_counts, N,
                                   k_min=10, k_max=None, step=5):
    """Compute E_C(k) at fine resolution using cumulative sums.

    Args:
        membership: (N, C) boolean array (genes sorted by TPM descending)
        cat_names: list of category names
        bg_counts: dict {category: count in full dataset}
        N: total number of genes
        k_min: minimum tier size (default: 10)
        k_max: maximum tier size (default: N//2)
        step: step between evaluation points (default: 5)

    Returns:
        k_values: 1D array of tier sizes
        enrichment_matrix: (K, C) array of fold changes
    """
    if k_max is None:
        k_max = N // 2
    k_max = min(k_max, N)

    k_values = np.arange(k_min, k_max + 1, step)

    # Background rates per category
    n_cats = len(cat_names)
    bg_rates = np.array([bg_counts.get(cat, 0) / N for cat in cat_names])

    # Cumulative category counts along the expression gradient
    cumsum = np.cumsum(membership, axis=0)  # (N, C)

    # Counts at each k value
    counts = cumsum[k_values - 1]  # (K, C) — k_values are 1-based sizes, cumsum 0-indexed

    # Enrichment: (observed_rate / expected_rate)
    observed_rates = counts / k_values[:, None]  # (K, C)

    # Avoid division by zero for categories with no background genes
    with np.errstate(divide='ignore', invalid='ignore'):
        enrichment_matrix = np.where(
            bg_rates[None, :] > 0,
            observed_rates / bg_rates[None, :],
            0.0)

    return k_values, enrichment_matrix


# ---------------------------------------------------------------------------
# Continuous D_KL computation
# ---------------------------------------------------------------------------

def compute_continuous_dkl(membership, cat_names, bg_counts, N,
                            k_min=10, k_max=None, step=5,
                            laplace_alpha=1e-6):
    """Compute D_KL(k) at fine resolution.

    D_KL measures how much the category distribution at tier k diverges
    from the background distribution. Uses Laplace smoothing.

    Args:
        membership: (N, C) boolean array
        cat_names: list of category names
        bg_counts: dict {category: count}
        N: total genes
        k_min, k_max, step: resolution parameters
        laplace_alpha: smoothing constant

    Returns:
        k_values: 1D array of tier sizes
        dkl_values: 1D array of D_KL at each k
    """
    if k_max is None:
        k_max = N // 2
    k_max = min(k_max, N)

    k_values = np.arange(k_min, k_max + 1, step)
    n_cats = len(cat_names)

    # Background distribution (Q)
    q = np.array([bg_counts.get(cat, 0) for cat in cat_names], dtype=float)
    q = q + laplace_alpha
    q = q / q.sum()

    # Cumulative counts
    cumsum = np.cumsum(membership, axis=0)
    counts = cumsum[k_values - 1].astype(float)  # (K, C)

    # Tier distributions (P) with smoothing
    counts_smooth = counts + laplace_alpha
    p = counts_smooth / counts_smooth.sum(axis=1, keepdims=True)  # (K, C)

    # D_KL = sum(p * log(p/q))
    dkl_values = np.sum(p * np.log(p / q[None, :]), axis=1)

    return k_values, dkl_values


# ---------------------------------------------------------------------------
# Profile shape classification
# ---------------------------------------------------------------------------

def classify_profile_shapes(k_values, enrichment_matrix, cat_names):
    """Classify enrichment profile shapes using linear trend.

    Returns:
        shape_stats: {category: {slope, r_value, shape_class}}
            shape_class: 'apex-concentrated' (negative slope),
                         'distributed' (positive slope),
                         'flat' (near-zero slope)
    """
    from scipy.stats import linregress

    shape_stats = {}
    # Normalise k to [0, 1] for comparable slopes
    k_norm = (k_values - k_values[0]) / (k_values[-1] - k_values[0])

    for i, cat in enumerate(cat_names):
        curve = enrichment_matrix[:, i]
        if np.all(curve == 0):
            shape_stats[cat] = {
                'slope': 0.0, 'r_value': 0.0,
                'max_enrichment': 0.0, 'shape_class': 'absent'}
            continue

        result = linregress(k_norm, curve)
        slope = result.slope
        max_fc = float(np.max(curve))

        if abs(slope) < 0.1:
            shape_class = 'flat'
        elif slope < 0:
            shape_class = 'apex-concentrated'
        else:
            shape_class = 'distributed'

        shape_stats[cat] = {
            'slope': round(float(slope), 4),
            'r_value': round(float(result.rvalue), 4),
            'max_enrichment': round(max_fc, 4),
            'shape_class': shape_class,
        }

    return shape_stats


# ---------------------------------------------------------------------------
# Permutation-based global profile test
# ---------------------------------------------------------------------------

def permutation_global_test(membership, cat_names, bg_counts, N,
                             k_min=10, k_max=None, step=5,
                             n_permutations=1000, seed=42):
    """Permutation-based global profile test.

    Tests whether each category's enrichment profile, as a whole, differs
    significantly from flat (E_C(k) = 1 for all k). Gene-category
    assignments are shuffled while preserving category sizes.

    Two test statistics per category:
        - supremum: max|E_C(k) - 1| across all k
        - integral: mean|E_C(k) - 1| across all k

    Args:
        membership: (N, C) boolean array
        cat_names: list of category names
        bg_counts: dict {category: count}
        N: total genes
        k_min, k_max, step: resolution parameters
        n_permutations: number of permutations (default: 1000)
        seed: random seed

    Returns:
        profile_stats: {category: {
            supremum_obs, integral_obs,
            supremum_p, integral_p,
            null_envelope_lower, null_envelope_upper,
            observed_curve
        }}
    """
    if k_max is None:
        k_max = N // 2
    k_max = min(k_max, N)

    k_values = np.arange(k_min, k_max + 1, step)
    n_k = len(k_values)
    n_cats = len(cat_names)

    bg_rates = np.array([bg_counts.get(cat, 0) / N for cat in cat_names])

    rng = np.random.default_rng(seed)

    # --- Observed enrichment curves ---
    cumsum_obs = np.cumsum(membership, axis=0)
    counts_obs = cumsum_obs[k_values - 1]  # (K, C)
    obs_rates = counts_obs / k_values[:, None]

    with np.errstate(divide='ignore', invalid='ignore'):
        obs_curves = np.where(bg_rates[None, :] > 0,
                              obs_rates / bg_rates[None, :], 0.0)  # (K, C)

    # Mask: only evaluate at k values where expected count >= 1
    # This prevents zero-count tiers from dominating the supremum
    expected_counts = k_values[:, None] * bg_rates[None, :]  # (K, C)
    valid_mask = expected_counts >= 1.0  # (K, C)

    # Observed statistics (per category)
    obs_deviation = np.abs(obs_curves - 1.0)  # (K, C)
    # Supremum: max deviation where expected >= 1
    obs_supremum = np.zeros(n_cats)
    for i in range(n_cats):
        if np.any(valid_mask[:, i]):
            obs_supremum[i] = np.max(obs_deviation[valid_mask[:, i], i])
    # Integral: mean deviation where expected >= 1
    obs_integral = np.zeros(n_cats)
    for i in range(n_cats):
        if np.any(valid_mask[:, i]):
            obs_integral[i] = np.mean(obs_deviation[valid_mask[:, i], i])

    # --- Permutation null ---
    null_supremum = np.zeros((n_permutations, n_cats))
    null_integral = np.zeros((n_permutations, n_cats))
    # Store null curves for envelope (collect pointwise values)
    null_curves = np.zeros((n_permutations, n_k, n_cats))

    for p in range(n_permutations):
        perm_idx = rng.permutation(N)
        perm_membership = membership[perm_idx]
        perm_cumsum = np.cumsum(perm_membership, axis=0)
        perm_counts = perm_cumsum[k_values - 1]
        perm_rates = perm_counts / k_values[:, None]

        with np.errstate(divide='ignore', invalid='ignore'):
            perm_enrichment = np.where(
                bg_rates[None, :] > 0,
                perm_rates / bg_rates[None, :], 0.0)

        null_curves[p] = perm_enrichment
        perm_deviation = np.abs(perm_enrichment - 1.0)
        # Apply same valid_mask as observed (expected >= 1)
        for i in range(n_cats):
            if np.any(valid_mask[:, i]):
                null_supremum[p, i] = np.max(perm_deviation[valid_mask[:, i], i])
                null_integral[p, i] = np.mean(perm_deviation[valid_mask[:, i], i])

        if (p + 1) % 200 == 0:
            logger.info(f"  Global profile test: {p + 1}/{n_permutations} permutations")

    # --- P-values and envelopes ---
    profile_stats = {}
    for i, cat in enumerate(cat_names):
        sup_p = (np.sum(null_supremum[:, i] >= obs_supremum[i]) + 1) / (n_permutations + 1)
        int_p = (np.sum(null_integral[:, i] >= obs_integral[i]) + 1) / (n_permutations + 1)

        # Pointwise 95% null envelope
        env_lower = np.percentile(null_curves[:, :, i], 2.5, axis=0)
        env_upper = np.percentile(null_curves[:, :, i], 97.5, axis=0)

        profile_stats[cat] = {
            'supremum_obs': round(float(obs_supremum[i]), 4),
            'integral_obs': round(float(obs_integral[i]), 4),
            'supremum_p': float(sup_p),
            'integral_p': float(int_p),
            'null_envelope_lower': env_lower,
            'null_envelope_upper': env_upper,
            'observed_curve': obs_curves[:, i],
        }

    return k_values, profile_stats


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_continuous_enrichment_tsv(k_values, enrichment_matrix, cat_names,
                                    prefix, analysis_type='BP_MF'):
    """Save continuous enrichment matrix as wide-format TSV.

    Columns: k, Category1, Category2, ...
    Each row is one k value with fold changes for all categories.
    """
    out_path = f"{prefix}_{analysis_type}_continuous_enrichment.tsv"

    header = ['k'] + list(cat_names)
    try:
        with open(out_path, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for j, k in enumerate(k_values):
                row = [str(int(k))]
                for i in range(len(cat_names)):
                    row.append(f"{enrichment_matrix[j, i]:.4f}")
                f.write('\t'.join(row) + '\n')
        logger.info(f"Saved continuous enrichment: {out_path}")
    except Exception as e:
        logger.error(f"Could not save continuous enrichment TSV: {e}")

    return out_path


def save_profile_test_tsv(profile_stats, cat_names, prefix,
                           analysis_type='BP_MF'):
    """Save global profile test results as TSV."""
    out_path = f"{prefix}_{analysis_type}_profile_test.tsv"

    header = ['Category', 'Supremum_Obs', 'Integral_Obs',
              'Supremum_P', 'Integral_P']

    try:
        with open(out_path, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for cat in cat_names:
                s = profile_stats[cat]
                row = [
                    cat,
                    f"{s['supremum_obs']:.4f}",
                    f"{s['integral_obs']:.4f}",
                    f"{s['supremum_p']:.6f}",
                    f"{s['integral_p']:.6f}",
                ]
                f.write('\t'.join(row) + '\n')
        logger.info(f"Saved profile test results: {out_path}")
    except Exception as e:
        logger.error(f"Could not save profile test TSV: {e}")

    return out_path


def save_continuous_dkl_tsv(k_values, dkl_values, prefix,
                             analysis_type='BP_MF'):
    """Save continuous D_KL values as TSV."""
    out_path = f"{prefix}_{analysis_type}_continuous_dkl.tsv"

    try:
        with open(out_path, 'w') as f:
            f.write('k\tD_KL\n')
            for k, dkl in zip(k_values, dkl_values):
                f.write(f"{int(k)}\t{dkl:.6f}\n")
        logger.info(f"Saved continuous D_KL: {out_path}")
    except Exception as e:
        logger.error(f"Could not save continuous D_KL TSV: {e}")

    return out_path


def save_shape_classification_tsv(shape_stats, prefix,
                                   analysis_type='BP_MF'):
    """Save profile shape classification results."""
    out_path = f"{prefix}_{analysis_type}_profile_shapes.tsv"

    header = ['Category', 'Slope', 'R_Value', 'Max_Enrichment', 'Shape_Class']

    try:
        with open(out_path, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for cat, s in shape_stats.items():
                row = [
                    cat,
                    f"{s['slope']:.4f}",
                    f"{s['r_value']:.4f}",
                    f"{s['max_enrichment']:.4f}",
                    s['shape_class'],
                ]
                f.write('\t'.join(row) + '\n')
        logger.info(f"Saved shape classification: {out_path}")
    except Exception as e:
        logger.error(f"Could not save shape classification TSV: {e}")

    return out_path
