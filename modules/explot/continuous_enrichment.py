#!/usr/bin/env python3
"""
Continuous enrichment function computation for SCEPTR.

Computes E_C(t) as a smooth, continuous enrichment function over the
normalised expression gradient t in (0, 1]. The discrete enrichment
E_C(k) = (|T_k n C| / k) / (|C| / N) is evaluated at every integer k
using cumulative sums for O(N*C) efficiency, then smoothed via Gaussian
kernel CDF estimation to produce a differentiable function that can be
evaluated at arbitrary resolution.

Includes permutation-based global profile tests and continuous D_KL.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import logging
import numpy as np
from scipy.ndimage import gaussian_filter1d
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
# Kernel smoothing
# ---------------------------------------------------------------------------

def _adaptive_bandwidth(N, n_category_genes):
    """Choose Gaussian kernel bandwidth (in gene-rank units).

    Bandwidth is proportional to the square root of the expected spacing
    between consecutive category members in the ranked gene list. This
    ensures smoothing scales with the natural length scale of the data:
    sparse categories (wide spacing) get broader smoothing while dense
    categories (narrow spacing) get minimal smoothing that preserves
    sharp features like apex-concentrated peaks.

    sigma = max(3, 0.5 * sqrt(N / |C|))

    Returns sigma in units of gene ranks (applied to the step-1 grid).
    """
    spacing = N / max(n_category_genes, 1)
    return max(3.0, 0.5 * np.sqrt(spacing))


def _smooth_enrichment(enrichment_raw, bg_counts, cat_names, N):
    """Apply Gaussian kernel smoothing to raw enrichment curves.

    Smooths each category's enrichment curve using an adaptive bandwidth
    that accounts for category size. This transforms the discrete step
    function E_C(k) into a smooth, differentiable approximation.

    Args:
        enrichment_raw: (K, C) raw enrichment values at step=1
        bg_counts: dict {category: gene count}
        cat_names: list of category names
        N: total genes

    Returns:
        enrichment_smooth: (K, C) smoothed enrichment values
    """
    n_cats = len(cat_names)
    enrichment_smooth = np.empty_like(enrichment_raw)

    for i, cat in enumerate(cat_names):
        n_genes = bg_counts.get(cat, 0)
        sigma = _adaptive_bandwidth(N, n_genes)
        enrichment_smooth[:, i] = gaussian_filter1d(
            enrichment_raw[:, i], sigma=sigma, mode='nearest')

    return enrichment_smooth


# ---------------------------------------------------------------------------
# Continuous enrichment computation
# ---------------------------------------------------------------------------

def _compute_raw_enrichment(membership, cat_names, bg_counts, N,
                            k_min=10, k_max=None):
    """Compute raw E_C(k) at every integer k (step=1).

    This is the discrete enrichment function before smoothing.

    Returns:
        k_all: 1D array of every integer k from k_min to k_max
        enrichment_raw: (K, C) array of raw fold changes
        bg_rates: 1D array of background rates per category
    """
    if k_max is None:
        k_max = N // 2
    k_max = min(k_max, N)

    k_all = np.arange(k_min, k_max + 1)

    bg_rates = np.array([bg_counts.get(cat, 0) / N for cat in cat_names])

    cumsum = np.cumsum(membership, axis=0)  # (N, C)
    counts = cumsum[k_all - 1]  # (K, C)
    observed_rates = counts / k_all[:, None]

    with np.errstate(divide='ignore', invalid='ignore'):
        enrichment_raw = np.where(
            bg_rates[None, :] > 0,
            observed_rates / bg_rates[None, :],
            0.0)

    return k_all, enrichment_raw, bg_rates


def compute_continuous_enrichment(membership, cat_names, bg_counts, N,
                                   k_min=10, k_max=None, step=5):
    """Compute smoothed continuous enrichment function E_C(t).

    Evaluates the discrete enrichment E_C(k) at every integer k, applies
    Gaussian kernel CDF smoothing to produce a differentiable function,
    then resamples at the requested step size for output.

    The smoothing bandwidth adapts to category size: sparse categories
    receive wider smoothing to suppress sampling noise while dense
    categories retain fine-grained signal.

    Args:
        membership: (N, C) boolean array (genes sorted by TPM descending)
        cat_names: list of category names
        bg_counts: dict {category: count in full dataset}
        N: total number of genes
        k_min: minimum gene rank (default: 10)
        k_max: maximum gene rank (default: N//2)
        step: step between output evaluation points (default: 5)

    Returns:
        k_values: 1D array of output gene ranks
        enrichment_matrix: (K, C) array of smoothed fold changes
    """
    # Compute at every integer k
    k_all, enrichment_raw, _ = _compute_raw_enrichment(
        membership, cat_names, bg_counts, N, k_min=k_min, k_max=k_max)

    # Smooth
    enrichment_smooth = _smooth_enrichment(
        enrichment_raw, bg_counts, cat_names, N)

    # Resample at requested step size
    indices = np.arange(0, len(k_all), step)
    k_values = k_all[indices]
    enrichment_matrix = enrichment_smooth[indices]

    return k_values, enrichment_matrix


# ---------------------------------------------------------------------------
# Continuous D_KL computation
# ---------------------------------------------------------------------------

def compute_continuous_dkl(membership, cat_names, bg_counts, N,
                            k_min=10, k_max=None, step=5,
                            laplace_alpha=1e-6):
    """Compute smoothed D_KL(t) at fine resolution.

    D_KL measures how much the category distribution at tier k diverges
    from the background distribution. Computed at every integer k, then
    Gaussian-smoothed and resampled at the requested step size.

    Args:
        membership: (N, C) boolean array
        cat_names: list of category names
        bg_counts: dict {category: count}
        N: total genes
        k_min, k_max, step: resolution parameters
        laplace_alpha: smoothing constant

    Returns:
        k_values: 1D array of gene ranks
        dkl_values: 1D array of smoothed D_KL at each k
    """
    if k_max is None:
        k_max = N // 2
    k_max = min(k_max, N)

    # Compute at every integer k
    k_all = np.arange(k_min, k_max + 1)
    n_cats = len(cat_names)

    # Background distribution (Q)
    q = np.array([bg_counts.get(cat, 0) for cat in cat_names], dtype=float)
    q = q + laplace_alpha
    q = q / q.sum()

    # Cumulative counts
    cumsum = np.cumsum(membership, axis=0)
    counts = cumsum[k_all - 1].astype(float)  # (K, C)

    # Tier distributions (P) with smoothing
    counts_smooth = counts + laplace_alpha
    p = counts_smooth / counts_smooth.sum(axis=1, keepdims=True)  # (K, C)

    # D_KL = sum(p * log(p/q))
    dkl_raw = np.sum(p * np.log(p / q[None, :]), axis=1)

    # Smooth D_KL curve using median category bandwidth
    cat_sizes = np.array([bg_counts.get(cat, 1) for cat in cat_names])
    median_spacing = N / np.median(cat_sizes[cat_sizes > 0])
    sigma = max(3.0, 0.5 * np.sqrt(median_spacing))
    dkl_smooth = gaussian_filter1d(dkl_raw, sigma=sigma, mode='nearest')

    # Resample at requested step
    indices = np.arange(0, len(k_all), step)
    k_values = k_all[indices]
    dkl_values = dkl_smooth[indices]

    return k_values, dkl_values


# ---------------------------------------------------------------------------
# Compositional analysis: Functional Allocation Profile and Apex Distance
# ---------------------------------------------------------------------------

def compute_functional_allocation(membership, cat_names, bg_counts, N,
                                   k_min=10, k_max=None, step=5,
                                   replace_zeros='multiplicative'):
    """Compute the Functional Allocation Profile (FAP).

    At each tier k, computes the proportion of annotated genes belonging
    to each functional category. This compositional vector traces how
    the cell distributes transcriptional resources across competing
    functional programmes as the expression window broadens.

    Unlike independent EC(k) curves, the FAP respects the compositional
    constraint: if 40% of the apex goes to Translation, that is 40% not
    available to other programmes.

    Args:
        membership: (N, C) boolean array (genes sorted by TPM descending)
        cat_names: list of category names
        bg_counts: dict {category: count in full dataset}
        N: total number of genes
        k_min, k_max, step: resolution parameters
        replace_zeros: strategy for zero proportions in log-ratio transforms
            'multiplicative' (default): multiplicative replacement (Martin-Fernandez et al. 2003)

    Returns:
        k_values: 1D array of gene ranks
        compositions: (K, C) array of proportions (rows sum to 1)
        background: 1D array of background proportions
    """
    if k_max is None:
        k_max = N // 2
    k_max = min(k_max, N)

    k_all = np.arange(k_min, k_max + 1)
    n_cats = len(cat_names)

    # Cumulative counts
    cumsum = np.cumsum(membership, axis=0)
    counts = cumsum[k_all - 1].astype(float)  # (K, C)

    # Compositional closure: proportions among annotated genes
    row_sums = counts.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1.0)  # avoid division by zero
    compositions_raw = counts / row_sums  # (K, C)

    # Background composition
    bg = np.array([bg_counts.get(cat, 0) for cat in cat_names], dtype=float)
    bg_sum = bg.sum()
    if bg_sum > 0:
        background = bg / bg_sum
    else:
        background = np.ones(n_cats) / n_cats

    # Zero replacement for log-ratio transforms
    if replace_zeros == 'multiplicative':
        compositions = _multiplicative_replacement(compositions_raw)
        background = _multiplicative_replacement(background.reshape(1, -1))[0]
    else:
        compositions = compositions_raw
        background = background

    # Resample
    indices = np.arange(0, len(k_all), step)
    k_values = k_all[indices]
    compositions = compositions[indices]

    return k_values, compositions, background


def _multiplicative_replacement(X, delta=None):
    """Replace zeros in compositional data (Martin-Fernandez et al. 2003).

    For each row, zeros are replaced with a small value delta and non-zero
    entries are adjusted to maintain the unit-sum constraint.

    Args:
        X: (K, C) array of compositions (rows should sum to 1)
        delta: replacement value (default: 0.65/C^2, following Palarea-Albaladejo and Martin-Fernandez 2015)

    Returns:
        X_replaced: (K, C) array with zeros replaced
    """
    X = np.array(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(1, -1)

    n_components = X.shape[1]
    if delta is None:
        delta = 0.65 / (n_components ** 2)

    X_out = X.copy()
    for i in range(X.shape[0]):
        row = X_out[i]
        zeros = row == 0
        n_zeros = zeros.sum()
        if n_zeros == 0:
            continue
        if n_zeros == n_components:
            X_out[i] = np.ones(n_components) / n_components
            continue
        # Replace zeros with delta
        row[zeros] = delta
        # Scale non-zeros to preserve sum
        non_zero_sum = row[~zeros].sum()
        correction = 1.0 - n_zeros * delta
        row[~zeros] = row[~zeros] * correction / non_zero_sum
        X_out[i] = row

    return X_out


def clr_transform(compositions):
    """Centered log-ratio transform for compositional data.

    Maps compositions from the simplex to real space, enabling standard
    statistical operations (distances, PCA) that respect compositional
    geometry.

    Args:
        compositions: (K, C) array of strictly positive compositions

    Returns:
        clr: (K, C) array of CLR-transformed values
    """
    log_comp = np.log(compositions)
    geometric_mean = log_comp.mean(axis=1, keepdims=True)
    return log_comp - geometric_mean


def aitchison_distance(comp_a, comp_b):
    """Aitchison distance between two compositions.

    The natural distance metric on the simplex. Equivalent to Euclidean
    distance in CLR-transformed space.

    Args:
        comp_a, comp_b: 1D arrays of strictly positive compositions (same length)

    Returns:
        float: Aitchison distance
    """
    log_ratio = np.log(comp_a / comp_b)
    n = len(comp_a)
    return np.sqrt(np.sum((log_ratio - log_ratio.mean()) ** 2))


def compute_compositional_apex_distance(membership, cat_names, bg_counts, N,
                                         apex_k=50, n_permutations=1000,
                                         seed=42):
    """Compute the Compositional Apex Distance (CAD).

    Measures how functionally specialised the expression apex is relative
    to the whole-transcriptome background, using the Aitchison distance
    on the compositional simplex.

    The Aitchison distance properly accounts for the compositional constraint
    (proportions sum to 1), unlike DKL which can be dominated by a single
    category with extreme enrichment.

    A permutation test assesses significance by shuffling gene-category
    assignments and recomputing CAD under the null.

    Args:
        membership: (N, C) boolean array
        cat_names: list of category names
        bg_counts: dict {category: count}
        N: total genes
        apex_k: tier size for the apex (default: 50)
        n_permutations: number of permutations
        seed: random seed

    Returns:
        dict with: cad (float), p_value (float),
                   apex_composition (1D array), background (1D array),
                   null_distribution (1D array of permuted CAD values)
    """
    n_cats = len(cat_names)

    # Observed compositions
    apex_counts = np.sum(membership[:apex_k], axis=0).astype(float)
    bg = np.array([bg_counts.get(cat, 0) for cat in cat_names], dtype=float)

    # Close to compositions (with zero replacement)
    apex_comp = _multiplicative_replacement(
        (apex_counts / max(apex_counts.sum(), 1)).reshape(1, -1))[0]
    bg_comp = _multiplicative_replacement(
        (bg / max(bg.sum(), 1)).reshape(1, -1))[0]

    # Observed CAD
    cad_obs = aitchison_distance(apex_comp, bg_comp)

    # Permutation null
    rng = np.random.default_rng(seed)
    null_cads = np.zeros(n_permutations)

    for p in range(n_permutations):
        perm_idx = rng.permutation(N)
        perm_membership = membership[perm_idx]
        perm_apex = np.sum(perm_membership[:apex_k], axis=0).astype(float)
        perm_comp = _multiplicative_replacement(
            (perm_apex / max(perm_apex.sum(), 1)).reshape(1, -1))[0]
        null_cads[p] = aitchison_distance(perm_comp, bg_comp)

    # P-value (Phipson-Smyth)
    p_value = (np.sum(null_cads >= cad_obs) + 1) / (n_permutations + 1)

    return {
        'cad': float(cad_obs),
        'p_value': float(p_value),
        'z_score': float((cad_obs - null_cads.mean()) / max(null_cads.std(), 1e-10)),
        'apex_composition': apex_comp,
        'background': bg_comp,
        'null_distribution': null_cads,
        'null_mean': float(null_cads.mean()),
        'null_sd': float(null_cads.std()),
    }


def compute_trajectory_distance(compositions_a, compositions_b):
    """Integrated Aitchison distance between two compositional trajectories.

    Measures how different two transcriptomes' functional allocation
    strategies are across the full expression gradient.

    Args:
        compositions_a, compositions_b: (K, C) arrays of compositions
            (must have the same K and C dimensions)

    Returns:
        float: mean Aitchison distance across all tiers
    """
    assert compositions_a.shape == compositions_b.shape
    n_tiers = compositions_a.shape[0]
    distances = np.array([
        aitchison_distance(compositions_a[i], compositions_b[i])
        for i in range(n_tiers)
    ])
    return float(distances.mean())


# ---------------------------------------------------------------------------
# Profile shape classification
# ---------------------------------------------------------------------------

def classify_profile_shapes(k_values, enrichment_matrix, cat_names,
                            slope_threshold=0.075):
    """Classify enrichment profile shapes using linear trend.

    A descriptive convenience for summarising each category's dominant
    expression mode. The continuous slope is always reported alongside
    the discrete label; all statistical claims should use slopes with
    confidence intervals, not the discrete classification.

    Args:
        k_values: 1D array of gene ranks
        enrichment_matrix: (K, C) array of enrichment values
        cat_names: list of category names
        slope_threshold: absolute slope below which a profile is classified
            as flat (default: 0.075). Cross-dataset validation confirmed
            84% consistency at this threshold.

    Returns:
        shape_stats: {category: {slope, r_value, max_enrichment, shape_class}}
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

        if abs(slope) < slope_threshold:
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

def _build_expression_strata(N, n_strata=10):
    """Partition gene ranks into coarse expression strata.

    Returns:
        strata: 1D int array of length N, assigning each gene to a stratum
        strata_ranges: list of (start, end) index pairs per stratum
    """
    strata = np.zeros(N, dtype=int)
    strata_ranges = []
    stratum_size = N // n_strata
    for s in range(n_strata):
        start = s * stratum_size
        end = (s + 1) * stratum_size if s < n_strata - 1 else N
        strata[start:end] = s
        strata_ranges.append((start, end))
    return strata, strata_ranges


def _stratified_permutation(membership_col, strata, strata_ranges, rng):
    """Generate one stratified-permutation null membership vector.

    Preserves the number of category members in each expression stratum,
    but randomises which genes within each stratum are members.
    """
    N = len(membership_col)
    perm = np.zeros(N, dtype=bool)
    for s, (start, end) in enumerate(strata_ranges):
        n_in_stratum = int(membership_col[start:end].sum())
        if n_in_stratum == 0:
            continue
        stratum_size = end - start
        chosen = rng.choice(stratum_size, size=min(n_in_stratum, stratum_size),
                            replace=False)
        perm[start + chosen] = True
    return perm


def _stratified_permutation_batch(membership, strata_ranges, rng):
    """Batch stratified permutation for all categories simultaneously."""
    N, n_cats = membership.shape
    perm = np.zeros_like(membership)
    for start, end in strata_ranges:
        stratum_size = end - start
        for i in range(n_cats):
            n_in = int(membership[start:end, i].sum())
            if n_in > 0:
                chosen = rng.choice(stratum_size, size=min(n_in, stratum_size), replace=False)
                perm[start + chosen, i] = True
    return perm


def permutation_global_test(membership, cat_names, bg_counts, N,
                             k_min=10, k_max=None, step=5,
                             n_permutations=1000, seed=42,
                             n_strata=10):
    """Expression-conditional permutation test with kernel smoothing.

    Tests whether each category's smoothed enrichment profile shows
    fine-scale rank organisation beyond what is expected from its
    coarse expression composition. The null preserves each category's
    member counts per expression decile while randomising positions
    within deciles. Each category is permuted independently.

    This is a conditional test:
        H0: no category-specific rank structure beyond coarse
            expression composition
    rather than the unconditional:
        H0: category members are uniformly random among all genes

    The conditional null mean curve mu_C(k) replaces the flat
    baseline EC(k) = 1. Test statistics measure deviation of the
    observed curve from mu_C(k), not from 1. This means the test
    detects second-order structure (within-stratum concentration)
    rather than first-order bias (tendency toward high expression).

    Two test statistics per category:
        - integral (primary): mean|EC(k) - mu_C(k)| across all valid k
        - supremum (complementary): max|EC(k) - mu_C(k)| across all valid k

    Categories are permuted independently; inter-category overlap is
    not preserved under the null. Significance for overlapping sets
    should not be treated as independent.

    Args:
        membership: (N, C) boolean array (genes sorted by TPM desc)
        cat_names: list of category names
        bg_counts: dict {category: count}
        N: total genes
        k_min, k_max, step: resolution parameters
        n_permutations: number of permutations (default: 1000)
        seed: random seed
        n_strata: number of expression strata for conditional null
            (default: 10 = deciles)

    Returns:
        k_values: 1D array of gene ranks
        profile_stats: {category: {
            supremum_obs, integral_obs,
            supremum_p, integral_p,
            null_envelope_lower, null_envelope_upper,
            null_mean_curve,
            observed_curve
        }}
    """
    if k_max is None:
        k_max = N // 2
    k_max = min(k_max, N)

    # Compute at every integer k, then resample
    k_all = np.arange(k_min, k_max + 1)
    output_indices = np.arange(0, len(k_all), step)
    k_values = k_all[output_indices]
    n_k = len(k_values)
    n_cats = len(cat_names)

    bg_rates = np.array([bg_counts.get(cat, 0) / N for cat in cat_names])

    rng = np.random.default_rng(seed)

    # Build expression strata (genes are already sorted by TPM descending)
    strata, strata_ranges = _build_expression_strata(N, n_strata)

    # --- Observed enrichment curves (smoothed) ---
    cumsum_obs = np.cumsum(membership, axis=0)
    counts_obs = cumsum_obs[k_all - 1]  # (K_all, C)
    obs_rates = counts_obs / k_all[:, None]

    with np.errstate(divide='ignore', invalid='ignore'):
        obs_raw = np.where(bg_rates[None, :] > 0,
                           obs_rates / bg_rates[None, :], 0.0)

    obs_smooth = _smooth_enrichment(obs_raw, bg_counts, cat_names, N)
    obs_curves = obs_smooth[output_indices]  # (K, C)

    # Mask: only evaluate at k values where expected count >= 1
    expected_counts = k_values[:, None] * bg_rates[None, :]  # (K, C)
    valid_mask = expected_counts >= 1.0  # (K, C)

    # --- Stratified permutation null ---
    null_curves = np.zeros((n_permutations, n_k, n_cats))

    for p in range(n_permutations):
        # Build permuted membership matrix (stratified per category)
        perm_membership = _stratified_permutation_batch(membership, strata_ranges, rng)

        perm_cumsum = np.cumsum(perm_membership, axis=0)
        perm_counts = perm_cumsum[k_all - 1]
        perm_rates = perm_counts / k_all[:, None]

        with np.errstate(divide='ignore', invalid='ignore'):
            perm_raw = np.where(
                bg_rates[None, :] > 0,
                perm_rates / bg_rates[None, :], 0.0)

        perm_smooth = _smooth_enrichment(perm_raw, bg_counts, cat_names, N)
        null_curves[p] = perm_smooth[output_indices]

        if (p + 1) % 200 == 0:
            logger.info(f"  Global profile test: {p + 1}/{n_permutations} permutations")
            print(f"    Permutation {p + 1}/{n_permutations}...", flush=True)

    # --- Conditional null mean curve mu_C(k) ---
    null_mean = np.mean(null_curves, axis=0)  # (K, C)

    # --- Observed statistics: deviation from conditional null mean ---
    obs_deviation = np.abs(obs_curves - null_mean)  # (K, C)
    obs_supremum = np.zeros(n_cats)
    obs_integral = np.zeros(n_cats)
    for i in range(n_cats):
        if np.any(valid_mask[:, i]):
            obs_supremum[i] = np.max(obs_deviation[valid_mask[:, i], i])
            obs_integral[i] = np.mean(obs_deviation[valid_mask[:, i], i])

    # --- Null statistics: deviation from conditional null mean ---
    null_supremum = np.zeros((n_permutations, n_cats))
    null_integral = np.zeros((n_permutations, n_cats))
    for p in range(n_permutations):
        perm_deviation = np.abs(null_curves[p] - null_mean)
        for i in range(n_cats):
            if np.any(valid_mask[:, i]):
                null_supremum[p, i] = np.max(perm_deviation[valid_mask[:, i], i])
                null_integral[p, i] = np.mean(perm_deviation[valid_mask[:, i], i])

    # --- P-values and envelopes ---
    profile_stats = {}
    for i, cat in enumerate(cat_names):
        sup_p = (np.sum(null_supremum[:, i] >= obs_supremum[i]) + 1) / (n_permutations + 1)
        int_p = (np.sum(null_integral[:, i] >= obs_integral[i]) + 1) / (n_permutations + 1)

        # Pointwise 95% null envelope (conditional)
        env_lower = np.percentile(null_curves[:, :, i], 2.5, axis=0)
        env_upper = np.percentile(null_curves[:, :, i], 97.5, axis=0)

        profile_stats[cat] = {
            'supremum_obs': round(float(obs_supremum[i]), 4),
            'integral_obs': round(float(obs_integral[i]), 4),
            'supremum_p': float(sup_p),
            'integral_p': float(int_p),
            'null_envelope_lower': env_lower,
            'null_envelope_upper': env_upper,
            'null_mean_curve': null_mean[:, i],
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
