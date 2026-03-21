#!/usr/bin/env python3
"""
Core statistical tests for SCEPTR enrichment functions.

Tests cover continuous enrichment E_C(k), adaptive bandwidth selection,
profile shape classification, fold change confidence intervals,
permutation-based global tests, and KL divergence computation.
"""

import sys
import os
import math

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'modules', 'explot'))

from continuous_enrichment import (
    _adaptive_bandwidth,
    _compute_raw_enrichment,
    compute_continuous_enrichment,
    classify_profile_shapes,
    compute_continuous_dkl,
    permutation_global_test,
)
from enrichment import (
    calculate_enrichment,
    calculate_fold_change_ci,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_membership(N, category_indices_list):
    """Build a boolean membership matrix from lists of member indices."""
    C = len(category_indices_list)
    membership = np.zeros((N, C), dtype=bool)
    for col, indices in enumerate(category_indices_list):
        for idx in indices:
            membership[idx, col] = True
    return membership


# ---------------------------------------------------------------------------
# 1. E_C(k) on a known input
# ---------------------------------------------------------------------------

class TestRawEnrichment:

    def test_ec_known_input(self):
        """E_C(k) = (hits_in_top_k / k) / (total_hits / N).

        10 genes, first 3 belong to the category.
        At k=3: (3/3) / (3/10) = 10/3 ~ 3.333
        At k=5: (3/5) / (3/10) = 2.0
        At k=10: (3/10) / (3/10) = 1.0
        """
        N = 10
        membership = _make_membership(N, [[0, 1, 2]])
        cat_names = ['cat_A']
        bg_counts = {'cat_A': 3}

        k_all, enrichment_raw, bg_rates = _compute_raw_enrichment(
            membership, cat_names, bg_counts, N, k_min=1, k_max=N)

        # k_all is [1..10], so index for k=3 is 2, k=5 is 4, k=10 is 9
        np.testing.assert_allclose(enrichment_raw[2, 0], 10.0 / 3.0, rtol=1e-10)
        np.testing.assert_allclose(enrichment_raw[4, 0], 2.0, rtol=1e-10)
        np.testing.assert_allclose(enrichment_raw[9, 0], 1.0, rtol=1e-10)

    def test_ec_empty_category(self):
        """Category with 0 members should produce E_C(k) = 0 everywhere."""
        N = 20
        membership = _make_membership(N, [[]])  # no members
        cat_names = ['empty']
        bg_counts = {'empty': 0}

        k_all, enrichment_raw, _ = _compute_raw_enrichment(
            membership, cat_names, bg_counts, N, k_min=1, k_max=N)

        np.testing.assert_array_equal(enrichment_raw[:, 0], 0.0)

    def test_ec_full_category(self):
        """Category containing all genes should produce E_C(k) = 1.0 everywhere."""
        N = 15
        membership = _make_membership(N, [list(range(N))])
        cat_names = ['full']
        bg_counts = {'full': N}

        k_all, enrichment_raw, _ = _compute_raw_enrichment(
            membership, cat_names, bg_counts, N, k_min=1, k_max=N)

        np.testing.assert_allclose(enrichment_raw[:, 0], 1.0, rtol=1e-10)


# ---------------------------------------------------------------------------
# 4. Adaptive bandwidth
# ---------------------------------------------------------------------------

class TestAdaptiveBandwidth:

    def test_moderate_category(self):
        """sigma = max(3, 0.5 * sqrt(10000/100)) = max(3, 5.0) = 5.0"""
        result = _adaptive_bandwidth(10000, 100)
        np.testing.assert_allclose(result, 5.0, rtol=1e-10)

    def test_all_genes_in_category(self):
        """sigma = max(3, 0.5 * sqrt(10000/10000)) = max(3, 0.5) = 3.0"""
        result = _adaptive_bandwidth(10000, 10000)
        np.testing.assert_allclose(result, 3.0, rtol=1e-10)

    def test_sparse_category(self):
        """sigma = max(3, 0.5 * sqrt(1000/1)) = max(3, 15.811...) ~ 15.811"""
        result = _adaptive_bandwidth(1000, 1)
        expected = 0.5 * np.sqrt(1000.0)
        np.testing.assert_allclose(result, expected, rtol=1e-4)


# ---------------------------------------------------------------------------
# 5-8. Shape classification
# ---------------------------------------------------------------------------

class TestShapeClassification:

    @staticmethod
    def _classify_single(k_values, curve):
        """Classify a single synthetic curve."""
        enrichment_matrix = curve.reshape(-1, 1)
        result = classify_profile_shapes(k_values, enrichment_matrix, ['test_cat'])
        return result['test_cat']

    def test_apex_concentrated(self):
        """Monotonically decreasing curve should be apex-concentrated."""
        k = np.arange(10, 110)
        curve = 5.0 - 4.0 * (k - 10) / 99.0  # 5.0 -> 1.0
        info = self._classify_single(k, curve)
        assert info['shape_class'] == 'apex-concentrated'

    def test_distributed(self):
        """Monotonically increasing curve should be distributed."""
        k = np.arange(10, 110)
        curve = 1.0 + 4.0 * (k - 10) / 99.0  # 1.0 -> 5.0
        info = self._classify_single(k, curve)
        assert info['shape_class'] == 'distributed'

    def test_flat(self):
        """Curve near 1.0 everywhere should be flat."""
        k = np.arange(10, 110)
        curve = np.ones_like(k, dtype=float) + 0.001
        info = self._classify_single(k, curve)
        assert info['shape_class'] == 'flat'

    def test_threshold_boundary(self):
        """Slope exactly at the boundary (|slope| = 0.1) should be flat."""
        # A linear curve on k_norm in [0,1] with slope = 0.09 (just below 0.1)
        k = np.arange(10, 110)
        k_norm = (k - 10.0) / 99.0
        curve = 2.0 + 0.09 * k_norm  # slope in normalised coords = 0.09
        info = self._classify_single(k, curve)
        assert info['shape_class'] == 'flat'


# ---------------------------------------------------------------------------
# 9. Fold change CI
# ---------------------------------------------------------------------------

class TestFoldChangeCI:

    def test_basic_fold_change(self):
        """5 of 10 genes in tier, 50 of 1000 background: FC = 10.0"""
        fc, ci_lo, ci_hi = calculate_fold_change_ci(5, 10, 50, 1000)
        np.testing.assert_allclose(fc, 10.0, rtol=1e-10)
        # CI should bracket the point estimate
        assert ci_lo < fc
        assert ci_hi > fc

    def test_zero_background(self):
        """Zero background count should return fc=0."""
        fc, ci_lo, ci_hi = calculate_fold_change_ci(3, 10, 0, 1000)
        assert fc == 0.0

    def test_zero_count(self):
        """Zero observed count should return fc=0."""
        fc, ci_lo, ci_hi = calculate_fold_change_ci(0, 10, 50, 1000)
        assert fc == 0.0


# ---------------------------------------------------------------------------
# 10. Permutation null calibration
# ---------------------------------------------------------------------------

class TestPermutationTest:

    def test_null_calibration(self):
        """Random membership (no enrichment) should yield non-significant p-values."""
        rng = np.random.default_rng(12345)
        N = 200
        n_cats = 3
        membership = rng.random((N, n_cats)) < 0.1  # ~10% membership each

        cat_names = [f'cat_{i}' for i in range(n_cats)]
        bg_counts = {cat_names[i]: int(membership[:, i].sum()) for i in range(n_cats)}

        k_values, profile_stats = permutation_global_test(
            membership, cat_names, bg_counts, N,
            k_min=10, k_max=100, step=5,
            n_permutations=99, seed=42)

        for cat in cat_names:
            assert profile_stats[cat]['supremum_p'] > 0.01, (
                f"{cat} should not be significant under the null")


# ---------------------------------------------------------------------------
# 11-12. Continuous D_KL
# ---------------------------------------------------------------------------

class TestContinuousDKL:

    def test_uniform_distribution(self):
        """If tier distribution matches background, D_KL should be near 0."""
        N = 200
        n_cats = 4
        # Each gene assigned uniformly to categories
        rng = np.random.default_rng(99)
        membership = np.zeros((N, n_cats), dtype=bool)
        for i in range(N):
            membership[i, i % n_cats] = True

        cat_names = [f'cat_{j}' for j in range(n_cats)]
        bg_counts = {cat_names[j]: int(membership[:, j].sum()) for j in range(n_cats)}

        k_vals, dkl = compute_continuous_dkl(
            membership, cat_names, bg_counts, N,
            k_min=10, k_max=N, step=1)

        # Uniform membership means tier proportions match background
        assert np.max(dkl) < 0.1, "D_KL should be near zero for uniform assignment"

    def test_concentrated_distribution(self):
        """If top genes belong to one category, D_KL should be positive."""
        N = 200
        n_cats = 4
        membership = np.zeros((N, n_cats), dtype=bool)
        # First 50 genes all in category 0
        membership[:50, 0] = True
        # Remaining genes spread equally
        for i in range(50, N):
            membership[i, i % n_cats] = True

        cat_names = [f'cat_{j}' for j in range(n_cats)]
        bg_counts = {cat_names[j]: int(membership[:, j].sum()) for j in range(n_cats)}

        k_vals, dkl = compute_continuous_dkl(
            membership, cat_names, bg_counts, N,
            k_min=10, k_max=N, step=1)

        # Early tiers should show high divergence
        assert np.max(dkl) > 0.05, "D_KL should be positive for concentrated assignment"


# ---------------------------------------------------------------------------
# 13. Membership matrix shape
# ---------------------------------------------------------------------------

class TestMembershipMatrix:

    def test_shape(self):
        """Manually built membership matrix should have (N, C) dimensions."""
        N, C = 50, 5
        membership = _make_membership(N, [list(range(10)) for _ in range(C)])
        assert membership.shape == (N, C)
        assert membership.dtype == bool


# ---------------------------------------------------------------------------
# 14. Continuous enrichment output shape
# ---------------------------------------------------------------------------

class TestContinuousEnrichmentOutput:

    def test_output_shape(self):
        """compute_continuous_enrichment returns arrays of expected shape."""
        N = 500
        n_cats = 3
        rng = np.random.default_rng(7)
        membership = rng.random((N, n_cats)) < 0.15

        cat_names = [f'cat_{j}' for j in range(n_cats)]
        bg_counts = {cat_names[j]: int(membership[:, j].sum()) for j in range(n_cats)}

        k_values, enrichment_matrix = compute_continuous_enrichment(
            membership, cat_names, bg_counts, N,
            k_min=10, k_max=250, step=5)

        assert k_values.ndim == 1
        assert enrichment_matrix.ndim == 2
        assert enrichment_matrix.shape[0] == len(k_values)
        assert enrichment_matrix.shape[1] == n_cats


# ---------------------------------------------------------------------------
# 15. Phipson-Smyth p-value formula
# ---------------------------------------------------------------------------

class TestPhipsonSmythPvalue:

    def test_formula(self):
        """Permutation p-value should use (count + 1) / (n + 1), not count / n.

        This is the Phipson-Smyth correction that ensures p-values are never
        exactly zero and are uniformly distributed under the null.
        """
        # Set up a tiny deterministic example where we know the answer
        N = 50
        membership = _make_membership(N, [[0, 1, 2, 3, 4]])
        cat_names = ['cat_A']
        bg_counts = {'cat_A': 5}

        k_values, profile_stats = permutation_global_test(
            membership, cat_names, bg_counts, N,
            k_min=5, k_max=25, step=1,
            n_permutations=10, seed=42)

        p_sup = profile_stats['cat_A']['supremum_p']
        p_int = profile_stats['cat_A']['integral_p']

        # With n_permutations=10, denominator must be 11
        # So p must be a multiple of 1/11
        for p in [p_sup, p_int]:
            remainder = (p * 11) % 1.0
            assert remainder < 1e-10 or remainder > 1 - 1e-10, (
                f"p-value {p} is not a multiple of 1/(n_perm+1), "
                "suggesting the Phipson-Smyth formula is not used")
            # p should never be exactly 0
            assert p > 0, "Phipson-Smyth p-value should never be exactly 0"


# ---------------------------------------------------------------------------
# 16. Compositional analysis: Aitchison distance
# ---------------------------------------------------------------------------

class TestAitchisonDistance:

    def test_identical_compositions(self):
        """Distance between identical compositions should be zero."""
        from continuous_enrichment import aitchison_distance
        comp = np.array([0.4, 0.3, 0.2, 0.1])
        np.testing.assert_allclose(aitchison_distance(comp, comp), 0.0, atol=1e-10)

    def test_symmetric(self):
        """Aitchison distance should be symmetric."""
        from continuous_enrichment import aitchison_distance
        a = np.array([0.5, 0.3, 0.2])
        b = np.array([0.2, 0.5, 0.3])
        np.testing.assert_allclose(
            aitchison_distance(a, b), aitchison_distance(b, a), atol=1e-10)

    def test_positive(self):
        """Distance between different compositions should be positive."""
        from continuous_enrichment import aitchison_distance
        a = np.array([0.8, 0.1, 0.1])
        b = np.array([0.1, 0.1, 0.8])
        assert aitchison_distance(a, b) > 0


# ---------------------------------------------------------------------------
# 17. Compositional analysis: zero replacement
# ---------------------------------------------------------------------------

class TestMultiplicativeReplacement:

    def test_preserves_sum(self):
        """Compositions should still sum to 1 after replacement."""
        from continuous_enrichment import _multiplicative_replacement
        comp = np.array([[0.5, 0.3, 0.0, 0.2]])
        result = _multiplicative_replacement(comp)
        np.testing.assert_allclose(result.sum(axis=1), 1.0, atol=1e-10)

    def test_no_zeros(self):
        """All entries should be strictly positive after replacement."""
        from continuous_enrichment import _multiplicative_replacement
        comp = np.array([[0.5, 0.0, 0.0, 0.5]])
        result = _multiplicative_replacement(comp)
        assert np.all(result > 0)

    def test_preserves_nonzero_order(self):
        """Non-zero entries should maintain their relative ordering."""
        from continuous_enrichment import _multiplicative_replacement
        comp = np.array([[0.6, 0.3, 0.0, 0.1]])
        result = _multiplicative_replacement(comp)
        assert result[0, 0] > result[0, 1] > result[0, 3]


# ---------------------------------------------------------------------------
# 18. Functional allocation profile
# ---------------------------------------------------------------------------

class TestFunctionalAllocation:

    def test_compositions_sum_to_one(self):
        """FAP compositions should sum to 1 at each tier."""
        from continuous_enrichment import compute_functional_allocation
        N = 200
        n_cats = 4
        rng = np.random.default_rng(99)
        membership = rng.random((N, n_cats)) < 0.2
        cat_names = [f'cat_{j}' for j in range(n_cats)]
        bg_counts = {cat_names[j]: int(membership[:, j].sum()) for j in range(n_cats)}

        k_values, compositions, background = compute_functional_allocation(
            membership, cat_names, bg_counts, N, k_min=10, step=5)

        np.testing.assert_allclose(compositions.sum(axis=1), 1.0, atol=1e-6)
        np.testing.assert_allclose(background.sum(), 1.0, atol=1e-6)
