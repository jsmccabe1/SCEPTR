#!/usr/bin/env python3
"""
Continuous enrichment curve visualisation for SCEPTR.

Produces publication-quality plots of continuous enrichment functions
with permutation null envelopes.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import logging
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

logger = logging.getLogger('sceptr.explot.continuous_vis')

# Publication colour palette (colourblind-safe)
CATEGORY_COLOURS = [
    '#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F',
    '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85',
    '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF', '#AEC7E8',
    '#FFBB78', '#98DF8A', '#FF9896', '#C5B0D5', '#C49C94',
]


def create_continuous_enrichment_plot(
    k_values, enrichment_matrix, cat_names, profile_stats,
    prefix, analysis_type='BP_MF',
    highlight_categories=None, default_tiers=(50, 100, 250, 500),
    max_categories_highlight=5, figsize=(12, 7)):
    """Create continuous enrichment curve plot with null envelope.

    Args:
        k_values: 1D array of gene ranks
        enrichment_matrix: (K, C) array of fold changes
        cat_names: list of category names
        profile_stats: dict from permutation_global_test
        prefix: output file prefix
        analysis_type: 'BP_MF' or 'CC'
        highlight_categories: list of categories to highlight (auto-selected if None)
        default_tiers: vertical reference lines at default gene ranks
        max_categories_highlight: max categories to highlight
        figsize: figure size
    """
    n_cats = len(cat_names)

    # Auto-select highlight categories: highest supremum statistic
    if highlight_categories is None:
        sup_scores = [(cat, profile_stats[cat]['supremum_obs'])
                      for cat in cat_names if cat in profile_stats]
        sup_scores.sort(key=lambda x: x[1], reverse=True)
        highlight_categories = [c for c, _ in sup_scores[:max_categories_highlight]]

    fig, ax = plt.subplots(figsize=figsize)

    # Null envelope (use first category's envelope as representative,
    # or compute global envelope)
    # Actually, draw per-category envelopes only for highlighted categories
    # For the global reference, shade a light band from the average envelope
    all_lower = []
    all_upper = []
    for cat in cat_names:
        if cat in profile_stats:
            all_lower.append(profile_stats[cat]['null_envelope_lower'])
            all_upper.append(profile_stats[cat]['null_envelope_upper'])
    if all_lower:
        mean_lower = np.mean(all_lower, axis=0)
        mean_upper = np.mean(all_upper, axis=0)
        ax.fill_between(k_values, mean_lower, mean_upper,
                        color='grey', alpha=0.12, label='95% null envelope (mean)')

    # Plot non-highlighted categories in grey
    for i, cat in enumerate(cat_names):
        if cat not in highlight_categories:
            ax.plot(k_values, enrichment_matrix[:, i],
                    color='#CCCCCC', linewidth=0.5, alpha=0.6, zorder=1)

    # Plot highlighted categories with colour
    highlight_handles = []
    for j, cat in enumerate(highlight_categories):
        if cat not in cat_names:
            continue
        i = cat_names.index(cat)
        colour = CATEGORY_COLOURS[j % len(CATEGORY_COLOURS)]

        ax.plot(k_values, enrichment_matrix[:, i],
                color=colour, linewidth=2.0, alpha=0.9, zorder=3)

        # Mark regions outside null envelope
        if cat in profile_stats:
            env_lo = profile_stats[cat]['null_envelope_lower']
            env_hi = profile_stats[cat]['null_envelope_upper']
            curve = enrichment_matrix[:, i]
            outside = (curve > env_hi) | (curve < env_lo)
            if np.any(outside):
                ax.scatter(k_values[outside], curve[outside],
                           color=colour, s=3, alpha=0.4, zorder=4)

        # Build legend label with p-value
        p_label = ''
        if cat in profile_stats:
            p_sup = profile_stats[cat]['supremum_p']
            if p_sup < 0.001:
                p_label = ' (p < 0.001)'
            elif p_sup < 0.05:
                p_label = f' (p = {p_sup:.3f})'
        highlight_handles.append(
            Line2D([0], [0], color=colour, linewidth=2, label=f'{cat}{p_label}'))

    # Reference line at FC = 1.0
    ax.axhline(y=1.0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)

    # Default tier markers
    for tier in default_tiers:
        if tier <= k_values[-1]:
            ax.axvline(x=tier, color='grey', linestyle=':', linewidth=0.5, alpha=0.4)

    # Axes
    ax.set_xlabel('Gene rank k', fontsize=12)
    ax.set_ylabel('Fold enrichment E_C(t)', fontsize=12)
    ax.set_title('Continuous Enrichment Profiles', fontsize=14)

    # Legend
    grey_handle = Line2D([0], [0], color='#CCCCCC', linewidth=1,
                         label=f'Other categories ({n_cats - len(highlight_categories)})')
    env_handle = plt.Rectangle((0, 0), 1, 1, fc='grey', alpha=0.15,
                               label='95% null envelope')
    all_handles = highlight_handles + [grey_handle, env_handle]
    ax.legend(handles=all_handles, fontsize=8, loc='upper right',
              framealpha=0.9, ncol=1)

    ax.set_xlim(k_values[0], k_values[-1])
    ax.tick_params(labelsize=10)

    plt.tight_layout()

    for ext in ['png', 'svg']:
        out = f"{prefix}_{analysis_type}_continuous_enrichment.{ext}"
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Saved: {out}")

    plt.close(fig)


def create_continuous_dkl_plot(
    k_values, dkl_values, prefix, analysis_type='BP_MF',
    default_tiers=(50, 100, 250, 500), figsize=(10, 5)):
    """Plot continuous D_KL gradient across the expression hierarchy."""
    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(k_values, dkl_values, color='#3C5488', linewidth=2)

    # Log-linear fit
    mask = dkl_values > 0
    if np.sum(mask) > 10:
        from scipy.stats import linregress
        log_k = np.log(k_values[mask])
        log_dkl = np.log(dkl_values[mask])
        result = linregress(log_k, log_dkl)
        fit_dkl = np.exp(result.intercept + result.slope * log_k)
        ax.plot(k_values[mask], fit_dkl, '--', color='#E64B35',
                linewidth=1.5, alpha=0.7,
                label=f'Log-linear fit (slope={result.slope:.2f}, '
                      f'R²={result.rvalue**2:.3f})')

    # Default tier markers
    for tier in default_tiers:
        if tier <= k_values[-1]:
            ax.axvline(x=tier, color='grey', linestyle=':', linewidth=0.5, alpha=0.4)

    ax.set_xlabel('Gene rank k', fontsize=12)
    ax.set_ylabel('D_KL (bits)', fontsize=12)
    ax.set_title('Continuous Functional Specialisation Gradient', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=10)
    ax.tick_params(labelsize=10)

    plt.tight_layout()

    for ext in ['png', 'svg']:
        out = f"{prefix}_{analysis_type}_continuous_dkl.{ext}"
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Saved: {out}")

    plt.close(fig)


def create_significance_landscape(
    k_values, enrichment_matrix, cat_names, profile_stats,
    prefix, analysis_type='BP_MF', figsize=(14, 6)):
    """Create significance landscape: heatmap of enrichment across k.

    Each row is a category, each column is a k value. Colour encodes
    fold change, and significance (outside null envelope) is overlaid.
    """
    n_cats = len(cat_names)

    # Sort categories by maximum fold change (descending)
    max_fc = [np.max(enrichment_matrix[:, i]) for i in range(n_cats)]
    sort_idx = np.argsort(max_fc)[::-1]
    sorted_cats = [cat_names[i] for i in sort_idx]
    sorted_matrix = enrichment_matrix[:, sort_idx].T  # (C, K)

    fig, ax = plt.subplots(figsize=figsize)

    # Diverging colourmap centred at FC=1
    vmax = min(np.percentile(sorted_matrix[sorted_matrix > 0], 98), 10)
    im = ax.imshow(sorted_matrix, aspect='auto', cmap='RdBu_r',
                   vmin=0, vmax=vmax,
                   extent=[k_values[0], k_values[-1], n_cats - 0.5, -0.5])

    # Mark significant regions (outside null envelope)
    for row_j, cat in enumerate(sorted_cats):
        if cat not in profile_stats:
            continue
        orig_i = cat_names.index(cat)
        env_lo = profile_stats[cat]['null_envelope_lower']
        env_hi = profile_stats[cat]['null_envelope_upper']
        curve = enrichment_matrix[:, orig_i]
        outside = (curve > env_hi) | (curve < env_lo)
        if np.any(outside):
            sig_k = k_values[outside]
            ax.scatter(sig_k, np.full_like(sig_k, row_j, dtype=float),
                       color='black', s=0.3, alpha=0.5, zorder=2)

    ax.set_yticks(range(n_cats))
    ax.set_yticklabels(sorted_cats, fontsize=7)
    ax.set_xlabel('Gene rank k', fontsize=11)
    ax.set_title('Enrichment Significance Landscape', fontsize=13)

    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Fold enrichment', fontsize=10)

    plt.tight_layout()

    for ext in ['png', 'svg']:
        out = f"{prefix}_{analysis_type}_significance_landscape.{ext}"
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Saved: {out}")

    plt.close(fig)


def create_profile_shape_plot(
    shape_stats, prefix, analysis_type='BP_MF', figsize=(8, 6)):
    """Plot profile shape classification: slope vs max enrichment."""
    cats = list(shape_stats.keys())
    slopes = [shape_stats[c]['slope'] for c in cats]
    max_fcs = [shape_stats[c]['max_enrichment'] for c in cats]
    shapes = [shape_stats[c]['shape_class'] for c in cats]

    colour_map = {
        'apex-concentrated': '#E64B35',
        'distributed': '#4DBBD5',
        'flat': '#7F7F7F',
        'absent': '#CCCCCC',
    }
    colours = [colour_map.get(s, '#7F7F7F') for s in shapes]

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(slopes, max_fcs, c=colours, s=60, edgecolors='black',
               linewidth=0.5, zorder=3)

    for i, cat in enumerate(cats):
        if max_fcs[i] > 1.5 or abs(slopes[i]) > 0.3:
            ax.annotate(cat, (slopes[i], max_fcs[i]),
                        fontsize=6, ha='center', va='bottom',
                        xytext=(0, 4), textcoords='offset points')

    ax.axvline(x=0, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axhline(y=1.0, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

    ax.set_xlabel('Profile slope (normalised)', fontsize=11)
    ax.set_ylabel('Maximum fold enrichment', fontsize=11)
    ax.set_title('Profile Shape Classification', fontsize=13)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E64B35', label='Apex-concentrated'),
        Patch(facecolor='#4DBBD5', label='Distributed'),
        Patch(facecolor='#7F7F7F', label='Flat'),
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc='upper right')

    plt.tight_layout()

    for ext in ['png', 'svg']:
        out = f"{prefix}_{analysis_type}_profile_shapes.{ext}"
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Saved: {out}")

    plt.close(fig)
