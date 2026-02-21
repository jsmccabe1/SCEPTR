#!/usr/bin/env python3
"""
Comparison visualisation for SCEPTR cross-sample analysis.

Generates publication-quality figures for comparing functional enrichment
profiles between two conditions, including differential heatmaps, grouped
bar plots, volcano plots, and gradient overlay figures.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import os
import logging
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from typing import Dict, List, Tuple, Any, Optional

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.visualisation.comparison_charts')

# Consistent tier colour scheme
TIER_COLORS = {
    'top_50': '#1f77b4',
    'top_100': '#ff7f0e',
    'top_250': '#2ca02c',
    'top_500': '#d62728',
}

# Condition colours
COND_COLORS = ['#2980b9', '#c0392b']


def _save_figure(fig, path_prefix: str, dpi: int = 300) -> Tuple[str, str]:
    """Save figure as both PNG and SVG."""
    png = f"{path_prefix}.png"
    svg = f"{path_prefix}.svg"
    fig.savefig(png, dpi=dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(svg, bbox_inches='tight', facecolor='white', format='svg')
    plt.close(fig)
    logger.info(f"Saved figure: {png}")
    return png, svg


def create_radar_overlay(
    results_a: Dict[str, Dict[str, Any]],
    results_b: Dict[str, Dict[str, Any]],
    categories: List[str],
    label_a: str,
    label_b: str,
    output_prefix: str,
    tiers: Optional[List[str]] = None,
) -> Tuple[str, str]:
    """
    Create side-by-side radar plots for two conditions.

    Each panel shows all tiers overlaid for one condition, making it easy
    to compare the overall functional investment profiles.
    """
    logger.info("Creating radar overlay comparison...")

    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 10,
        'figure.facecolor': 'white',
    })

    if tiers is None:
        tiers = sorted(results_a.keys(), key=lambda x: int(x.split('_')[-1]))

    num_cats = len(categories)
    angles = np.linspace(0, 2 * np.pi, num_cats, endpoint=False).tolist()
    angles += angles[:1]

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(16, 7), subplot_kw=dict(projection='polar'),
        facecolor='white')

    tier_cols = list(TIER_COLORS.values())

    def _plot_panel(ax, results, title):
        max_pct = 0
        for i, tn in enumerate(tiers):
            tier_data = results.get(tn, {})
            counts = tier_data.get('counts', {})
            total = sum(counts.values()) or 1
            pcts = [counts.get(c, 0) / total * 100 for c in categories]
            max_pct = max(max_pct, max(pcts) if pcts else 0)
            pcts_closed = pcts + [pcts[0]]
            color = tier_cols[i % len(tier_cols)]
            ax.plot(angles, pcts_closed, linewidth=2.2, color=color,
                    marker='o', markersize=4, markeredgecolor='white',
                    markeredgewidth=0.8,
                    label=tn.replace('_', ' ').title(), zorder=10 - i)
            ax.fill(angles, pcts_closed, color=color, alpha=0.10)

        max_pct = max(20, ((int(max_pct) // 20) + 1) * 20)
        max_pct = min(100, max_pct)

        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_xticks(angles[:-1])
        fmt = []
        for c in categories:
            if len(c) > 20:
                w = c.split()
                c = ' '.join(w[:len(w) // 2]) + '\n' + ' '.join(w[len(w) // 2:])
            fmt.append(c)
        ax.set_xticklabels(fmt, fontsize=8, color='#333333')
        yticks = list(range(0, max_pct + 20, 20))
        ax.set_yticks(yticks)
        ax.set_yticklabels([f'{y}%' for y in yticks], fontsize=8,
                           color='#666666', alpha=0.7)
        ax.set_ylim(0, max_pct)
        ax.grid(True, linestyle='-', alpha=0.3, linewidth=0.5)
        ax.set_facecolor('white')
        ax.set_title(title, fontsize=12, fontweight='bold', pad=15)

    _plot_panel(ax_a, results_a, label_a)
    _plot_panel(ax_b, results_b, label_b)

    # Legend from first panel
    handles, labels = ax_a.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=len(tiers),
               frameon=True, fontsize=9, facecolor='white', edgecolor='#cccccc',
               bbox_to_anchor=(0.5, -0.02))

    fig.suptitle('Cross-Sample Functional Profile Comparison',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    return _save_figure(fig, f"{output_prefix}_radar_overlay")


def create_differential_heatmap(
    diff_df,
    categories: List[str],
    tiers: List[int],
    label_a: str,
    label_b: str,
    output_prefix: str,
) -> Tuple[str, str]:
    """
    Create a heatmap of fold-change differences (Category x Tier).

    Blue = enriched in condition A, Red = enriched in condition B.
    Significance marked with asterisks.
    """
    logger.info("Creating differential enrichment heatmap...")

    plt.style.use('default')

    # Build matrix
    cat_order = sorted(categories)
    tier_labels = [f'Tier {t}' for t in tiers]
    matrix = np.zeros((len(cat_order), len(tiers)))
    sig_mask = np.zeros_like(matrix, dtype=bool)

    for ci, cat in enumerate(cat_order):
        for ti, t in enumerate(tiers):
            row = diff_df[(diff_df['Category'] == cat) & (diff_df['Tier'] == t)]
            if len(row) > 0:
                matrix[ci, ti] = row['FC_Diff'].values[0]
                sig_mask[ci, ti] = row['Perm_FDR'].values[0] < 0.05

    vmax = max(abs(matrix.min()), abs(matrix.max()), 0.5)

    fig, ax = plt.subplots(figsize=(8, max(6, len(cat_order) * 0.45)),
                           facecolor='white')

    cmap = plt.cm.RdBu_r
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.imshow(matrix, aspect='auto', cmap=cmap, norm=norm)

    ax.set_xticks(range(len(tiers)))
    ax.set_xticklabels(tier_labels, fontsize=10)
    ax.set_yticks(range(len(cat_order)))
    ax.set_yticklabels(cat_order, fontsize=9)

    # Annotate cells
    for ci in range(len(cat_order)):
        for ti in range(len(tiers)):
            val = matrix[ci, ti]
            sig = sig_mask[ci, ti]
            text_color = 'white' if abs(val) > vmax * 0.6 else 'black'
            marker = f'{val:+.2f}' + (' *' if sig else '')
            ax.text(ti, ci, marker, ha='center', va='center',
                    fontsize=7.5, color=text_color, fontweight='bold' if sig else 'normal')

    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label(f'FC Difference ({label_b} \u2212 {label_a})', fontsize=10)

    ax.set_title('Differential Enrichment Heatmap', fontsize=13, fontweight='bold',
                 pad=12)
    ax.set_xlabel('Expression Tier', fontsize=11)

    plt.figtext(0.5, -0.01,
                f'Blue = enriched in {label_a}  |  Red = enriched in {label_b}  |  * FDR < 0.05',
                ha='center', fontsize=8, color='#666666', style='italic')

    plt.tight_layout()
    return _save_figure(fig, f"{output_prefix}_differential_heatmap")


def create_fc_barplot(
    diff_df,
    categories: List[str],
    tiers: List[int],
    label_a: str,
    label_b: str,
    output_prefix: str,
) -> Tuple[str, str]:
    """
    Create grouped bar plot showing FC for both conditions side by side.
    """
    logger.info("Creating grouped FC bar plot...")

    plt.style.use('default')

    # Use top categories by max absolute difference
    cat_max = {}
    for cat in categories:
        rows = diff_df[diff_df['Category'] == cat]
        cat_max[cat] = rows['FC_Diff'].abs().max() if len(rows) > 0 else 0
    top_cats = sorted(cat_max.keys(), key=lambda c: cat_max[c], reverse=True)[:10]

    n_cats = len(top_cats)
    n_tiers = len(tiers)

    fig, axes = plt.subplots(1, n_tiers, figsize=(5 * n_tiers, max(5, n_cats * 0.4)),
                             sharey=True, facecolor='white')
    if n_tiers == 1:
        axes = [axes]

    for ti, (ax, t) in enumerate(zip(axes, tiers)):
        y = np.arange(n_cats)
        bar_h = 0.35

        fc_a_vals = []
        fc_b_vals = []
        sig_vals = []
        for cat in top_cats:
            row = diff_df[(diff_df['Category'] == cat) & (diff_df['Tier'] == t)]
            if len(row) > 0:
                fc_a_vals.append(row['FC_A'].values[0])
                fc_b_vals.append(row['FC_B'].values[0])
                sig_vals.append(row['Perm_FDR'].values[0] < 0.05)
            else:
                fc_a_vals.append(0)
                fc_b_vals.append(0)
                sig_vals.append(False)

        ax.barh(y + bar_h / 2, fc_a_vals, bar_h, color=COND_COLORS[0],
                alpha=0.85, label=label_a, edgecolor='white', linewidth=0.5)
        ax.barh(y - bar_h / 2, fc_b_vals, bar_h, color=COND_COLORS[1],
                alpha=0.85, label=label_b, edgecolor='white', linewidth=0.5)

        # Mark significant
        for j, sig in enumerate(sig_vals):
            if sig:
                ax.text(max(fc_a_vals[j], fc_b_vals[j]) + 0.05, j, '*',
                        fontsize=14, fontweight='bold', color='#e74c3c',
                        va='center')

        ax.axvline(x=1.0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
        ax.set_yticks(y)
        if ti == 0:
            ax.set_yticklabels([c[:25] for c in top_cats], fontsize=9)
        ax.set_xlabel('Fold Change', fontsize=10)
        ax.set_title(f'Tier {t}', fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='x', linestyle='--', alpha=0.3)
        if ti == 0:
            ax.legend(fontsize=9, loc='lower right')

    fig.suptitle(f'Fold Enrichment: {label_a} vs {label_b}',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    return _save_figure(fig, f"{output_prefix}_fc_barplot")


def create_volcano_plot(
    diff_df,
    label_a: str,
    label_b: str,
    output_prefix: str,
    fdr_threshold: float = 0.05,
) -> Tuple[str, str]:
    """
    Volcano-style plot: FC_Diff vs -log10(perm_p_value).

    One point per category x tier, significant points labelled.
    """
    logger.info("Creating volcano plot...")

    plt.style.use('default')

    fig, ax = plt.subplots(figsize=(10, 7), facecolor='white')

    x = diff_df['FC_Diff'].values
    y = -np.log10(diff_df['Perm_P_Value'].clip(lower=1e-10).values)
    sig = diff_df['Perm_FDR'].values < fdr_threshold
    tiers_col = diff_df['Tier'].values

    tier_vals = sorted(diff_df['Tier'].unique())
    tier_color_map = {}
    tier_cols_list = list(TIER_COLORS.values())
    for i, t in enumerate(tier_vals):
        tier_color_map[t] = tier_cols_list[i % len(tier_cols_list)]

    # Plot non-significant
    for t in tier_vals:
        mask = (tiers_col == t) & (~sig)
        ax.scatter(x[mask], y[mask], c=tier_color_map[t], alpha=0.35, s=40,
                   edgecolors='none', label=f'Tier {t}')

    # Plot significant
    for t in tier_vals:
        mask = (tiers_col == t) & sig
        ax.scatter(x[mask], y[mask], c=tier_color_map[t], alpha=1.0, s=80,
                   edgecolors='black', linewidths=0.8, zorder=5)

    # Label significant points
    sig_df = diff_df[diff_df['Perm_FDR'] < fdr_threshold]
    for _, row in sig_df.iterrows():
        ax.annotate(f"{row['Category'][:18]}\n(T{row['Tier']})",
                    (row['FC_Diff'], -np.log10(max(row['Perm_P_Value'], 1e-10))),
                    fontsize=7, ha='center', va='bottom',
                    textcoords='offset points', xytext=(0, 6),
                    arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

    # Threshold line
    ax.axhline(y=-np.log10(fdr_threshold), color='gray', linestyle='--',
               alpha=0.5, linewidth=0.8)
    ax.axvline(x=0, color='gray', linestyle='-', alpha=0.3, linewidth=0.8)

    ax.set_xlabel(f'FC Difference ({label_b} \u2212 {label_a})', fontsize=11)
    ax.set_ylabel('-log10(Permutation p-value)', fontsize=11)
    ax.set_title(f'Volcano Plot: {label_a} vs {label_b}',
                 fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.2)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=9, loc='upper right',
              frameon=True, facecolor='white', edgecolor='#cccccc')

    plt.tight_layout()
    return _save_figure(fig, f"{output_prefix}_volcano")


def create_gradient_overlay(
    diff_df,
    categories: List[str],
    tiers: List[int],
    label_a: str,
    label_b: str,
    output_prefix: str,
    fdr_threshold: float = 0.05,
) -> Tuple[str, str]:
    """
    Gradient overlay: line plot showing FC across tiers for both conditions.

    Each significant category gets a subplot showing how enrichment changes
    across the expression gradient for both conditions.
    """
    logger.info("Creating gradient overlay plot...")

    plt.style.use('default')

    # Select categories with at least one significant tier
    sig_cats = diff_df[diff_df['Perm_FDR'] < fdr_threshold]['Category'].unique()
    if len(sig_cats) == 0:
        # Fall back to top categories by max FC difference
        cat_max = {}
        for cat in categories:
            rows = diff_df[diff_df['Category'] == cat]
            cat_max[cat] = rows['FC_Diff'].abs().max() if len(rows) > 0 else 0
        sig_cats = sorted(cat_max.keys(), key=lambda c: cat_max[c], reverse=True)[:6]

    n_cats = len(sig_cats)
    if n_cats == 0:
        logger.warning("No categories to plot for gradient overlay")
        # Create minimal placeholder
        fig, ax = plt.subplots(figsize=(8, 4), facecolor='white')
        ax.text(0.5, 0.5, 'No significant differential enrichment detected',
                ha='center', va='center', fontsize=12, color='#666666')
        ax.axis('off')
        return _save_figure(fig, f"{output_prefix}_gradient_overlay")

    ncols = min(3, n_cats)
    nrows = (n_cats + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4 * nrows),
                             facecolor='white', squeeze=False)

    for idx, cat in enumerate(sig_cats):
        row_i = idx // ncols
        col_i = idx % ncols
        ax = axes[row_i][col_i]

        cat_data = diff_df[diff_df['Category'] == cat].sort_values('Tier')
        tier_x = cat_data['Tier'].values

        fc_a = cat_data['FC_A'].values
        fc_b = cat_data['FC_B'].values
        sigs = cat_data['Perm_FDR'].values < fdr_threshold

        ax.plot(tier_x, fc_a, '-o', color=COND_COLORS[0], linewidth=2,
                markersize=7, label=label_a, markeredgecolor='white',
                markeredgewidth=1)
        ax.plot(tier_x, fc_b, '-s', color=COND_COLORS[1], linewidth=2,
                markersize=7, label=label_b, markeredgecolor='white',
                markeredgewidth=1)

        # Highlight significant tiers
        for i, (t, s) in enumerate(zip(tier_x, sigs)):
            if s:
                ax.axvspan(t - 15, t + 15, alpha=0.1, color='#e74c3c')
                ax.text(t, max(fc_a[i], fc_b[i]) * 1.05, '*',
                        fontsize=14, ha='center', color='#e74c3c',
                        fontweight='bold')

        ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
        ax.set_xlabel('Expression Tier (top N genes)', fontsize=9)
        ax.set_ylabel('Fold Change', fontsize=9)
        ax.set_title(cat[:30], fontsize=10, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(alpha=0.2)
        if idx == 0:
            ax.legend(fontsize=8, loc='best')

    # Hide unused axes
    for idx in range(n_cats, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis('off')

    fig.suptitle(f'Expression Gradient Comparison: {label_a} vs {label_b}',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    return _save_figure(fig, f"{output_prefix}_gradient_overlay")
