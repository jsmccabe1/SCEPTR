#!/usr/bin/env python3
"""
Interactive HTML report for SCEPTR cross-sample comparison.

Generates a self-contained Plotly.js dashboard with continuous enrichment
overlays, differential heatmap, concordance metrics, and methods
documentation. Matches the navy/cyan design of single-sample reports.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import json
import logging
import os
from datetime import datetime
from typing import Dict, Any, Optional, List

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.reporting.comparison_report')

# ---------------------------------------------------------------------------
# Colourblind-safe palette (same as interactive_report.py)
# ---------------------------------------------------------------------------
PALETTE = [
    '#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F',
    '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85',
    '#E377C2', '#BCBD22', '#17BECF', '#AEC7E8', '#FFBB78',
]

_MIN_GENES_CONTINUOUS = 25


def _format_pval(p: float) -> str:
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def _np_to_list(arr):
    if hasattr(arr, 'tolist'):
        return arr.tolist()
    return list(arr)


def _cat_colour_map(cat_names):
    return {cat: PALETTE[i % len(PALETTE)] for i, cat in enumerate(cat_names)}


# ---------------------------------------------------------------------------
# CSS - navy/cyan theme matching interactive_report.py
# ---------------------------------------------------------------------------

_CSS = """
:root {
    --primary: #0f2b46;
    --primary-light: #1a3d5c;
    --primary-dark: #091e33;
    --accent: #22d3ee;
    --accent-dim: #0ea5c9;
    --bg: #f0f2f5;
    --card: #ffffff;
    --text: #1a2332;
    --text-light: #546478;
    --border: #dce1e8;
    --enriched: #E64B35;
    --depleted: #4DBBD5;
    --sig-bg: #f0fdfa;
    --cond-a: #3C5488;
    --cond-b: #E64B35;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: 'Segoe UI', -apple-system, BlinkMacSystemFont, sans-serif;
    line-height: 1.6; color: var(--text);
    max-width: 1300px; margin: 0 auto; padding: 24px;
    background: var(--bg);
}
.header {
    background: linear-gradient(135deg, var(--primary), var(--primary-dark));
    color: white; padding: 36px 40px; border-radius: 12px;
    margin-bottom: 28px; position: relative; overflow: hidden;
}
.header::before {
    content: ''; position: absolute; top: -60%; right: -8%;
    width: 350px; height: 350px; border-radius: 50%;
    background: rgba(34,211,238,0.08);
}
.header::after {
    content: ''; position: absolute; bottom: -50%; left: -5%;
    width: 250px; height: 250px; border-radius: 50%;
    background: rgba(34,211,238,0.05);
}
.header h1 {
    font-size: 1.9em; font-weight: 700; margin-bottom: 6px;
    position: relative; z-index: 1;
}
.header .subtitle {
    opacity: 0.9; font-size: 1.1em; position: relative; z-index: 1;
    color: var(--accent);
}
.header .meta {
    opacity: 0.65; font-size: 0.85em; margin-top: 10px;
    position: relative; z-index: 1;
}

.section {
    background: var(--card); padding: 28px 32px;
    margin-bottom: 24px; border-radius: 10px;
    box-shadow: 0 1px 6px rgba(0,0,0,0.06);
}
.section h2 {
    color: var(--primary); font-size: 1.35em; font-weight: 700;
    border-bottom: 2px solid var(--accent-dim);
    padding-bottom: 10px; margin-bottom: 20px;
}
.section h3 { color: var(--text); font-size: 1.1em; margin: 24px 0 12px 0; }
.section p { margin-bottom: 12px; }

.stats-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 16px; margin: 16px 0 8px 0;
}
.stat-card {
    background: var(--bg); padding: 18px 16px;
    border-radius: 8px; text-align: center;
    border-left: 4px solid var(--primary);
    transition: transform 0.15s;
}
.stat-card:hover { transform: translateY(-2px); }
.stat-card.accent { border-left-color: var(--accent-dim); }
.stat-card.cond-a { border-left-color: var(--cond-a); }
.stat-card.cond-b { border-left-color: var(--cond-b); }
.stat-card .label {
    font-size: 0.78em; text-transform: uppercase;
    letter-spacing: 0.5px; color: var(--text-light); font-weight: 600;
}
.stat-card .value {
    font-size: 1.6em; font-weight: 700; color: var(--text);
    margin: 4px 0 2px 0;
}
.stat-card .detail { font-size: 0.78em; color: var(--text-light); }

.plot-container { margin: 20px 0; min-height: 400px; }
.export-btn {
    display: inline-block; padding: 6px 14px; margin: 4px 2px;
    background: var(--primary); color: white; border: none;
    border-radius: 5px; cursor: pointer; font-size: 0.8em;
    transition: background 0.15s;
}
.export-btn:hover { background: var(--primary-dark); }

table {
    width: 100%; border-collapse: collapse;
    margin: 14px 0; font-size: 0.85em;
}
th, td { padding: 8px 12px; text-align: left; border-bottom: 1px solid var(--border); }
th {
    background: var(--primary); color: white;
    font-weight: 600; font-size: 0.82em;
    text-transform: uppercase; letter-spacing: 0.3px;
}
tr:hover { background: #f8fafb; }
tr.sig { background: var(--sig-bg); }
.enriched-a { color: var(--cond-a); font-weight: 700; }
.enriched-b { color: var(--cond-b); font-weight: 700; }
.table-note { font-size: 0.8em; color: var(--text-light); margin-top: 8px; }

.methods-box {
    background: #f8fafb; padding: 18px 22px; border-radius: 8px;
    border-left: 3px solid var(--accent-dim);
    font-size: 0.88em; line-height: 1.7; color: var(--text-light);
}
.methods-box p { margin-bottom: 8px; }
.methods-box strong { color: var(--text); }
.footer {
    text-align: center; padding: 16px; color: var(--text-light);
    font-size: 0.78em;
}
.legend-inline {
    display: inline-block; width: 12px; height: 3px;
    vertical-align: middle; margin: 0 4px;
}
.legend-inline.solid { border-top: 3px solid; }
.legend-inline.dashed { border-top: 3px dashed; }
@media print {
    .export-btn { display: none; }
    .plot-container { page-break-inside: avoid; }
}
@media (max-width: 850px) {
    .stats-grid { grid-template-columns: 1fr 1fr; }
}
"""


# ---------------------------------------------------------------------------
# Plotly data builders
# ---------------------------------------------------------------------------

def _build_continuous_overlay(cont_a, cont_b, cat_names, colour_map,
                               label_a, label_b, cat_gene_counts_a=None,
                               cat_gene_counts_b=None):
    """Build Plotly traces for continuous enrichment overlay.

    Condition A: solid lines. Condition B: dashed lines. Same colour per category.
    Top 5 categories visible by default, rest legendonly.
    """
    from scipy.ndimage import gaussian_filter1d

    traces = []

    if cat_gene_counts_a is None:
        cat_gene_counts_a = {}
    if cat_gene_counts_b is None:
        cat_gene_counts_b = {}

    # Determine which categories have enough genes in at least one condition
    sparse_cats = set()
    for cat in cat_names:
        count_a = cat_gene_counts_a.get(cat, 9999)
        count_b = cat_gene_counts_b.get(cat, 9999)
        if count_a < _MIN_GENES_CONTINUOUS and count_b < _MIN_GENES_CONTINUOUS:
            sparse_cats.add(cat)

    if sparse_cats:
        logger.info(f"Suppressing continuous curves for {len(sparse_cats)} sparse "
                     f"categories (<{_MIN_GENES_CONTINUOUS} genes in both conditions)")

    # Rank categories by max enrichment across both conditions
    rankings = []
    for cat in cat_names:
        if cat in sparse_cats:
            continue
        max_fc = 0
        if cont_a is not None and cat in cont_a.columns:
            max_fc = max(max_fc, cont_a[cat].max())
        if cont_b is not None and cat in cont_b.columns:
            max_fc = max(max_fc, cont_b[cat].max())
        rankings.append((cat, max_fc))
    rankings.sort(key=lambda x: x[1], reverse=True)

    sigma = 2.0

    # Baseline = 1.0 line
    k_min = 10
    k_max = 500
    if cont_a is not None and 'k' in cont_a.columns:
        k_max = max(k_max, cont_a['k'].max())
    if cont_b is not None and 'k' in cont_b.columns:
        k_max = max(k_max, cont_b['k'].max())

    for rank, (cat, max_fc) in enumerate(rankings):
        colour = colour_map.get(cat, '#CCCCCC')
        visible = True if rank < 5 else 'legendonly'

        # Condition A (solid)
        if cont_a is not None and cat in cont_a.columns:
            k_a = cont_a['k'].values
            y_a = gaussian_filter1d(cont_a[cat].values.astype(float), sigma=sigma)
            traces.append({
                'x': _np_to_list(k_a), 'y': _np_to_list(y_a),
                'type': 'scatter', 'mode': 'lines',
                'name': f'{cat}',
                'legendgroup': cat,
                'line': {'color': colour, 'width': 2.5 if rank < 5 else 1.5,
                         'dash': 'solid'},
                'visible': visible,
                'hovertemplate': (f'<b>{cat}</b> ({label_a})<br>'
                                  f'k = %{{x:,.0f}}<br>'
                                  f'E_C(k) = %{{y:.2f}}x<extra></extra>'),
            })

        # Condition B (dashed)
        if cont_b is not None and cat in cont_b.columns:
            k_b = cont_b['k'].values
            y_b = gaussian_filter1d(cont_b[cat].values.astype(float), sigma=sigma)
            traces.append({
                'x': _np_to_list(k_b), 'y': _np_to_list(y_b),
                'type': 'scatter', 'mode': 'lines',
                'name': f'{cat} ({label_b})',
                'legendgroup': cat,
                'showlegend': False,
                'line': {'color': colour, 'width': 2.5 if rank < 5 else 1.5,
                         'dash': 'dash'},
                'visible': visible,
                'hovertemplate': (f'<b>{cat}</b> ({label_b})<br>'
                                  f'k = %{{x:,.0f}}<br>'
                                  f'E_C(k) = %{{y:.2f}}x<extra></extra>'),
            })

    # Baseline
    traces.append({
        'x': [k_min, k_max], 'y': [1.0, 1.0],
        'type': 'scatter', 'mode': 'lines',
        'line': {'color': 'rgba(0,0,0,0.3)', 'dash': 'dash', 'width': 1},
        'showlegend': False, 'hoverinfo': 'skip',
    })

    layout = {
        'xaxis': {'title': 'Gene rank k', 'type': 'log', 'gridcolor': '#f0f0f0'},
        'yaxis': {'title': 'Fold enrichment E<sub>C</sub>(k)', 'gridcolor': '#f0f0f0'},
        'legend': {'orientation': 'v', 'x': 1.02, 'y': 1, 'font': {'size': 11}},
        'hovermode': 'closest',
        'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        'margin': {'l': 60, 'r': 200, 't': 30, 'b': 50},
        'shapes': [{'type': 'line', 'x0': t, 'x1': t,
                    'y0': 0, 'y1': 1, 'yref': 'paper',
                    'line': {'color': '#ddd', 'width': 1, 'dash': 'dot'}}
                   for t in [50, 100, 250, 500]],
    }
    return traces, layout


def _build_heatmap_data(diff_df, cat_names, tiers, label_a, label_b):
    """Build Plotly heatmap traces for differential enrichment."""
    # Matrix: rows=categories, cols=tiers
    tier_labels = [f'Top {t}' for t in tiers]
    cats_sorted = sorted(cat_names)

    z_vals = []
    hover_text = []
    for cat in cats_sorted:
        row_z = []
        row_hover = []
        for tier in tiers:
            match = diff_df[(diff_df['Category'] == cat) & (diff_df['Tier'] == tier)]
            if len(match) == 0:
                row_z.append(0)
                row_hover.append(f'{cat}<br>Tier {tier}<br>No data')
            else:
                r = match.iloc[0]
                fc_diff = r['FC_Diff']
                fdr = r['Perm_FDR']
                sig = '*' if fdr < 0.05 else ''
                row_z.append(fc_diff)
                row_hover.append(
                    f'<b>{cat}</b><br>'
                    f'Tier {tier}<br>'
                    f'{label_a}: {r["FC_A"]:.2f}x<br>'
                    f'{label_b}: {r["FC_B"]:.2f}x<br>'
                    f'Delta: {fc_diff:+.2f}{sig}<br>'
                    f'FDR: {_format_pval(fdr)}'
                )
        z_vals.append(row_z)
        hover_text.append(row_hover)

    # Determine symmetric colour scale
    max_abs = max(abs(v) for row in z_vals for v in row) if z_vals else 1
    max_abs = max(max_abs, 0.5)

    traces = [{
        'z': z_vals,
        'x': tier_labels,
        'y': cats_sorted,
        'type': 'heatmap',
        'colorscale': [
            [0, '#3C5488'], [0.5, '#FFFFFF'], [1, '#E64B35']
        ],
        'zmin': -max_abs, 'zmax': max_abs,
        'text': hover_text,
        'hovertemplate': '%{text}<extra></extra>',
        'colorbar': {
            'title': {'text': f'Delta FC<br>({label_b} - {label_a})', 'side': 'right'},
            'thickness': 15,
        },
    }]

    layout = {
        'yaxis': {'autorange': 'reversed', 'tickfont': {'size': 11}},
        'xaxis': {'side': 'top', 'tickfont': {'size': 12}},
        'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        'margin': {'l': 280, 'r': 80, 't': 50, 'b': 30},
    }
    return traces, layout


# ---------------------------------------------------------------------------
# HTML table builders
# ---------------------------------------------------------------------------

def _build_apex_displacement_data(diff_df, cat_names, label_a, label_b,
                                   apex_tier=50):
    """Build Plotly traces and layout for the apex displacement chart.

    Shows which categories gained or lost share of the expression apex
    between two conditions  -  the resource reallocation view.
    """
    # Get fold changes at the apex tier
    tier_data = diff_df[diff_df['Tier'] == apex_tier].copy()
    if len(tier_data) == 0:
        # Try smallest available tier
        apex_tier = diff_df['Tier'].min()
        tier_data = diff_df[diff_df['Tier'] == apex_tier].copy()

    if len(tier_data) == 0:
        return None, None

    # Compute budget shares from fold changes
    fc_a = {}
    fc_b = {}
    for _, row in tier_data.iterrows():
        cat = row['Category']
        fc_a[cat] = max(float(row.get('FC_A', 0)), 0)
        fc_b[cat] = max(float(row.get('FC_B', 0)), 0)

    total_a = sum(fc_a.values()) or 1.0
    total_b = sum(fc_b.values()) or 1.0

    displacements = []
    for cat in fc_a:
        share_a = fc_a[cat] / total_a * 100
        share_b = fc_b[cat] / total_b * 100
        delta = share_b - share_a
        if abs(delta) > 0.1:
            displacements.append({
                'category': cat, 'share_a': share_a,
                'share_b': share_b, 'delta': delta
            })

    if not displacements:
        return None, None

    # Sort by delta
    displacements.sort(key=lambda x: x['delta'])

    cats = [d['category'] for d in displacements]
    deltas = [d['delta'] for d in displacements]
    colors = ['#c0392b' if d < 0 else '#2980b9' for d in deltas]
    hover = [
        f"{d['category']}<br>"
        f"{label_a}: {d['share_a']:.1f}%<br>"
        f"{label_b}: {d['share_b']:.1f}%<br>"
        f"Change: {d['delta']:+.1f}%"
        for d in displacements
    ]

    traces = [{
        'type': 'bar',
        'y': cats,
        'x': deltas,
        'orientation': 'h',
        'marker': {'color': colors, 'opacity': 0.85},
        'text': [f"{d:+.1f}%" for d in deltas],
        'textposition': 'outside',
        'textfont': {'size': 10},
        'hovertext': hover,
        'hoverinfo': 'text',
    }]

    layout = {
        'title': {
            'text': f'Apex Budget Reallocation (top {apex_tier} genes)',
            'font': {'size': 14},
        },
        'xaxis': {
            'title': f'Change in apex budget share (%)',
            'zeroline': True, 'zerolinewidth': 2, 'zerolinecolor': '#333',
        },
        'yaxis': {'automargin': True},
        'height': max(350, len(cats) * 28 + 120),
        'margin': {'l': 200, 'r': 80, 't': 60, 'b': 60},
        'showlegend': False,
        'plot_bgcolor': 'white',
        'annotations': [{
            'text': f'\u2190 Lost apex share | Gained apex share \u2192',
            'xref': 'paper', 'yref': 'paper',
            'x': 0.5, 'y': -0.12, 'showarrow': False,
            'font': {'size': 10, 'color': '#888'},
        }],
    }

    return traces, layout


def _build_concordance_table(concordance_data):
    if not concordance_data:
        return "<p>No concordance data available.</p>"

    rows = []
    for entry in concordance_data:
        rho = entry.get('Spearman_Rho', 0)
        ci_lo = entry.get('Spearman_CI_Lower', 0)
        ci_hi = entry.get('Spearman_CI_Upper', 0)
        jaccard = entry.get('Jaccard_Similarity', 0)
        n_shared = entry.get('N_Shared_Genes', 0)

        if rho > 0.7:
            rho_style = ' style="color:#00A087;font-weight:700;"'
        elif rho < 0.3:
            rho_style = ' style="color:#E64B35;font-weight:700;"'
        else:
            rho_style = ''

        rows.append(
            f'<tr>'
            f'<td><strong>Top {entry["Tier"]}</strong></td>'
            f'<td{rho_style}>{rho:.3f}</td>'
            f'<td>{ci_lo:.3f} - {ci_hi:.3f}</td>'
            f'<td>{jaccard:.3f}</td>'
            f'<td>{n_shared:,}</td>'
            f'</tr>'
        )

    return f"""
    <table>
        <thead><tr>
            <th>Tier</th><th>Spearman \u03C1</th><th>95% CI</th>
            <th>Jaccard</th><th>Shared Genes</th>
        </tr></thead>
        <tbody>{''.join(rows)}</tbody>
    </table>
    <p class="table-note">
        Spearman \u03C1 measures rank correlation of fold-change profiles.
        Jaccard index measures overlap of top-k gene sets.
    </p>
    """


def _build_differential_table(diff_df, tiers, label_a, label_b):
    if diff_df is None or len(diff_df) == 0:
        return "<p>No differential enrichment data available.</p>"

    categories = sorted(diff_df['Category'].unique())

    header = '<th>Category</th>'
    for t in tiers:
        header += f'<th>Top {t}</th>'

    rows = []
    for cat in categories:
        cells = f'<td><strong>{cat}</strong></td>'
        any_sig = False

        for t in tiers:
            row = diff_df[(diff_df['Category'] == cat) & (diff_df['Tier'] == t)]
            if len(row) == 0:
                cells += '<td style="color:var(--text-light)">-</td>'
                continue

            r = row.iloc[0]
            fc_a = r['FC_A']
            fc_b = r['FC_B']
            fc_diff = r['FC_Diff']
            fdr = r['Perm_FDR']
            sig = fdr < 0.05

            if sig:
                any_sig = True

            if fc_diff > 0.1:
                arrow = '<span style="color:var(--cond-b);">\u25B2</span>'
            elif fc_diff < -0.1:
                arrow = '<span style="color:var(--cond-a);">\u25BC</span>'
            else:
                arrow = '\u2022'

            sig_marker = ' *' if sig else ''

            cells += (f'<td>'
                      f'<span class="enriched-a">{fc_a:.2f}x</span> | '
                      f'<span class="enriched-b">{fc_b:.2f}x</span><br>'
                      f'<small>{arrow} {fc_diff:+.2f} ({_format_pval(fdr)}){sig_marker}</small>'
                      f'</td>')

        row_cls = ' class="sig"' if any_sig else ''
        rows.append(f'<tr{row_cls}>{cells}</tr>')

    return f"""
    <table>
        <thead><tr>{header}</tr></thead>
        <tbody>{''.join(rows)}</tbody>
    </table>
    <p class="table-note">
        * Significant (BH-adjusted permutation p &lt; 0.05).
        \u0394 = fold-change difference ({label_b} - {label_a}).
        <span style="color:var(--cond-b);">\u25B2</span> enriched in {label_b},
        <span style="color:var(--cond-a);">\u25BC</span> enriched in {label_a}.
        Values show <span class="enriched-a">{label_a}</span> |
        <span class="enriched-b">{label_b}</span>.
    </p>
    """


# ---------------------------------------------------------------------------
# Hero summary
# ---------------------------------------------------------------------------

def _generate_hero_summary(diff_df, concordance_data, label_a, label_b,
                            cont_a, cont_b, cat_names):
    """Generate a short plain-English summary of the comparison."""
    lines = []

    # Concordance
    if concordance_data:
        rho_250 = None
        for entry in concordance_data:
            if entry['Tier'] == 250:
                rho_250 = entry['Spearman_Rho']
        if rho_250 is None:
            rho_250 = concordance_data[-1]['Spearman_Rho']

        if rho_250 > 0.7:
            lines.append(f"{label_a} and {label_b} share broadly similar "
                         f"functional profiles (Spearman \u03C1 = {rho_250:.2f}).")
        elif rho_250 > 0.3:
            lines.append(f"{label_a} and {label_b} show moderately different "
                         f"functional profiles (Spearman \u03C1 = {rho_250:.2f}).")
        else:
            lines.append(f"{label_a} and {label_b} have substantially divergent "
                         f"functional profiles (Spearman \u03C1 = {rho_250:.2f}).")

    # Significant categories
    if diff_df is not None and len(diff_df) > 0:
        sig = diff_df[diff_df['Perm_FDR'] < 0.05]
        if len(sig) > 0:
            sig_cats = sorted(sig['Category'].unique())
            if len(sig_cats) <= 3:
                cat_str = ', '.join(f'<strong>{c}</strong>' for c in sig_cats)
            else:
                cat_str = (', '.join(f'<strong>{c}</strong>' for c in sig_cats[:3])
                           + f' and {len(sig_cats) - 3} others')
            lines.append(f"Significant differences detected in {cat_str}.")

            # Strongest difference
            best = sig.loc[sig['FC_Diff'].abs().idxmax()]
            direction = label_b if best['FC_Diff'] > 0 else label_a
            lines.append(
                f"The largest shift is in <strong>{best['Category']}</strong> "
                f"(enriched in {direction}, \u0394 = {best['FC_Diff']:+.2f} "
                f"at top {int(best['Tier'])}).")
        else:
            lines.append("No categories show statistically significant differences "
                         "between conditions at FDR < 0.05.")

    return " ".join(lines)


# ---------------------------------------------------------------------------
# Main report generator
# ---------------------------------------------------------------------------

def generate_comparison_report(
    diff_df,
    concordance_data: List[Dict[str, Any]],
    label_a: str,
    label_b: str,
    n_genes_a: int,
    n_genes_b: int,
    n_shared: int,
    n_permutations: int,
    tiers: List[int],
    output_prefix: str,
    figures_dir: Optional[str] = None,
    continuous_a: Optional[pd.DataFrame] = None,
    continuous_b: Optional[pd.DataFrame] = None,
    cat_gene_counts_a: Optional[Dict[str, int]] = None,
    cat_gene_counts_b: Optional[Dict[str, int]] = None,
) -> Optional[str]:
    """
    Generate self-contained interactive HTML comparison report.

    Args:
        diff_df: DataFrame with differential enrichment results
        concordance_data: List of per-tier concordance metric dicts
        label_a, label_b: Condition labels
        n_genes_a, n_genes_b: Gene counts per condition
        n_shared: Shared genes between conditions
        n_permutations: Number of permutations used
        tiers: List of tier sizes
        output_prefix: Output path prefix
        figures_dir: Directory containing figure files (legacy, unused)
        continuous_a: DataFrame with continuous enrichment for condition A
        continuous_b: DataFrame with continuous enrichment for condition B
        cat_gene_counts_a: Dict of {category: gene_count} for condition A
        cat_gene_counts_b: Dict of {category: gene_count} for condition B

    Returns:
        Path to generated HTML file, or None on error.
    """
    html_path = f"{output_prefix}_comparison_report.html"

    # Count significant results
    n_sig = int((diff_df['Perm_FDR'] < 0.05).sum()) if len(diff_df) > 0 else 0
    n_categories = len(diff_df['Category'].unique()) if len(diff_df) > 0 else 0
    cat_names = sorted(diff_df['Category'].unique()) if len(diff_df) > 0 else []
    colour_map = _cat_colour_map(cat_names)

    # Hero summary
    hero_html = _generate_hero_summary(
        diff_df, concordance_data, label_a, label_b,
        continuous_a, continuous_b, cat_names)

    # Concordance table
    concordance_table = _build_concordance_table(concordance_data)

    # Differential table
    differential_table = _build_differential_table(diff_df, tiers, label_a, label_b)

    # Build Plotly chart data
    has_continuous = continuous_a is not None or continuous_b is not None

    chart_js_blocks = []
    continuous_section = ''
    if has_continuous:
        cont_traces, cont_layout = _build_continuous_overlay(
            continuous_a, continuous_b, cat_names, colour_map,
            label_a, label_b, cat_gene_counts_a, cat_gene_counts_b)
        chart_js_blocks.append(
            f"Plotly.newPlot('cont-overlay', "
            f"{json.dumps(cont_traces)}, "
            f"{json.dumps(cont_layout)}, plotlyConfig);")
        continuous_section = f"""
<div class="section">
    <h2>Continuous Enrichment Comparison</h2>
    <p>Enrichment fold-change E<sub>C</sub>(k) computed continuously across all
    gene ranks. <span class="legend-inline solid" style="border-color:var(--cond-a);"></span>
    Solid lines = <strong>{label_a}</strong>,
    <span class="legend-inline dashed" style="border-color:var(--cond-b);"></span>
    dashed lines = <strong>{label_b}</strong>. Same colour = same category.
    Click legend entries to toggle visibility.</p>
    <div class="plot-container" id="cont-overlay"></div>
    <div style="text-align:right;">
        <button class="export-btn" onclick="exportPlot('cont-overlay','sceptr_comparison_enrichment','svg')">Export SVG</button>
        <button class="export-btn" onclick="exportPlot('cont-overlay','sceptr_comparison_enrichment','png')">Export PNG</button>
    </div>
</div>
"""

    # Apex displacement chart
    apex_traces, apex_layout = _build_apex_displacement_data(
        diff_df, cat_names, label_a, label_b)
    apex_section = ''
    if apex_traces is not None:
        chart_js_blocks.append(
            f"Plotly.newPlot('apex-displacement', "
            f"{json.dumps(apex_traces)}, "
            f"{json.dumps(apex_layout)}, plotlyConfig);")
        apex_section = f"""
<div class="section">
    <h2>Apex Budget Reallocation</h2>
    <p>How the cell's transcriptional investment at the expression apex shifted
    between <strong>{label_a}</strong> and <strong>{label_b}</strong>. Categories
    on the right gained apex share; categories on the left lost it. This shows
    what the cell traded away to invest in the activated programmes &mdash; a
    resource competition view that scalar enrichment scores cannot provide.</p>
    <div class="plot-container" id="apex-displacement"></div>
    <div style="text-align:right;">
        <button class="export-btn" onclick="exportPlot('apex-displacement','sceptr_apex_displacement','svg')">Export SVG</button>
        <button class="export-btn" onclick="exportPlot('apex-displacement','sceptr_apex_displacement','png')">Export PNG</button>
    </div>
</div>
"""

    # Heatmap
    heatmap_traces, heatmap_layout = _build_heatmap_data(
        diff_df, cat_names, tiers, label_a, label_b)
    # Dynamic height based on number of categories
    heatmap_height = max(400, len(cat_names) * 35 + 100)
    heatmap_layout['height'] = heatmap_height
    chart_js_blocks.append(
        f"Plotly.newPlot('diff-heatmap', "
        f"{json.dumps(heatmap_traces)}, "
        f"{json.dumps(heatmap_layout)}, plotlyConfig);")

    # JavaScript
    chart_renders = '\n    '.join(chart_js_blocks)

    try:
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SCEPTR Comparison: {label_a} vs {label_b}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>{_CSS}</style>
</head>
<body>

<div class="header">
    <h1>SCEPTR Cross-Sample Comparison</h1>
    <div class="subtitle">{label_a} vs {label_b}</div>
    <div class="meta">Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}
    &middot; {n_permutations:,} permutations &middot; {n_categories} categories
    &middot; SCEPTR Compare</div>
</div>

<div class="section">
    <h2>Summary</h2>
    <p>{hero_html}</p>
    <div class="stats-grid">
        <div class="stat-card cond-a">
            <div class="label">{label_a}</div>
            <div class="value">{n_genes_a:,}</div>
            <div class="detail">annotated genes</div>
        </div>
        <div class="stat-card cond-b">
            <div class="label">{label_b}</div>
            <div class="value">{n_genes_b:,}</div>
            <div class="detail">annotated genes</div>
        </div>
        <div class="stat-card">
            <div class="label">Shared Genes</div>
            <div class="value">{n_shared:,}</div>
            <div class="detail">{n_shared / max(n_genes_a, 1) * 100:.0f}% of {label_a}</div>
        </div>
        <div class="stat-card accent">
            <div class="label">Significant Differences</div>
            <div class="value">{n_sig}</div>
            <div class="detail">of {len(diff_df)} tests (FDR &lt; 0.05)</div>
        </div>
    </div>
</div>

{continuous_section}

{apex_section}

<div class="section">
    <h2>Concordance Overview</h2>
    <p>How similar are the two conditions' functional enrichment profiles?
    High Spearman \u03C1 indicates conserved functional investment patterns;
    low \u03C1 indicates divergent functional programmes.</p>
    {concordance_table}
</div>

<div class="section">
    <h2>Differential Enrichment</h2>
    <p>Gene-label permutation test ({n_permutations:,} permutations) with per-tier
    Benjamini-Hochberg FDR correction. Tests whether fold-change differences exceed
    those expected by chance when gene labels are randomly reassigned.</p>
    {differential_table}
    <h3>Differential Heatmap</h3>
    <p>Blue = enriched in <strong>{label_a}</strong>, red = enriched in
    <strong>{label_b}</strong>. Colour intensity reflects the magnitude of the
    fold-change difference.</p>
    <div class="plot-container" id="diff-heatmap"></div>
    <div style="text-align:right;">
        <button class="export-btn" onclick="exportPlot('diff-heatmap','sceptr_comparison_heatmap','svg')">Export SVG</button>
        <button class="export-btn" onclick="exportPlot('diff-heatmap','sceptr_comparison_heatmap','png')">Export PNG</button>
    </div>
</div>

<div class="section">
    <h2>Methods</h2>
    <div class="methods-box">
        <p><strong>Permutation test:</strong> Genes from both conditions are pooled
        and randomly re-partitioned into two groups matching original sizes. For each
        permutation, genes are re-ranked by TPM within each group, tier memberships
        are recomputed, and category fold changes are recalculated. The two-sided
        p-value is the fraction of permutations where the absolute fold-change
        difference equals or exceeds the observed absolute difference, computed as
        (count + 1) / (n_permutations + 1).</p>

        <p><strong>Multiple testing correction:</strong> Per-tier Benjamini-Hochberg
        FDR correction across all categories within each expression tier.</p>

        <p><strong>Concordance metrics:</strong> Spearman rank correlation of
        category fold-change profiles between conditions, with Fisher z-transform
        95% confidence intervals. Jaccard similarity of top-k gene sets measures
        overlap of the most highly expressed genes between conditions.</p>

        <p><strong>Continuous enrichment:</strong> E<sub>C</sub>(k) is computed at
        every gene rank k by expanding the top-k set one gene at a time. Curves are
        smoothed with a Gaussian filter (sigma = 2.0) for visual clarity. Solid lines
        represent {label_a}; dashed lines represent {label_b}.</p>

        <p><strong>Categorisation:</strong> Genes are assigned to functional categories
        using SCEPTR's dual-method approach (keyword matching + GO ID hierarchy overlap).
        Category assignments are computed independently per condition.</p>
    </div>
</div>

<div class="section">
    <h2>Interpretation Guide</h2>
    <div class="methods-box">
        <p><strong>Significant positive delta:</strong> Category is more enriched in
        {label_b} than {label_a} at that expression tier - {label_b} invests
        more transcriptional resources in this function among its most highly
        expressed genes.</p>

        <p><strong>Significant negative delta:</strong> Category is more enriched in
        {label_a} than {label_b}.</p>

        <p><strong>High concordance (Spearman \u03C1 > 0.7):</strong> The two
        conditions share a broadly similar functional investment pattern, differing
        in specific categories rather than globally.</p>

        <p><strong>Low concordance (Spearman \u03C1 < 0.3):</strong> The two
        conditions have fundamentally different functional investment profiles,
        suggesting a major biological shift.</p>
    </div>
</div>

<div class="footer">
    SCEPTR v1.0.0 &middot; Cross-Sample Comparison &middot;
    {datetime.now().strftime('%Y-%m-%d')}
</div>

<script>
document.addEventListener('DOMContentLoaded', function() {{
    var plotlyConfig = {{
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['lasso2d', 'select2d'],
        toImageButtonOptions: {{ format: 'svg', filename: 'sceptr_comparison', scale: 2 }}
    }};
    {chart_renders}
}});
function exportPlot(divId, filename, format) {{
    var opts = {{format: format, width: 1200, height: 700, scale: 3}};
    if (format === 'svg') {{ opts.scale = 1; }}
    Plotly.downloadImage(divId, Object.assign(opts, {{filename: filename}}));
}}
</script>

</body>
</html>"""

        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html)

        logger.info(f"Generated comparison report: {html_path}")
        return html_path

    except Exception as e:
        logger.error(f"Error generating comparison report: {e}")
        import traceback
        traceback.print_exc()
        return None
