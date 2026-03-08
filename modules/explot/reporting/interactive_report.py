#!/usr/bin/env python3
"""
Interactive HTML report generation for SCEPTR ExPlot.

Generates self-contained interactive dashboard with Plotly.js charts,
continuous enrichment curves, D_KL specialisation, category report cards,
and publication-ready figure export.

Supports combined reports: functional + cellular in one HTML file.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import json
import logging
import os
from datetime import datetime
from typing import Dict, Any, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger('sceptr.explot.reporting.interactive')

# ---------------------------------------------------------------------------
# Colourblind-safe palette (15 colours, consistent across all panels)
# ---------------------------------------------------------------------------
PALETTE = [
    '#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F',
    '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85',
    '#E377C2', '#BCBD22', '#17BECF', '#AEC7E8', '#FFBB78',
]


def _format_pval(p: float) -> str:
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def _np_to_list(arr):
    """Safely convert numpy array to Python list for JSON serialisation."""
    if hasattr(arr, 'tolist'):
        return arr.tolist()
    return list(arr)


def _categorisation_context(pct_cat, n_categories):
    """Return a short contextual note for the categorisation percentage."""
    if n_categories <= 20:
        if pct_cat < 65:
            return (f" -expected for a targeted {n_categories}-category set; "
                    "enrichment analysis uses tier-level coverage, not genome-wide")
        return ""
    else:
        if pct_cat < 70:
            return (" -lower than typical for a broad category set; "
                    "consider whether annotations are sparse for this organism")
        return ""


def _cat_colour(cat_names, idx):
    return PALETTE[idx % len(PALETTE)]


def _cat_colour_map(cat_names):
    return {cat: PALETTE[i % len(PALETTE)] for i, cat in enumerate(cat_names)}


# ---------------------------------------------------------------------------
# Report data serialisation (for combined reports)
# ---------------------------------------------------------------------------
class _NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)


def save_report_data(results, all_results, cont_results, output_path):
    """Save analysis results to JSON for later combined report generation."""
    data = {
        'results': results,
        'all_results': {
            'total_categorised': all_results.get('total_categorised', 0),
            'uncategorised': all_results.get('uncategorised', 0),
            'multi_category_genes': all_results.get('multi_category_genes', 0),
            'method_summary': all_results.get('method_summary', {}),
            'gene_categories': all_results.get('gene_categories', {}),
        },
    }
    if cont_results:
        data['cont_results'] = {
            'k_values': cont_results['k_values'],
            'enrichment_matrix': cont_results['enrichment_matrix'],
            'dkl_values': cont_results['dkl_values'],
            'cat_names': cont_results['cat_names'],
            'profile_stats': cont_results.get('profile_stats', {}),
            'shape_stats': cont_results.get('shape_stats', {}),
        }

    with open(output_path, 'w') as f:
        json.dump(data, f, cls=_NumpyEncoder)
    logger.info(f"Saved report data: {output_path}")


def load_report_data(path):
    """Load analysis results from JSON, converting arrays back to numpy."""
    with open(path) as f:
        data = json.load(f)

    if 'cont_results' in data and data['cont_results']:
        cr = data['cont_results']
        cr['k_values'] = np.array(cr['k_values'])
        cr['enrichment_matrix'] = np.array(cr['enrichment_matrix'])
        cr['dkl_values'] = np.array(cr['dkl_values'])
        # Convert null envelopes back to numpy
        for cat, ps in cr.get('profile_stats', {}).items():
            for key in ['null_envelope_lower', 'null_envelope_upper']:
                if key in ps:
                    ps[key] = np.array(ps[key])

    return data


# ---------------------------------------------------------------------------
# Hero summary generator
# ---------------------------------------------------------------------------
def _generate_hero_summary(analysis_blocks):
    """Generate a plain-English summary across all analysis blocks."""
    all_top_cats = []

    for block in analysis_blocks:
        cont_results = block.get('cont_results')
        results = block['results']
        tier_names = block['tier_names']
        cat_names = block['cat_names']
        label = block['section_title']

        if cont_results and 'shape_stats' in cont_results:
            for cat in cat_names:
                ss = cont_results['shape_stats'].get(cat, {})
                max_fc = ss.get('max_enrichment', 0)
                ps = cont_results.get('profile_stats', {}).get(cat, {})
                sup_p = ps.get('supremum_p', 1.0)
                if max_fc > 1.3 and sup_p < 0.05:
                    all_top_cats.append((cat, max_fc, label,
                                        ss.get('shape_class', 'unknown')))
        else:
            default_tier = 'top_250' if 'top_250' in tier_names else tier_names[-1]
            enrichment = results.get(default_tier, {}).get('enrichment', {})
            for c, s in enrichment.items():
                if s.get('significant', False) and 'Uncharacterised' not in c:
                    all_top_cats.append((c, s.get('fold_change', 0), label, 'unknown'))

    all_top_cats.sort(key=lambda x: x[1], reverse=True)

    if not all_top_cats:
        return ("No categories show strong significant enrichment -this transcriptome "
                "has a relatively uniform functional distribution.")

    lines = []
    top3 = all_top_cats[:3]
    names = [f"<strong>{c[0]}</strong> ({c[1]:.1f}×)" for c in top3]
    if len(names) == 3:
        cat_str = f"{names[0]}, {names[1]}, and {names[2]}"
    elif len(names) == 2:
        cat_str = f"{names[0]} and {names[1]}"
    else:
        cat_str = names[0]
    lines.append(f"The most prominent programmes in this transcriptome are {cat_str}.")

    apex_count = sum(1 for c in all_top_cats[:5] if c[3] == 'apex-concentrated')
    dist_count = sum(1 for c in all_top_cats[:5] if c[3] == 'distributed')
    if apex_count > dist_count:
        lines.append(
            "Enrichment is concentrated among the most highly expressed genes.")
    elif dist_count > apex_count:
        lines.append(
            "Enrichment is broadly distributed across expression levels.")

    return " ".join(lines)


# ---------------------------------------------------------------------------
# Transcriptome Overview (landscape data)
# ---------------------------------------------------------------------------
def _compute_landscape_data(df_sorted):
    """Compute landscape statistics from the sorted expression dataframe."""
    tpm = df_sorted['TPM'].values.copy()
    n = len(tpm)
    total = tpm.sum()
    if total == 0 or n == 0:
        return None

    # Gini coefficient
    sorted_tpm = np.sort(tpm)
    index = np.arange(1, n + 1)
    gini = (2 * np.sum(index * sorted_tpm) - (n + 1) * np.sum(sorted_tpm)) / (n * np.sum(sorted_tpm))

    # Lorenz curve (ascending sort for cumulative)
    cumulative = np.cumsum(sorted_tpm) / total
    x_lorenz = np.arange(1, n + 1) / n

    # Top-gene concentration (descending)
    tpm_desc = np.sort(tpm)[::-1]
    cum_desc = np.cumsum(tpm_desc) / total
    top10_pct = cum_desc[min(9, n - 1)] * 100
    top50_pct = cum_desc[min(49, n - 1)] * 100
    top100_pct = cum_desc[min(99, n - 1)] * 100

    # Expression histogram (log10 of nonzero)
    tpm_nonzero = tpm[tpm > 0]
    log_tpm = np.log10(tpm_nonzero) if len(tpm_nonzero) > 0 else np.array([])
    hist_counts, hist_edges = np.histogram(log_tpm, bins=50) if len(log_tpm) > 0 else ([], [])

    # Annotation quality by decile
    n_per_decile = max(1, n // 10)
    decile_stats = []
    has_uniprot = 'uniprot_id' in df_sorted.columns
    go_cols = [c for c in df_sorted.columns if c.startswith('GO_')]

    for i in range(10):
        start = i * n_per_decile
        end = min((i + 1) * n_per_decile, n)
        subset = df_sorted.iloc[start:end]
        n_sub = len(subset)
        if n_sub == 0:
            continue

        diamond_pct = 0.0
        if has_uniprot:
            has_hit = subset['uniprot_id'].notna() & (subset['uniprot_id'].astype(str).str.strip() != '')
            diamond_pct = round(has_hit.sum() / n_sub * 100, 1)

        go_pct = 0.0
        if go_cols:
            has_go = pd.Series([False] * n_sub, index=subset.index)
            for col in go_cols:
                has_go = has_go | (subset[col].notna() & (subset[col].astype(str).str.strip() != '') &
                                   (subset[col].astype(str).str.strip() != 'nan'))
            go_pct = round(has_go.sum() / n_sub * 100, 1)

        decile_stats.append({
            'decile': i + 1, 'diamond_pct': diamond_pct, 'go_pct': go_pct,
        })

    # Overall annotation stats
    overall_diamond_pct = 0.0
    if has_uniprot:
        has_hit = df_sorted['uniprot_id'].notna() & (df_sorted['uniprot_id'].astype(str).str.strip() != '')
        overall_diamond_pct = round(has_hit.sum() / n * 100, 1)
    overall_go_pct = 0.0
    if go_cols:
        has_go = pd.Series([False] * n, index=df_sorted.index)
        for col in go_cols:
            has_go = has_go | (df_sorted[col].notna() & (df_sorted[col].astype(str).str.strip() != '') &
                               (df_sorted[col].astype(str).str.strip() != 'nan'))
        overall_go_pct = round(has_go.sum() / n * 100, 1)

    # Taxonomic distribution
    tax_data = None
    if 'organism' in df_sorted.columns:
        org = df_sorted[df_sorted['organism'].notna() &
                        (df_sorted['organism'].astype(str).str.strip() != '')]
        if len(org) > 0:
            org_counts = org['organism'].value_counts()
            top_n = 15
            top_orgs = org_counts.head(top_n)
            other_count = org_counts.iloc[top_n:].sum() if len(org_counts) > top_n else 0
            other_n_species = len(org_counts) - top_n if len(org_counts) > top_n else 0
            tax_data = {
                'names': list(top_orgs.index),
                'counts': [int(c) for c in top_orgs.values],
                'total_hits': len(org),
                'unique_organisms': len(org_counts),
                'other_count': int(other_count),
                'other_n_species': other_n_species,
            }

    # Subsample Lorenz curve for Plotly (every 20th point to keep it responsive)
    step = max(1, n // 500)
    lorenz_x = x_lorenz[::step].tolist()
    lorenz_y = cumulative[::step].tolist()
    # Ensure endpoints
    if lorenz_x[0] != 0:
        lorenz_x.insert(0, 0.0)
        lorenz_y.insert(0, 0.0)
    if lorenz_x[-1] != 1.0:
        lorenz_x.append(1.0)
        lorenz_y.append(1.0)

    return {
        'gini': round(gini, 3),
        'top10_pct': round(top10_pct, 1),
        'top50_pct': round(top50_pct, 1),
        'top100_pct': round(top100_pct, 1),
        'median_tpm': round(float(np.median(tpm)), 2),
        'mean_tpm': round(float(np.mean(tpm)), 2),
        'total_genes': n,
        'expressed_genes': int((tpm > 0).sum()),
        'lorenz_x': lorenz_x,
        'lorenz_y': lorenz_y,
        'hist_counts': [int(c) for c in hist_counts],
        'hist_edges': [round(float(e), 3) for e in hist_edges],
        'decile_stats': decile_stats,
        'overall_diamond_pct': overall_diamond_pct,
        'overall_go_pct': overall_go_pct,
        'tax_data': tax_data,
    }


def _gini_interpretation(gini):
    if gini > 0.85:
        return ("Highly concentrated -a small number of genes dominate expression. "
                "Common in specialised cells or organisms under stress.")
    elif gini > 0.7:
        return ("Moderately concentrated -typical for differentiated tissues "
                "with active translation and metabolism.")
    elif gini > 0.5:
        return "Moderate distribution -relatively diverse transcriptional activity."
    else:
        return "Broadly distributed -many genes contribute to total expression."


def _build_landscape_section(landscape_data, chart_render_calls,
                              chart_tab_map=None, tab_id=None):
    """Build the Transcriptome Overview section HTML with Plotly charts."""
    if landscape_data is None:
        return ""

    ld = landscape_data
    html = ""

    # Stats cards
    html += f"""
    <div class="section">
        <h2>Expression Landscape</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="label">Total Genes</div>
                <div class="value">{ld['total_genes']:,}</div>
            </div>
            <div class="stat-card">
                <div class="label">Expressed (TPM &gt; 0)</div>
                <div class="value">{ld['expressed_genes']:,}</div>
                <div class="detail">{round(ld['expressed_genes']/ld['total_genes']*100)}% of total</div>
            </div>
            <div class="stat-card accent">
                <div class="label">Gini Coefficient</div>
                <div class="value">{ld['gini']:.3f}</div>
                <div class="detail">Expression concentration</div>
            </div>
            <div class="stat-card accent">
                <div class="label">Top 10 Genes</div>
                <div class="value">{ld['top10_pct']:.0f}%</div>
                <div class="detail">of total expression</div>
            </div>
            <div class="stat-card">
                <div class="label">Top 50 Genes</div>
                <div class="value">{ld['top50_pct']:.0f}%</div>
                <div class="detail">of total expression</div>
            </div>
            <div class="stat-card">
                <div class="label">Top 100 Genes</div>
                <div class="value">{ld['top100_pct']:.0f}%</div>
                <div class="detail">of total expression</div>
            </div>
        </div>
        <div class="dkl-interp">
            <strong>Gini = {ld['gini']:.3f}:</strong> {_gini_interpretation(ld['gini'])}
            Median TPM: {ld['median_tpm']:.1f} &middot; Mean TPM: {ld['mean_tpm']:.1f}
        </div>
    </div>"""

    # Lorenz curve (interactive)
    lorenz_traces = [
        {
            'x': ld['lorenz_x'], 'y': ld['lorenz_y'],
            'type': 'scatter', 'mode': 'lines',
            'name': f"Expression (Gini = {ld['gini']:.3f})",
            'line': {'color': '#0f2b46', 'width': 2.5},
            'fill': 'tozeroy', 'fillcolor': 'rgba(15,43,70,0.08)',
            'hovertemplate': ('Top %{customdata[0]:.1f}% of genes<br>'
                              'account for %{customdata[1]:.1f}% of expression<extra></extra>'),
            'customdata': [[round((1 - x) * 100, 1), round((1 - y) * 100, 1)]
                           for x, y in zip(ld['lorenz_x'], ld['lorenz_y'])],
        },
        {
            'x': [0, 1], 'y': [0, 1],
            'type': 'scatter', 'mode': 'lines',
            'name': 'Perfect equality',
            'line': {'color': '#999', 'dash': 'dash', 'width': 1},
            'hoverinfo': 'skip',
        },
    ]
    lorenz_layout = {
        'xaxis': {'title': 'Cumulative fraction of genes (ranked by TPM)',
                  'range': [0, 1], 'gridcolor': '#f0f0f0'},
        'yaxis': {'title': 'Cumulative fraction of total expression',
                  'range': [0, 1], 'gridcolor': '#f0f0f0'},
        'legend': {'x': 0.02, 'y': 0.98, 'font': {'size': 11}},
        'hovermode': 'closest',
        'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        'margin': {'l': 60, 'r': 20, 't': 20, 'b': 50},
        'annotations': [],
    }
    # Annotate key top-gene markers on Lorenz curve
    # Position annotations along right edge with correct expression percentages
    for i, (pct, label, expr_pct) in enumerate([
            (10, 'Top 10', ld['top10_pct']),
            (50, 'Top 50', ld['top50_pct']),
            (100, 'Top 100', ld['top100_pct'])]):
        if pct <= ld['total_genes']:
            # Place at y positions spaced down from top
            y_pos = 0.95 - i * 0.08
            lorenz_layout['annotations'].append({
                'x': 1.0, 'y': y_pos, 'xref': 'paper', 'yref': 'paper',
                'text': f'{label}: {expr_pct:.0f}%',
                'showarrow': False,
                'font': {'size': 11, 'color': '#0e7490', 'weight': 'bold'},
                'xanchor': 'right',
            })

    chart_id = 'landscape-lorenz'
    html += f"""
    <div class="section">
        <h2>Expression Concentration</h2>
        <p>The Lorenz curve shows how expression is distributed across genes.
        The further the curve bows from the diagonal, the more concentrated expression is
        among a few genes. Hover to explore.</p>
        <div id="{chart_id}" class="plot-container" style="min-height:400px;"></div>
        <div>
            <button class="export-btn" onclick="exportPlot('{chart_id}','lorenz_curve','svg')">Export SVG</button>
            <button class="export-btn" onclick="exportPlot('{chart_id}','lorenz_curve','png')">Export PNG</button>
        </div>
    </div>"""
    chart_render_calls.append(
        f"Plotly.newPlot('{chart_id}', {json.dumps(lorenz_traces)}, "
        f"{json.dumps(lorenz_layout)}, plotlyConfig);")

    # Expression histogram
    if ld['hist_counts']:
        bin_centres = [(ld['hist_edges'][i] + ld['hist_edges'][i+1]) / 2
                       for i in range(len(ld['hist_counts']))]
        bin_widths = [ld['hist_edges'][i+1] - ld['hist_edges'][i]
                      for i in range(len(ld['hist_counts']))]
        hist_traces = [{
            'x': bin_centres,
            'y': ld['hist_counts'],
            'type': 'bar',
            'width': bin_widths,
            'marker': {'color': '#0f2b46', 'opacity': 0.7,
                       'line': {'color': 'white', 'width': 0.5}},
            'name': 'Gene count',
            'hovertemplate': ('log₁₀(TPM) ≈ %{x:.1f}<br>'
                              '%{y} genes<extra></extra>'),
        }]
        # Median line
        median_log = np.log10(ld['median_tpm']) if ld['median_tpm'] > 0 else 0
        hist_layout = {
            'xaxis': {'title': 'log₁₀(TPM)', 'gridcolor': '#f0f0f0'},
            'yaxis': {'title': 'Number of genes', 'gridcolor': '#f0f0f0'},
            'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
            'margin': {'l': 60, 'r': 20, 't': 20, 'b': 50},
            'showlegend': False,
            'shapes': [{
                'type': 'line', 'x0': median_log, 'x1': median_log,
                'y0': 0, 'y1': 1, 'yref': 'paper',
                'line': {'color': '#0e7490', 'width': 2, 'dash': 'dash'},
            }],
            'annotations': [{
                'x': median_log, 'y': 1, 'yref': 'paper',
                'text': f'Median: {ld["median_tpm"]:.1f} TPM',
                'showarrow': False,
                'font': {'size': 11, 'color': '#0e7490'},
                'yanchor': 'bottom',
            }],
        }
        chart_id = 'landscape-histogram'
        html += f"""
        <div class="section">
            <h2>Expression Distribution</h2>
            <p>Distribution of expression values across all genes (log₁₀ TPM).
            The dashed line marks the median.</p>
            <div id="{chart_id}" class="plot-container" style="min-height:350px;"></div>
            <div>
                <button class="export-btn" onclick="exportPlot('{chart_id}','expression_distribution','svg')">Export SVG</button>
                <button class="export-btn" onclick="exportPlot('{chart_id}','expression_distribution','png')">Export PNG</button>
            </div>
        </div>"""
        chart_render_calls.append(
            f"Plotly.newPlot('{chart_id}', {json.dumps(hist_traces)}, "
            f"{json.dumps(hist_layout)}, plotlyConfig);")

    # Annotation quality by decile
    if ld['decile_stats']:
        deciles = [f"D{d['decile']}" for d in ld['decile_stats']]
        diamond_rates = [d['diamond_pct'] for d in ld['decile_stats']]
        go_rates = [d['go_pct'] for d in ld['decile_stats']]

        ann_traces = [
            {
                'x': deciles, 'y': diamond_rates,
                'type': 'bar', 'name': 'DIAMOND hit',
                'marker': {'color': '#0f2b46', 'opacity': 0.8},
                'hovertemplate': '%{x}<br>DIAMOND: %{y:.1f}%<extra></extra>',
            },
            {
                'x': deciles, 'y': go_rates,
                'type': 'bar', 'name': 'GO annotation',
                'marker': {'color': '#0ea5c9', 'opacity': 0.8},
                'hovertemplate': '%{x}<br>GO: %{y:.1f}%<extra></extra>',
            },
        ]
        ann_layout = {
            'barmode': 'group',
            'xaxis': {'title': 'Expression decile (1 = highest)', 'gridcolor': '#f0f0f0'},
            'yaxis': {'title': 'Annotation rate (%)', 'range': [0, 105], 'gridcolor': '#f0f0f0'},
            'legend': {'x': 0.02, 'y': 0.98, 'font': {'size': 11}},
            'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
            'margin': {'l': 60, 'r': 20, 't': 20, 'b': 50},
        }
        chart_id = 'landscape-annotation'
        html += f"""
        <div class="section">
            <h2>Annotation Quality</h2>
            <p>Annotation completeness across expression levels. Bias in annotation rate
            can affect enrichment results -highly expressed genes are typically better annotated.
            Overall: {ld['overall_diamond_pct']:.0f}% DIAMOND hits, {ld['overall_go_pct']:.0f}% GO annotated.</p>
            <div id="{chart_id}" class="plot-container" style="min-height:350px;"></div>
            <div>
                <button class="export-btn" onclick="exportPlot('{chart_id}','annotation_quality','svg')">Export SVG</button>
                <button class="export-btn" onclick="exportPlot('{chart_id}','annotation_quality','png')">Export PNG</button>
            </div>
        </div>"""
        chart_render_calls.append(
            f"Plotly.newPlot('{chart_id}', {json.dumps(ann_traces)}, "
            f"{json.dumps(ann_layout)}, plotlyConfig);")

    # Taxonomic distribution
    tax = ld.get('tax_data')
    if tax:
        names = list(tax['names'])
        counts = list(tax['counts'])
        total_hits = tax['total_hits']
        pcts = [round(c / total_hits * 100, 1) for c in counts]

        if tax['other_count'] > 0:
            names.append(f"Other ({tax['other_n_species']} species)")
            counts.append(tax['other_count'])
            pcts.append(round(tax['other_count'] / total_hits * 100, 1))

        # Reverse for horizontal bar (top at top)
        names_r = names[::-1]
        pcts_r = pcts[::-1]
        colors = ['#0f2b46'] * len(names_r)
        if tax['other_count'] > 0:
            colors[0] = '#bbb'  # "Other" bar

        tax_traces = [{
            'y': names_r, 'x': pcts_r,
            'type': 'bar', 'orientation': 'h',
            'marker': {'color': colors, 'opacity': 0.8},
            'text': [f'{p:.1f}%' if p > 2 else '' for p in pcts_r],
            'textposition': 'outside',
            'textfont': {'size': 10, 'color': '#546478'},
            'hovertemplate': '%{y}<br>%{x:.1f}% of hits<extra></extra>',
        }]
        tax_layout = {
            'xaxis': {'title': 'Percentage of DIAMOND hits', 'gridcolor': '#f0f0f0'},
            'yaxis': {'automargin': True},
            'showlegend': False,
            'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
            'margin': {'l': 250, 'r': 60, 't': 20, 'b': 50},
            'height': max(350, len(names_r) * 22),
        }
        chart_id = 'landscape-taxonomy'
        html += f"""
        <div class="section">
            <h2>Taxonomic Distribution</h2>
            <p>Species distribution of best DIAMOND hits reveals phylogenetic context -
            which reference organisms are closest to the study species and whose literature
            is most relevant for biological interpretation.
            {tax['unique_organisms']:,} unique species detected across {total_hits:,} annotated genes.</p>
            <div id="{chart_id}" class="plot-container" style="min-height:{max(350, len(names_r)*22)}px;"></div>
            <div>
                <button class="export-btn" onclick="exportPlot('{chart_id}','taxonomic_distribution','svg')">Export SVG</button>
                <button class="export-btn" onclick="exportPlot('{chart_id}','taxonomic_distribution','png')">Export PNG</button>
            </div>
        </div>"""
        chart_render_calls.append(
            f"Plotly.newPlot('{chart_id}', {json.dumps(tax_traces)}, "
            f"{json.dumps(tax_layout)}, plotlyConfig);")

    return html


# ---------------------------------------------------------------------------
# D_KL interpretation
# ---------------------------------------------------------------------------
def _generate_dkl_interpretation(k_values, dkl_values, shape_stats=None):
    if len(k_values) == 0 or len(dkl_values) == 0:
        return ""

    idx50 = np.argmin(np.abs(k_values - 50))
    idx500 = np.argmin(np.abs(k_values - 500))
    dkl_50 = dkl_values[idx50]
    dkl_500 = dkl_values[idx500]

    mask = dkl_values > 0
    slope = None
    r2 = None
    if np.sum(mask) > 10:
        from scipy.stats import linregress
        log_k = np.log(k_values[mask])
        log_dkl = np.log(dkl_values[mask])
        result = linregress(log_k, log_dkl)
        slope = result.slope
        r2 = result.rvalue ** 2

    lines = []
    lines.append(
        f"D<sub>KL</sub> = <strong>{dkl_50:.3f}</strong> at the top 50 genes, "
        f"declining to <strong>{dkl_500:.3f}</strong> at the top 500.")

    if dkl_50 > 0.5:
        lines.append(
            "<strong>High functional specialisation</strong> at the expression apex -"
            "the most highly expressed genes are dominated by a small number of "
            "functional programmes.")
    elif dkl_50 > 0.2:
        lines.append(
            "<strong>Moderate functional specialisation</strong> at the expression apex.")
    else:
        lines.append(
            "<strong>Relatively uniform functional distribution</strong> even among "
            "the most highly expressed genes.")

    if slope is not None:
        lines.append(
            f"Power-law decay: slope = {slope:.3f}, R² = {r2:.3f}.")
        if slope < -0.2:
            lines.append(
                "A steep slope -specialisation is <strong>concentrated at the "
                "expression apex</strong> and dilutes rapidly down the ranking.")
        elif slope < -0.1:
            lines.append(
                "A moderate slope -functional programmes are <strong>sustained "
                "across a broad expression range</strong>.")
        else:
            lines.append(
                "A shallow slope -enrichment is <strong>broadly distributed</strong> "
                "with no sharp boundary between functional programmes and background.")

    return " ".join(lines)


# ---------------------------------------------------------------------------
# Category report card
# ---------------------------------------------------------------------------
def _generate_category_card(cat, cat_idx, cat_names, cont_results,
                            results, all_results, tier_names, df_sorted=None):
    colour = _cat_colour(cat_names, cat_idx)

    tier_fcs = {}
    for tn in tier_names:
        stats = results.get(tn, {}).get('enrichment', {}).get(cat, {})
        tier_fcs[tn] = {
            'fc': stats.get('fold_change', 0),
            'p': stats.get('adjusted_p_value', 1.0),
            'sig': stats.get('significant', False),
        }

    shape_info = ""
    profile_p = ""
    if cont_results and 'shape_stats' in cont_results:
        ss = cont_results['shape_stats'].get(cat, {})
        shape_class = ss.get('shape_class', 'unknown')
        shape_icon = {'apex-concentrated': '⬆', 'distributed': '⬇',
                      'flat': '➡', 'absent': '-'}.get(shape_class, '?')
        shape_info = (f'<span class="shape-badge shape-{shape_class}">'
                      f'{shape_icon} {shape_class.replace("-", " ").title()}</span>')

    if cont_results and 'profile_stats' in cont_results:
        ps = cont_results['profile_stats'].get(cat, {})
        sup_p = ps.get('supremum_p', 1.0)
        if sup_p < 0.05:
            profile_p = f'<span class="sig-badge">p = {_format_pval(sup_p)}</span>'

    method_html = ""
    method_summary = all_results.get('method_summary', {}).get(cat, {})
    if method_summary:
        kw = method_summary.get('keyword_only', 0)
        go = method_summary.get('go_id_only', 0)
        both = method_summary.get('both', 0)
        total = kw + go + both
        if total > 0:
            kw_pct = kw / total * 100
            go_pct = go / total * 100
            both_pct = both / total * 100
            method_html = f"""
            <div class="method-bar">
                <div class="method-seg kw" style="width:{kw_pct}%"
                     title="Keyword: {kw}"></div>
                <div class="method-seg go" style="width:{go_pct}%"
                     title="GO: {go}"></div>
                <div class="method-seg both" style="width:{both_pct}%"
                     title="Both: {both}"></div>
            </div>
            <div class="method-legend">
                <span class="ml-kw">Keyword {kw_pct:.0f}%</span>
                <span class="ml-go">GO {go_pct:.0f}%</span>
                <span class="ml-both">Both {both_pct:.0f}%</span>
            </div>"""

    top_genes_html = ""
    gene_cats = all_results.get('gene_categories', {})
    if df_sorted is not None and gene_cats:
        cat_gene_indices = [int(gid) for gid, cats in gene_cats.items()
                           if cat in cats]
        if cat_gene_indices:
            cat_gene_indices.sort()
            top_n = min(8, len(cat_gene_indices))
            top_indices = cat_gene_indices[:top_n]
            gene_rows = []
            for idx in top_indices:
                if idx < len(df_sorted):
                    row = df_sorted.iloc[idx]
                    name = str(row.get('gene_names', row.get('protein_name', 'unknown')))
                    if len(name) > 50:
                        name = name[:47] + '...'
                    tpm = row.get('TPM', 0)
                    gene_rows.append(
                        f'<tr><td>#{idx+1}</td><td>{name}</td>'
                        f'<td>{tpm:,.0f}</td></tr>')
            if gene_rows:
                top_genes_html = f"""
                <table class="gene-table">
                    <thead><tr><th>Rank</th><th>Gene</th><th>TPM</th></tr></thead>
                    <tbody>{''.join(gene_rows)}</tbody>
                </table>"""

    interp = _build_category_interpretation(
        cat, cont_results, method_summary, tier_fcs, tier_names)

    return f"""
    <div class="category-card" style="border-left-color:{colour}">
        <div class="card-header">
            <h3 style="color:{colour}">{cat}</h3>
            <div class="card-badges">{shape_info} {profile_p}</div>
        </div>
        <p class="card-interp">{interp}</p>
        {method_html}
        {top_genes_html}
    </div>"""


def _build_category_interpretation(cat, cont_results, method_summary,
                                    tier_fcs, tier_names):
    parts = []
    if cont_results and 'shape_stats' in cont_results:
        ss = cont_results['shape_stats'].get(cat, {})
        max_fc = ss.get('max_enrichment', 0)
        shape_class = ss.get('shape_class', 'unknown')
        if max_fc > 1.5:
            if shape_class == 'apex-concentrated':
                parts.append(f"Strongly enriched at the expression apex (peak {max_fc:.1f}×) "
                             f"- concentrated among the most highly expressed genes.")
            elif shape_class == 'distributed':
                parts.append(f"Enrichment increases at broader tiers (peak {max_fc:.1f}×) "
                             f"- distributed across the expression gradient.")
            elif shape_class == 'flat':
                parts.append(f"Consistently enriched across tiers (peak {max_fc:.1f}×) "
                             f"- stable functional representation.")
        elif max_fc > 1.0:
            parts.append(f"Modestly enriched (peak {max_fc:.1f}×).")
        else:
            parts.append(f"Not enriched above background.")
    elif tier_fcs:
        max_tier = max(tier_fcs.items(), key=lambda x: x[1]['fc'])
        fc = max_tier[1]['fc']
        parts.append(f"Peak {fc:.1f}× enrichment at {max_tier[0]}.")
    return " ".join(parts) if parts else ""


# ---------------------------------------------------------------------------
# Plotly data builders
# ---------------------------------------------------------------------------
_MIN_GENES_CONTINUOUS = 25  # Suppress continuous curves for categories below this


def _build_enrichment_plotly_data(cont_results, colour_map, id_prefix='',
                                   cat_gene_counts=None):
    k_values = cont_results['k_values']
    enrichment_matrix = cont_results['enrichment_matrix']
    cat_names = cont_results['cat_names']
    profile_stats = cont_results.get('profile_stats', {})
    traces = []

    # Identify categories with enough genes for stable continuous curves
    if cat_gene_counts is None:
        cat_gene_counts = {}
    sparse_cats = {cat for cat in cat_names
                   if cat_gene_counts.get(cat, 9999) < _MIN_GENES_CONTINUOUS}
    if sparse_cats:
        logger.info(f"Suppressing continuous curves for {len(sparse_cats)} "
                     f"sparse categories (<{_MIN_GENES_CONTINUOUS} genes): "
                     f"{', '.join(sorted(sparse_cats))}")

    # Null envelope (only from non-sparse categories)
    all_lower, all_upper = [], []
    for cat in cat_names:
        if cat in profile_stats and cat not in sparse_cats:
            all_lower.append(profile_stats[cat]['null_envelope_lower'])
            all_upper.append(profile_stats[cat]['null_envelope_upper'])
    if all_lower:
        mean_lower = np.mean(all_lower, axis=0)
        mean_upper = np.mean(all_upper, axis=0)
        traces.append({
            'x': _np_to_list(k_values), 'y': _np_to_list(mean_upper),
            'type': 'scatter', 'mode': 'lines',
            'line': {'width': 0}, 'showlegend': False, 'hoverinfo': 'skip',
        })
        traces.append({
            'x': _np_to_list(k_values), 'y': _np_to_list(mean_lower),
            'type': 'scatter', 'mode': 'lines',
            'line': {'width': 0}, 'fill': 'tonexty',
            'fillcolor': 'rgba(180,180,180,0.15)',
            'name': '95% null envelope', 'hoverinfo': 'skip',
        })

    # Filter out sparse categories before ranking
    plottable = [(i, cat, np.max(enrichment_matrix[:, i]))
                 for i, cat in enumerate(cat_names) if cat not in sparse_cats]
    plottable.sort(key=lambda x: x[2], reverse=True)

    for rank, (i, cat, max_fc) in enumerate(plottable):
        colour = colour_map.get(cat, '#CCCCCC')
        p_info = ''
        if cat in profile_stats:
            sp = profile_stats[cat]['supremum_p']
            p_info = f', p={_format_pval(sp)}'
        visible = True if rank < 5 else 'legendonly'
        traces.append({
            'x': _np_to_list(k_values),
            'y': _np_to_list(enrichment_matrix[:, i]),
            'type': 'scatter', 'mode': 'lines', 'name': cat,
            'line': {'color': colour, 'width': 2.5 if rank < 5 else 1.2},
            'visible': visible,
            'hovertemplate': (f'<b>{cat}</b><br>k = %{{x:,.0f}}<br>'
                              f'E(k) = %{{y:.2f}}x{p_info}<extra></extra>'),
        })

    traces.append({
        'x': [_np_to_list(k_values)[0], _np_to_list(k_values)[-1]],
        'y': [1.0, 1.0],
        'type': 'scatter', 'mode': 'lines',
        'line': {'color': 'rgba(0,0,0,0.3)', 'dash': 'dash', 'width': 1},
        'showlegend': False, 'hoverinfo': 'skip',
    })

    layout = {
        'xaxis': {'title': 'Tier size k (genes)', 'type': 'log', 'gridcolor': '#f0f0f0'},
        'yaxis': {'title': 'Fold enrichment E<sub>C</sub>(k)', 'gridcolor': '#f0f0f0'},
        'legend': {'orientation': 'v', 'x': 1.02, 'y': 1, 'font': {'size': 11}},
        'hovermode': 'closest',
        'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        'margin': {'l': 60, 'r': 200, 't': 20, 'b': 50},
        'shapes': [{'type': 'line', 'x0': t, 'x1': t,
                    'y0': 0, 'y1': 1, 'yref': 'paper',
                    'line': {'color': '#ddd', 'width': 1, 'dash': 'dot'}}
                   for t in [50, 100, 250, 500]],
    }
    return traces, layout


def _build_dkl_plotly_data(cont_results):
    k_values = cont_results['k_values']
    dkl_values = cont_results['dkl_values']
    traces = [{
        'x': _np_to_list(k_values), 'y': _np_to_list(dkl_values),
        'type': 'scatter', 'mode': 'lines',
        'name': 'D<sub>KL</sub>',
        'line': {'color': '#3C5488', 'width': 2.5},
        'hovertemplate': 'k = %{x:,.0f}<br>D<sub>KL</sub> = %{y:.4f}<extra></extra>',
    }]
    mask = dkl_values > 0
    if np.sum(mask) > 10:
        from scipy.stats import linregress
        log_k = np.log(k_values[mask])
        log_dkl = np.log(dkl_values[mask])
        result = linregress(log_k, log_dkl)
        fit_dkl = np.exp(result.intercept + result.slope * log_k)
        traces.append({
            'x': _np_to_list(k_values[mask]), 'y': _np_to_list(fit_dkl),
            'type': 'scatter', 'mode': 'lines',
            'name': f'Power-law fit (\u03b2={result.slope:.3f}, R\u00b2={result.rvalue**2:.3f})',
            'line': {'color': '#E64B35', 'width': 1.5, 'dash': 'dash'},
            'hoverinfo': 'skip',
        })
    layout = {
        'xaxis': {'title': 'Tier size k', 'type': 'log', 'gridcolor': '#f0f0f0'},
        'yaxis': {'title': 'D<sub>KL</sub> (bits)', 'type': 'log', 'gridcolor': '#f0f0f0'},
        'hovermode': 'closest',
        'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        'legend': {'font': {'size': 10}},
        'margin': {'l': 60, 'r': 40, 't': 20, 'b': 50},
    }
    return traces, layout


def _build_bar_plotly_data(results, tier_name, colour_map):
    enrichment = results.get(tier_name, {}).get('enrichment', {})
    if not enrichment:
        return [], {}
    cats = sorted(enrichment.keys(),
                  key=lambda c: enrichment[c].get('fold_change', 0), reverse=True)
    cats = [c for c in cats if 'Uncharacterised' not in c and 'uncharacterised' not in c]
    fcs = [enrichment[c].get('fold_change', 0) for c in cats]
    sigs = [enrichment[c].get('significant', False) for c in cats]
    pvals = [enrichment[c].get('adjusted_p_value', 1.0) for c in cats]
    colours = [colour_map.get(c, '#CCCCCC') for c in cats]
    opacities = [0.9 if s else 0.3 for s in sigs]
    trace = {
        'y': cats, 'x': fcs, 'type': 'bar', 'orientation': 'h',
        'marker': {'color': colours, 'opacity': opacities,
                   'line': {'width': 0.5, 'color': 'white'}},
        'hovertemplate': [
            f'<b>{c}</b><br>FC = {fc:.2f}×<br>FDR = {_format_pval(p)}'
            f'{"<br><b>Significant</b>" if s else ""}<extra></extra>'
            for c, fc, p, s in zip(cats, fcs, pvals, sigs)
        ],
    }
    tier_label = tier_name.replace('top_', 'Top ')
    layout = {
        'title': {'text': f'{tier_label} genes', 'font': {'size': 15}},
        'xaxis': {'title': 'Fold enrichment', 'gridcolor': '#f0f0f0'},
        'yaxis': {'autorange': 'reversed', 'tickfont': {'size': 11}},
        'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        'margin': {'l': 250, 'r': 30, 't': 50, 'b': 50},
        'showlegend': False,
        'shapes': [{'type': 'line', 'x0': 1, 'x1': 1,
                    'y0': -0.5, 'y1': len(cats) - 0.5,
                    'line': {'color': '#aaa', 'width': 1, 'dash': 'dash'}}],
    }
    return [trace], layout


def _build_significance_landscape_data(cont_results, colour_map,
                                        cat_gene_counts=None):
    k_values = cont_results['k_values']
    enrichment_matrix = cont_results['enrichment_matrix']
    cat_names = cont_results['cat_names']
    profile_stats = cont_results.get('profile_stats', {})
    if cat_gene_counts is None:
        cat_gene_counts = {}
    sparse_cats = {cat for cat in cat_names
                   if cat_gene_counts.get(cat, 9999) < _MIN_GENES_CONTINUOUS}
    integrals = [profile_stats[cat]['integral_obs'] if cat in profile_stats else 0
                 for cat in cat_names]
    sort_idx = np.argsort(integrals)[::-1]
    sorted_cats = [cat_names[i] for i in sort_idx]
    traces = []
    for cat in sorted_cats:
        if cat in sparse_cats:
            continue
        i = cat_names.index(cat)
        colour = colour_map.get(cat, '#CCCCCC')
        if cat not in profile_stats:
            continue
        env_lo = profile_stats[cat]['null_envelope_lower']
        env_hi = profile_stats[cat]['null_envelope_upper']
        curve = enrichment_matrix[:, i]
        outside = (curve > env_hi) | (curve < env_lo)
        if not np.any(outside):
            continue
        sig_k = k_values[outside]
        if len(sig_k) == 0:
            continue
        traces.append({
            'x': [float(sig_k[0]), float(sig_k[-1])], 'y': [cat, cat],
            'type': 'scatter', 'mode': 'lines',
            'line': {'color': colour, 'width': 12}, 'showlegend': False,
            'hovertemplate': (f'<b>{cat}</b><br>Significant from k={sig_k[0]:.0f} '
                              f'to k={sig_k[-1]:.0f}<extra></extra>'),
        })
    layout = {
        'xaxis': {'title': 'Tier size k', 'type': 'log', 'gridcolor': '#f0f0f0'},
        'yaxis': {'tickfont': {'size': 10}},
        'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        'margin': {'l': 250, 'r': 30, 't': 20, 'b': 50},
        'height': max(250, len(traces) * 28 + 100),
    }
    return traces, layout


def _build_enrichment_table_html(results, tier_names, cont_results=None):
    cats = []
    if tier_names and tier_names[0] in results:
        cats = list(results[tier_names[0]].get('enrichment', {}).keys())
    cats = [c for c in cats if 'Uncharacterised' not in c]
    if not cats:
        return "<p>No enrichment data available.</p>"

    shape_stats = cont_results.get('shape_stats', {}) if cont_results else {}
    trend_icons = {
        'apex-concentrated': '↘ Declining', 'distributed': '↗ Rising',
        'flat': '→ Stable', 'absent': '- Absent',
    }

    header = '<th>Category</th>'
    for tn in tier_names:
        header += f'<th>{tn.replace("top_", "Top ")}<br><small>FC (FDR)</small></th>'
    if shape_stats:
        header += '<th>Profile Trend</th>'

    rows = []
    for cat in cats:
        cells = f'<td><strong>{cat}</strong></td>'
        any_sig = False
        for tn in tier_names:
            stats = results.get(tn, {}).get('enrichment', {}).get(cat, {})
            fc = stats.get('fold_change', 0)
            adj_p = stats.get('adjusted_p_value', stats.get('p_value', 1.0))
            sig = stats.get('significant', False)
            if sig:
                any_sig = True
            if fc == 0:
                cells += '<td class="na-cell">-</td>'
            else:
                cls = 'enriched' if fc > 1.5 and sig else (
                    'depleted' if fc < 0.67 and sig else '')
                sig_mark = ' *' if sig else ''
                cells += (f'<td class="{cls}">{fc:.2f}×{sig_mark}<br>'
                          f'<small>{_format_pval(adj_p)}</small></td>')
        if shape_stats:
            ss = shape_stats.get(cat, {})
            shape_class = ss.get('shape_class', 'absent')
            cells += f'<td class="trend-cell trend-{shape_class}">{trend_icons.get(shape_class, "?")}</td>'
        row_cls = ' class="sig"' if any_sig else ''
        rows.append(f'<tr{row_cls}>{cells}</tr>')

    return f"""
    <table class="enrichment-table">
        <thead><tr>{header}</tr></thead>
        <tbody>{''.join(rows)}</tbody>
    </table>
    <p class="table-note">* Significant (BH-adjusted p &lt; 0.05). FC = fold change.
    ↘ = apex-concentrated, ↗ = distributed, → = stable across tiers.</p>"""


# ---------------------------------------------------------------------------
# Build all sections for one analysis type
# ---------------------------------------------------------------------------
def _build_analysis_sections(block, chart_render_calls, chart_counter):
    """Build HTML sections for one analysis block (functional or cellular).
    Returns (sections_html, chart_counter)."""
    results = block['results']
    all_results = block['all_results']
    cont_results = block.get('cont_results')
    tier_names = block['tier_names']
    cat_names = block['cat_names']
    colour_map = block['colour_map']
    df_sorted = block.get('df_sorted')
    section_title = block['section_title']
    section_desc = block['section_desc']
    prefix = block['id_prefix']

    html = ""

    # Count genes per category for sparse-category filtering
    gene_cats = all_results.get('gene_categories', {})
    cat_gene_counts = {}
    for gid, cats in gene_cats.items():
        for c in cats:
            cat_gene_counts[c] = cat_gene_counts.get(c, 0) + 1

    # Section header
    html += f"""
    <div class="section-divider">
        <h2 class="analysis-title">{section_title}</h2>
        <p class="analysis-desc">{section_desc}</p>
    </div>"""

    # Enrichment curves
    if cont_results:
        traces, layout = _build_enrichment_plotly_data(
            cont_results, colour_map, cat_gene_counts=cat_gene_counts)
        chart_id = f'{prefix}-enrichment'
        html += f"""
        <div class="section">
            <h2>Expression Profile</h2>
            <p>How each category is distributed across your expression ranking.
            Curves above the dashed line are enriched; the grey band is the 95% null envelope.
            Click legend entries to show/hide. Hover for values.</p>
            <div id="{chart_id}" class="plot-container"></div>
            <div>
                <button class="export-btn" onclick="exportPlot('{chart_id}','{prefix}_enrichment_curves','svg')">Export SVG</button>
                <button class="export-btn" onclick="exportPlot('{chart_id}','{prefix}_enrichment_curves','png')">Export PNG</button>
            </div>
        </div>"""
        chart_render_calls.append(
            f"Plotly.newPlot('{chart_id}', {json.dumps(traces, default=str)}, "
            f"{json.dumps(layout)}, plotlyConfig);")

    # Bar chart
    default_tier = 'top_250' if 'top_250' in tier_names else tier_names[-1]
    bar_traces, bar_layout = _build_bar_plotly_data(results, default_tier, colour_map)
    if bar_traces:
        chart_id = f'{prefix}-bar'
        html += f"""
        <div class="section">
            <h2>What's Enriched?</h2>
            <p>Fold enrichment at the {default_tier.replace('top_', 'top ')} most
            highly expressed genes. Bold = significant (FDR &lt; 0.05); faded = not.</p>
            <div id="{chart_id}" class="plot-container" style="min-height:350px;"></div>
            <div>
                <button class="export-btn" onclick="exportPlot('{chart_id}','{prefix}_enrichment_bar','svg')">Export SVG</button>
                <button class="export-btn" onclick="exportPlot('{chart_id}','{prefix}_enrichment_bar','png')">Export PNG</button>
            </div>
        </div>"""
        chart_render_calls.append(
            f"Plotly.newPlot('{chart_id}', {json.dumps(bar_traces, default=str)}, "
            f"{json.dumps(bar_layout)}, plotlyConfig);")

    # Category cards
    if cont_results:
        enriched_cats = []
        for cat in cat_names:
            ss = cont_results.get('shape_stats', {}).get(cat, {})
            max_fc = ss.get('max_enrichment', 0)
            ps = cont_results.get('profile_stats', {}).get(cat, {})
            sup_p = ps.get('supremum_p', 1.0)
            int_p = ps.get('integral_p', 1.0)
            if max_fc > 1.2 or sup_p < 0.05 or int_p < 0.05:
                enriched_cats.append((cat, max_fc))
        enriched_cats.sort(key=lambda x: x[1], reverse=True)
        if enriched_cats:
            card_blocks = []
            for cat, _ in enriched_cats:
                cat_idx = cat_names.index(cat)
                card_blocks.append(_generate_category_card(
                    cat, cat_idx, cat_names, cont_results,
                    results, all_results, tier_names, df_sorted))
            html += f"""
            <div class="section">
                <h2>Category Breakdown</h2>
                <p>Each enriched category with its profile shape, top genes, and how
                genes were assigned.</p>
                <div class="category-cards">{''.join(card_blocks)}</div>
            </div>"""

    # D_KL
    if cont_results and 'dkl_values' in cont_results:
        dkl_traces, dkl_layout = _build_dkl_plotly_data(cont_results)
        dkl_interp = _generate_dkl_interpretation(
            cont_results['k_values'], cont_results['dkl_values'],
            cont_results.get('shape_stats'))
        chart_id = f'{prefix}-dkl'
        html += f"""
        <div class="section">
            <h2>Specialisation Gradient</h2>
            <div class="dkl-interp">{dkl_interp}</div>
            <div id="{chart_id}" class="plot-container" style="min-height:350px;"></div>
            <div>
                <button class="export-btn" onclick="exportPlot('{chart_id}','{prefix}_dkl','svg')">Export SVG</button>
                <button class="export-btn" onclick="exportPlot('{chart_id}','{prefix}_dkl','png')">Export PNG</button>
            </div>
        </div>"""
        chart_render_calls.append(
            f"Plotly.newPlot('{chart_id}', {json.dumps(dkl_traces, default=str)}, "
            f"{json.dumps(dkl_layout)}, plotlyConfig);")

    # Significance landscape
    if cont_results and 'profile_stats' in cont_results:
        sig_traces, sig_layout = _build_significance_landscape_data(
            cont_results, colour_map, cat_gene_counts=cat_gene_counts)
        if sig_traces:
            chart_id = f'{prefix}-sig'
            html += f"""
            <div class="section">
                <h2>Significance Landscape</h2>
                <p>Where along the expression ranking each category's enrichment
                exceeds the 95% null envelope.</p>
                <div id="{chart_id}" class="plot-container"></div>
            </div>"""
            chart_render_calls.append(
                f"Plotly.newPlot('{chart_id}', {json.dumps(sig_traces, default=str)}, "
                f"{json.dumps(sig_layout)}, plotlyConfig);")

    # Enrichment table
    table_html = _build_enrichment_table_html(results, tier_names, cont_results)
    html += f"""
    <div class="section">
        <h2>Enrichment Summary</h2>
        <p>Fold enrichment at each tier (Fisher's exact test, BH-adjusted)
        {" with continuous profile trend" if cont_results else ""}.</p>
        {table_html}
    </div>"""

    return html


# ---------------------------------------------------------------------------
# CSS
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
    --enriched: #0ea5c9;
    --depleted: #e74c3c;
    --sig-bg: #f0fafb;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    line-height: 1.6; color: var(--text);
    max-width: 1400px; margin: 0 auto; padding: 20px;
    background: var(--bg);
}
.header {
    background: linear-gradient(135deg, var(--primary), var(--primary-dark));
    color: white; padding: 32px 36px; border-radius: 12px;
    margin-bottom: 24px; position: relative; overflow: hidden;
}
.header::after {
    content: ''; position: absolute; top: -50%; right: -10%;
    width: 350px; height: 350px; border-radius: 50%;
    background: rgba(255,255,255,0.05);
}
.header h1 { font-size: 1.8em; font-weight: 700; margin-bottom: 4px; position: relative; z-index: 1; }
.header .subtitle { opacity: 0.85; font-size: 1.0em; position: relative; z-index: 1; }
.header .meta { opacity: 0.6; font-size: 0.82em; margin-top: 8px; position: relative; z-index: 1; }

/* --- Tab navigation --- */
.tab-bar {
    display: flex; gap: 0;
    background: var(--card);
    border-radius: 12px 12px 0 0;
    box-shadow: 0 1px 4px rgba(0,0,0,0.06);
    margin-bottom: 0;
    position: sticky; top: 0; z-index: 100;
    overflow: hidden;
}
.tab-btn {
    flex: 1;
    padding: 18px 24px;
    background: transparent;
    border: none;
    cursor: pointer;
    font-family: inherit;
    font-size: 0.95em;
    font-weight: 600;
    color: var(--text-light);
    position: relative;
    transition: color 0.2s, background 0.2s;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 4px;
}
.tab-btn:hover {
    background: rgba(15,43,70,0.04);
    color: var(--text);
}
.tab-btn.active {
    color: var(--primary);
    background: rgba(15,43,70,0.06);
}
.tab-btn::after {
    content: '';
    position: absolute;
    bottom: 0; left: 20%; right: 20%;
    height: 3px;
    background: transparent;
    border-radius: 3px 3px 0 0;
    transition: background 0.25s, left 0.25s, right 0.25s;
}
.tab-btn.active::after {
    background: var(--accent);
    left: 10%; right: 10%;
}
.tab-icon {
    font-size: 1.4em;
    line-height: 1;
}
.tab-label {
    font-size: 0.82em;
    letter-spacing: 0.3px;
    text-transform: uppercase;
}
.tab-hint {
    font-size: 0.68em;
    font-weight: 400;
    color: var(--text-light);
    opacity: 0.7;
}
.tab-btn.active .tab-hint { opacity: 1; color: var(--accent-dim); }
.tab-panel {
    display: none;
    animation: tabFadeIn 0.3s ease;
}
.tab-panel.active { display: block; }
@keyframes tabFadeIn {
    from { opacity: 0; transform: translateY(6px); }
    to { opacity: 1; transform: translateY(0); }
}
.tab-content-area {
    background: var(--card);
    border-radius: 0 0 12px 12px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.06);
    padding: 4px 0 0 0;
    margin-bottom: 24px;
}

.section {
    background: var(--card); padding: 24px 28px;
    margin-bottom: 20px; border-radius: 10px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.06);
}
.section h2 {
    color: var(--primary); font-size: 1.3em; font-weight: 700;
    border-bottom: 2px solid var(--primary-light);
    padding-bottom: 8px; margin-bottom: 16px;
}
.section p { margin-bottom: 10px; font-size: 0.92em; }
.hero-summary {
    background: linear-gradient(135deg, #f8fafb, #eef3f7);
    padding: 20px 24px; border-radius: 10px;
    border-left: 4px solid var(--primary);
    font-size: 1.05em; line-height: 1.8;
    margin-bottom: 20px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.06);
}
.hero-summary strong { color: var(--primary-dark); }
.section-divider {
    background: linear-gradient(135deg, var(--primary-light), var(--primary));
    color: white; padding: 18px 28px; border-radius: 10px;
    margin: 32px 0 20px 0;
}
.section-divider .analysis-title {
    font-size: 1.4em; font-weight: 700; margin: 0;
    border: none; padding: 0; color: white;
}
.section-divider .analysis-desc {
    opacity: 0.85; font-size: 0.92em; margin: 4px 0 0 0;
}
.stats-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 14px; margin: 14px 0 6px 0;
}
.stat-card {
    background: var(--bg); padding: 16px 14px;
    border-radius: 8px; text-align: center;
    border-left: 3px solid var(--primary);
}
.stat-card.accent { border-left-color: var(--accent); }
.stat-card .label {
    font-size: 0.72em; text-transform: uppercase;
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
.enrichment-table {
    width: 100%; border-collapse: collapse; margin: 14px 0; font-size: 0.85em;
}
.enrichment-table th, .enrichment-table td {
    padding: 8px 12px; text-align: left; border-bottom: 1px solid var(--border);
}
.enrichment-table th {
    background: var(--primary); color: white;
    font-weight: 600; font-size: 0.82em;
    text-transform: uppercase; letter-spacing: 0.3px;
}
.enrichment-table tr:hover { background: #f8fafb; }
.enrichment-table tr.sig { background: var(--sig-bg); }
.enriched { color: var(--enriched); font-weight: 700; }
.depleted { color: var(--depleted); font-weight: 700; }
.na-cell { color: var(--text-light); }
.table-note { font-size: 0.8em; color: var(--text-light); margin-top: 8px; }
.trend-cell { font-weight: 600; white-space: nowrap; }
.trend-apex-concentrated { color: #E64B35; }
.trend-distributed { color: #4DBBD5; }
.trend-flat { color: #7F7F7F; }
.trend-absent { color: #CCCCCC; }
.dkl-interp {
    background: #f8fafb; padding: 16px 20px;
    border-radius: 8px; border-left: 3px solid var(--accent);
    font-size: 0.9em; line-height: 1.7; color: var(--text-light);
    margin: 16px 0;
}
.dkl-interp strong { color: var(--text); }
.category-cards {
    display: grid; grid-template-columns: repeat(auto-fill, minmax(380px, 1fr));
    gap: 16px; margin: 16px 0;
}
.category-card {
    background: var(--bg); padding: 16px 18px;
    border-radius: 8px; border-left: 4px solid #ccc;
}
.card-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px; }
.card-header h3 { font-size: 1.0em; margin: 0; }
.card-badges { display: flex; gap: 6px; flex-wrap: wrap; }
.shape-badge { font-size: 0.72em; padding: 2px 8px; border-radius: 12px; background: #e8e8e8; color: var(--text); }
.shape-apex-concentrated { background: #fde8e4; color: #E64B35; }
.shape-distributed { background: #e4f4f8; color: #2980b9; }
.shape-flat { background: #f0f0f0; color: #666; }
.sig-badge { font-size: 0.72em; padding: 2px 8px; border-radius: 12px; background: #e0f7fa; color: #0e7490; font-weight: 600; }
.card-interp { font-size: 0.85em; color: var(--text-light); margin-bottom: 10px; }
.method-bar { display: flex; height: 8px; border-radius: 4px; overflow: hidden; margin: 6px 0 4px 0; }
.method-seg { height: 100%; }
.method-seg.kw { background: #3C5488; }
.method-seg.go { background: #00A087; }
.method-seg.both { background: #F39B7F; }
.method-legend { font-size: 0.72em; color: var(--text-light); display: flex; gap: 10px; }
.ml-kw::before { content: '●'; color: #3C5488; margin-right: 2px; }
.ml-go::before { content: '●'; color: #00A087; margin-right: 2px; }
.ml-both::before { content: '●'; color: #F39B7F; margin-right: 2px; }
.gene-table { width: 100%; border-collapse: collapse; margin: 8px 0; font-size: 0.78em; }
.gene-table th, .gene-table td { padding: 3px 8px; text-align: left; border-bottom: 1px solid #eee; }
.gene-table th { background: var(--primary-light); color: white; font-size: 0.82em; font-weight: 600; }
.dominance-table { width: 100%; border-collapse: collapse; margin: 14px 0; font-size: 0.85em; }
.dominance-table th, .dominance-table td {
    padding: 8px 12px; text-align: left; border-bottom: 1px solid var(--border);
}
.dominance-table th {
    background: var(--primary); color: white;
    font-weight: 600; font-size: 0.82em;
}
.dominance-table tr:hover { background: #f8fafb; }
.methods-box {
    background: #f8fafb; padding: 18px 22px; border-radius: 8px;
    border-left: 3px solid var(--accent);
    font-size: 0.88em; line-height: 1.7; color: var(--text-light);
}
.methods-box p { margin-bottom: 8px; }
.methods-box strong { color: var(--text); }
.footer { text-align: center; padding: 16px; color: var(--text-light); font-size: 0.78em; }
@media print {
    .export-btn { display: none; }
    .tab-bar { display: none; }
    .tab-panel { display: block !important; }
    .plot-container { page-break-inside: avoid; }
}
@media (max-width: 850px) {
    .category-cards { grid-template-columns: 1fr; }
    .tab-btn { padding: 12px 8px; }
    .tab-label { font-size: 0.72em; }
    .tab-hint { display: none; }
}
"""

_JS_TEMPLATE = """
<script>
document.addEventListener('DOMContentLoaded', function() {
    var plotlyConfig = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['lasso2d', 'select2d'],
        toImageButtonOptions: { format: 'svg', filename: 'sceptr_figure', scale: 2 }
    };

    // Tab switching
    var tabBtns = document.querySelectorAll('.tab-btn');
    var tabPanels = document.querySelectorAll('.tab-panel');
    var pendingRender = {};

    TAB_CHART_MAP

    function switchTab(tabId) {
        tabBtns.forEach(function(b) { b.classList.remove('active'); });
        tabPanels.forEach(function(p) { p.classList.remove('active'); });
        var btn = document.querySelector('[data-tab="' + tabId + '"]');
        var panel = document.getElementById(tabId);
        if (btn) btn.classList.add('active');
        if (panel) panel.classList.add('active');
        // Render charts for this tab on first visit
        if (pendingRender[tabId]) {
            pendingRender[tabId].forEach(function(fn) { fn(); });
            delete pendingRender[tabId];
            // Relayout after render so Plotly sizes correctly
            setTimeout(function() {
                panel.querySelectorAll('.plot-container > div').forEach(function(el) {
                    if (el.id) Plotly.Plots.resize(el);
                });
            }, 50);
        }
    }

    tabBtns.forEach(function(btn) {
        btn.addEventListener('click', function() {
            switchTab(this.getAttribute('data-tab'));
        });
    });

    // Register chart render functions by tab
    CHART_RENDER_CALLS

    // Activate first tab and render its charts
    if (tabBtns.length > 0) {
        switchTab(tabBtns[0].getAttribute('data-tab'));
    }
});
function exportPlot(divId, filename, format) {
    var opts = {format: format, width: 1200, height: 700, scale: 3};
    if (format === 'svg') { opts.scale = 1; }
    Plotly.downloadImage(divId, Object.assign(opts, {filename: filename}));
}
</script>
"""


# ---------------------------------------------------------------------------
# Dominance table builder (for landscape tab)
# ---------------------------------------------------------------------------
def _build_dominance_table(df_sorted, top_n=30):
    """Build HTML table showing the most highly expressed genes."""
    if df_sorted is None or len(df_sorted) == 0:
        return ""
    tpm = df_sorted['TPM'].values
    total_expr = tpm.sum()
    if total_expr == 0:
        return ""

    top = df_sorted.head(top_n).copy()
    rows = []
    cum_pct = 0.0
    for i, (_, row) in enumerate(top.iterrows()):
        tpm_val = row.get('TPM', 0)
        pct = tpm_val / total_expr * 100
        cum_pct += pct
        name = str(row.get('protein_name', row.get('gene_names', 'unknown')))
        if len(name) > 55:
            name = name[:52] + '...'
        rows.append(
            f'<tr><td>{i+1}</td><td>{name}</td>'
            f'<td><strong>{tpm_val:,.1f}</strong></td>'
            f'<td>{pct:.1f}%</td><td>{cum_pct:.1f}%</td></tr>')

    return f"""
    <div class="section">
        <h2>Transcriptome Dominance</h2>
        <p>The {min(top_n, len(top))} most highly expressed genes and their
        cumulative contribution to total expression.</p>
        <table class="dominance-table">
            <thead><tr><th>Rank</th><th>Protein</th><th>TPM</th>
            <th>% Total</th><th>Cumul. %</th></tr></thead>
            <tbody>{''.join(rows)}</tbody>
        </table>
    </div>"""


# ---------------------------------------------------------------------------
# Main report generator - supports single or combined analysis
# ---------------------------------------------------------------------------
def generate_interactive_report(
    results: Dict[str, Dict[str, Any]],
    all_results: Dict[str, Any],
    output_prefix: str,
    total_genes: int,
    cont_results: Optional[Dict] = None,
    chart_type: str = 'BP_MF',
    report_title: str = 'Functional Profiling',
    description: str = 'biological process and molecular function',
    figures_dir: Optional[str] = None,
    df_sorted=None,
) -> Optional[str]:
    """Generate interactive report for a single analysis type."""
    block = _make_block(results, all_results, cont_results, chart_type,
                        report_title, description, df_sorted)
    return _generate_report(
        [block], output_prefix, total_genes,
        f'{report_title} Report', chart_type)


def generate_combined_report(
    functional_data: Dict,
    cellular_data: Dict,
    output_path: str,
    total_genes: int,
    df_sorted=None,
    landscape_data: Optional[Dict] = None,
) -> Optional[str]:
    """Generate a unified tabbed report combining functional + cellular + landscape."""
    blocks = []

    func_block = _make_block(
        functional_data['results'], functional_data['all_results'],
        functional_data.get('cont_results'),
        'func', 'What Your Transcriptome Is Doing',
        'Functional categories: biological process and molecular function',
        df_sorted)
    blocks.append(func_block)

    cell_block = _make_block(
        cellular_data['results'], cellular_data['all_results'],
        cellular_data.get('cont_results'),
        'cell', 'Where It\'s Doing It',
        'Cellular component localisation',
        df_sorted)
    blocks.append(cell_block)

    return _generate_report(
        blocks, output_path, total_genes,
        'Expression Profile Report', 'combined',
        is_path=True, landscape_data=landscape_data)


def _make_block(results, all_results, cont_results, id_prefix,
                section_title, section_desc, df_sorted=None):
    """Create an analysis block dict."""
    tier_names = sorted(results.keys(), key=lambda x: int(x.split('_')[-1]))
    cat_names = cont_results['cat_names'] if cont_results else list(
        results.get(tier_names[0], {}).get('enrichment', {}).keys())
    cat_names = [c for c in cat_names if 'Uncharacterised' not in c]
    colour_map = _cat_colour_map(cat_names)

    return {
        'results': results,
        'all_results': all_results,
        'cont_results': cont_results,
        'tier_names': tier_names,
        'cat_names': cat_names,
        'colour_map': colour_map,
        'df_sorted': df_sorted,
        'id_prefix': id_prefix,
        'section_title': section_title,
        'section_desc': section_desc,
    }


def _generate_report(analysis_blocks, output_prefix, total_genes,
                     report_title, chart_type, is_path=False,
                     landscape_data=None):
    """Core report generator supporting one or multiple analysis blocks with tabs."""
    html_path = output_prefix if is_path else f"{output_prefix}_{chart_type}_report.html"

    # Aggregate stats across blocks
    all_cat_names = []
    for b in analysis_blocks:
        all_cat_names.extend(b['cat_names'])

    # Use first block's stats for the summary cards (functional is primary)
    primary = analysis_blocks[0]
    p_all = primary['all_results']
    pct_cat = round(p_all.get('total_categorised', 0) / total_genes * 100, 1) if total_genes > 0 else 0

    sig_count = 0
    profile_sig = 0
    for b in analysis_blocks:
        for tn in b['tier_names']:
            for stats in b['results'].get(tn, {}).get('enrichment', {}).values():
                if stats.get('significant', False):
                    sig_count += 1
        cr = b.get('cont_results')
        if cr and 'profile_stats' in cr:
            profile_sig += sum(1 for s in cr['profile_stats'].values()
                               if s.get('supremum_p', 1) < 0.05)

    has_continuous = any(b.get('cont_results') for b in analysis_blocks)
    n_tiers = len(analysis_blocks[0]['tier_names'])

    # Hero summary
    hero_text = _generate_hero_summary(analysis_blocks)

    # Compute landscape data from expression dataframe if not provided
    df_sorted = primary.get('df_sorted')
    if landscape_data is None and df_sorted is not None and len(df_sorted) > 0:
        landscape_data = _compute_landscape_data(df_sorted)

    # Determine which tabs to show
    # tab_id -> (icon, label, hint)
    # SVG icons for tabs (inline, 20x20)
    _icon_func = ('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" '
                  'viewBox="0 0 24 24" fill="none" stroke="currentColor" '
                  'stroke-width="2" stroke-linecap="round" stroke-linejoin="round">'
                  '<polyline points="22 12 18 12 15 21 9 3 6 12 2 12"/></svg>')
    _icon_cell = ('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" '
                  'viewBox="0 0 24 24" fill="none" stroke="currentColor" '
                  'stroke-width="2" stroke-linecap="round" stroke-linejoin="round">'
                  '<circle cx="12" cy="12" r="10"/>'
                  '<circle cx="12" cy="12" r="4"/>'
                  '<line x1="12" y1="2" x2="12" y2="8"/>'
                  '<line x1="12" y1="16" x2="12" y2="22"/>'
                  '<line x1="2" y1="12" x2="8" y2="12"/>'
                  '<line x1="16" y1="12" x2="22" y2="12"/></svg>')
    _icon_land = ('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" '
                  'viewBox="0 0 24 24" fill="none" stroke="currentColor" '
                  'stroke-width="2" stroke-linecap="round" stroke-linejoin="round">'
                  '<rect x="3" y="3" width="18" height="18" rx="2"/>'
                  '<line x1="3" y1="9" x2="21" y2="9"/>'
                  '<line x1="9" y1="21" x2="9" y2="9"/></svg>')

    tabs = []
    tabs.append(('tab-functional', _icon_func,
                 'Functional Profiling', 'Biological process & molecular function'))
    if len(analysis_blocks) > 1:
        tabs.append(('tab-cellular', _icon_cell,
                     'Localisation Profiling', 'Cellular component'))
    if landscape_data is not None:
        tabs.append(('tab-landscape', _icon_land,
                     'Transcriptome Landscape', 'Expression, annotation & taxonomy'))

    # Build tab bar HTML
    tab_bar_html = '<div class="tab-bar">\n'
    for tab_id, icon, label, hint in tabs:
        tab_bar_html += (
            f'    <button class="tab-btn" data-tab="{tab_id}">'
            f'<span class="tab-icon">{icon}</span>'
            f'<span class="tab-label">{label}</span>'
            f'<span class="tab-hint">{hint}</span>'
            f'</button>\n')
    tab_bar_html += '</div>\n'

    # Build each tab panel's content
    # chart_render_calls is now list of (tab_id, js_code)
    chart_render_calls = []

    # Functional tab
    func_calls = []
    func_html = _build_analysis_sections(
        analysis_blocks[0], func_calls, 0)
    for js in func_calls:
        chart_render_calls.append(('tab-functional', js))

    # Cellular tab (if present)
    cell_html = ""
    if len(analysis_blocks) > 1:
        cell_calls = []
        cell_html = _build_analysis_sections(
            analysis_blocks[1], cell_calls, 1)
        for js in cell_calls:
            chart_render_calls.append(('tab-cellular', js))

    # Landscape tab (if data available)
    land_html = ""
    if landscape_data is not None:
        land_calls = []
        land_html = _build_landscape_section(landscape_data, land_calls)
        # Add dominance table
        land_html += _build_dominance_table(df_sorted)
        for js in land_calls:
            chart_render_calls.append(('tab-landscape', js))

    # Build tab panels HTML
    panels_html = f'<div id="tab-functional" class="tab-panel">{func_html}</div>\n'
    if cell_html:
        panels_html += f'<div id="tab-cellular" class="tab-panel">{cell_html}</div>\n'
    if land_html:
        panels_html += f'<div id="tab-landscape" class="tab-panel">{land_html}</div>\n'

    # Methods (shared, after tabs)
    methods_html = f"""
    <div class="section">
        <h2>Methods</h2>
        <div class="methods-box">
            <p><strong>Categorisation:</strong> Genes were assigned to functional categories
            using a dual-method approach: keyword matching against protein annotations
            (word-boundary regular expressions) and Gene Ontology (GO) ID overlap with
            curated anchor terms expanded through the GO hierarchy.</p>

            <p><strong>Discrete enrichment:</strong> For each expression tier, category
            enrichment was calculated as the fold change of the observed proportion
            relative to the full-dataset background. Statistical significance was assessed
            using Fisher's exact test with Benjamini-Hochberg FDR correction.</p>

            {"<p><strong>Continuous enrichment:</strong> The enrichment function E<sub>C</sub>(k) was computed at regular intervals from k=10 to k=N/2, producing a dense enrichment curve for each category. A permutation-based global profile test assessed whether each category's enrichment profile differs significantly from the null expectation (E=1 for all k), using supremum and integral test statistics.</p>" if has_continuous else ""}

            {"<p><strong>D<sub>KL</sub> gradient:</strong> The Kullback-Leibler divergence from the uniform expectation was computed at each tier size, quantifying overall functional specialisation. A log-linear fit captures the power-law decay of D<sub>KL</sub> with increasing tier size.</p>" if has_continuous else ""}

            <p><strong>Expression tiers:</strong> Genes ranked by TPM (Transcripts Per Million)
            from Salmon quantification. Each tier represents the top N most highly expressed genes.</p>
        </div>
    </div>"""

    # Build JS: group chart renders by tab for lazy loading
    tab_js_map = {}
    for tab_id, js_code in chart_render_calls:
        if tab_id not in tab_js_map:
            tab_js_map[tab_id] = []
        tab_js_map[tab_id].append(js_code)

    # Build the pendingRender map
    pending_lines = []
    for tab_id, js_list in tab_js_map.items():
        fns = ', '.join(
            [f'function() {{ {js} }}' for js in js_list])
        pending_lines.append(f"pendingRender['{tab_id}'] = [{fns}];")
    tab_chart_map_js = '\n    '.join(pending_lines)

    render_js = _JS_TEMPLATE.replace(
        'TAB_CHART_MAP', tab_chart_map_js).replace(
        'CHART_RENDER_CALLS', '// Charts registered via pendingRender above')

    try:
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SCEPTR {report_title}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>{_CSS}</style>
</head>
<body>

<div class="header">
    <h1>SCEPTR {report_title}</h1>
    <div class="subtitle">Continuous enrichment profiling</div>
    <div class="meta">Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}
    &middot; SCEPTR v1.0</div>
</div>

<div class="hero-summary">{hero_text}</div>

<div class="section">
    <div class="stats-grid">
        <div class="stat-card">
            <div class="label">Total Genes</div>
            <div class="value">{total_genes:,}</div>
        </div>
        <div class="stat-card">
            <div class="label">Categorised</div>
            <div class="value">{p_all.get('total_categorised', 0):,}</div>
            <div class="detail">{pct_cat}% of total{_categorisation_context(pct_cat, len(primary['cat_names']))}</div>
        </div>
        <div class="stat-card accent">
            <div class="label">Significant Enrichments</div>
            <div class="value">{sig_count}</div>
            <div class="detail">across {n_tiers} tiers{f" x {len(analysis_blocks)} analyses" if len(analysis_blocks) > 1 else ""}</div>
        </div>
        <div class="stat-card accent">
            <div class="label">{"Significant Profiles" if has_continuous else "Analysis Types"}</div>
            <div class="value">{profile_sig if has_continuous else len(analysis_blocks)}</div>
            <div class="detail">{"global profile test" if has_continuous else "functional + cellular"}</div>
        </div>
    </div>
</div>

{tab_bar_html}

<div class="tab-content-area">
{panels_html}
</div>

{methods_html}

<div class="footer">
    SCEPTR &middot; Statistical Characterisation of Expression Profiles Across Transcriptomic Resolution
    &middot; {datetime.now().strftime('%Y-%m-%d')}
</div>

{render_js}

</body>
</html>"""

        os.makedirs(os.path.dirname(html_path) or '.', exist_ok=True)
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html)

        logger.info(f"Generated interactive report: {html_path}")
        return html_path

    except Exception as e:
        logger.error(f"Error generating interactive report: {e}")
        import traceback
        traceback.print_exc()
        return None
