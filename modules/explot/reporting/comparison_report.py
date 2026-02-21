#!/usr/bin/env python3
"""
HTML report generation for SCEPTR cross-sample comparison.

Generates a self-contained dashboard-style report with base64-embedded figures,
summary statistics, concordance metrics, differential enrichment tables, and
methods documentation. Follows the same design language as the ExPlot reports.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import os
import base64
import logging
from datetime import datetime
from typing import Dict, Any, Optional, List

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.reporting.comparison_report')


def _img_to_base64(filepath: str) -> Optional[str]:
    """Read an image file and return base64 string, or None if missing."""
    if not os.path.exists(filepath):
        return None
    try:
        with open(filepath, 'rb') as f:
            return base64.b64encode(f.read()).decode('utf-8')
    except Exception as e:
        logger.warning(f"Could not read {filepath}: {e}")
        return None


def _format_pval(p: float) -> str:
    """Format p-value for display."""
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


# ---------------------------------------------------------------------------
# CSS (extended from html_report.py with comparison-specific additions)
# ---------------------------------------------------------------------------

_CSS = """
:root {
    --primary: #1a6b54;
    --primary-light: #2d9b78;
    --primary-dark: #124d3c;
    --accent: #e8913a;
    --bg: #f0f2f5;
    --card: #ffffff;
    --text: #2c3e50;
    --text-light: #5a6c7d;
    --border: #e0e4e8;
    --enriched-a: #2980b9;
    --enriched-b: #c0392b;
    --sig-bg: #fef9e7;
    --concordant: #27ae60;
    --divergent: #e74c3c;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: 'Segoe UI', -apple-system, BlinkMacSystemFont, sans-serif;
    line-height: 1.6; color: var(--text);
    max-width: 1300px; margin: 0 auto; padding: 24px;
    background: var(--bg);
}
.header {
    background: linear-gradient(135deg, #2c3e50, #1a252f);
    color: white; padding: 36px 40px; border-radius: 12px;
    margin-bottom: 28px; position: relative; overflow: hidden;
}
.header::after {
    content: ''; position: absolute; top: -40%; right: -10%;
    width: 300px; height: 300px; border-radius: 50%;
    background: rgba(255,255,255,0.06);
}
.header h1 { font-size: 1.9em; font-weight: 700; margin-bottom: 6px; position: relative; z-index: 1; }
.header .subtitle { opacity: 0.85; font-size: 1.05em; position: relative; z-index: 1; }
.header .meta { opacity: 0.65; font-size: 0.85em; margin-top: 10px; position: relative; z-index: 1; }

.section {
    background: var(--card); padding: 28px 32px;
    margin-bottom: 24px; border-radius: 10px;
    box-shadow: 0 1px 6px rgba(0,0,0,0.06);
}
.section h2 {
    color: var(--primary); font-size: 1.35em; font-weight: 700;
    border-bottom: 2px solid var(--primary-light);
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
.stat-card.accent { border-left-color: var(--accent); }
.stat-card.cond-a { border-left-color: var(--enriched-a); }
.stat-card.cond-b { border-left-color: var(--enriched-b); }
.stat-card .label {
    font-size: 0.78em; text-transform: uppercase;
    letter-spacing: 0.5px; color: var(--text-light); font-weight: 600;
}
.stat-card .value {
    font-size: 1.8em; font-weight: 700; color: var(--text);
    margin: 6px 0 2px 0;
}
.stat-card .detail { font-size: 0.82em; color: var(--text-light); }

table {
    width: 100%; border-collapse: collapse;
    margin: 16px 0; font-size: 0.88em;
}
th, td { padding: 10px 14px; text-align: left; border-bottom: 1px solid var(--border); }
th {
    background: var(--primary); color: white;
    font-weight: 600; font-size: 0.85em;
    text-transform: uppercase; letter-spacing: 0.3px;
}
tr:hover { background: #f8fafb; }
tr.sig { background: var(--sig-bg); }
.enriched-a { color: var(--enriched-a); font-weight: 700; }
.enriched-b { color: var(--enriched-b); font-weight: 700; }

.fig-container { text-align: center; margin: 24px 0; }
.fig-container img {
    max-width: 100%; height: auto; border-radius: 8px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08);
}
.fig-caption {
    color: var(--text-light); font-size: 0.88em;
    font-style: italic; margin-top: 8px;
}
.two-col {
    display: grid; grid-template-columns: 1fr 1fr; gap: 24px;
}
@media (max-width: 850px) { .two-col { grid-template-columns: 1fr; } }

.methods-box {
    background: #f8fafb; padding: 20px 24px;
    border-radius: 8px; border-left: 4px solid var(--accent);
    font-size: 0.9em; line-height: 1.7; color: var(--text-light);
}
.methods-box p { margin-bottom: 10px; }
.methods-box strong { color: var(--text); }
.footer {
    text-align: center; padding: 20px; color: var(--text-light);
    font-size: 0.82em;
}
.arrow-up { color: var(--enriched-b); }
.arrow-down { color: var(--enriched-a); }
"""


def _build_figure_html(b64: Optional[str], caption: str, fig_num: int) -> str:
    if not b64:
        return ""
    return f"""
    <div class="fig-container">
        <img src="data:image/png;base64,{b64}" alt="{caption}">
        <p class="fig-caption">Figure {fig_num}: {caption}</p>
    </div>
    """


def _build_concordance_table(concordance_data: List[Dict[str, Any]]) -> str:
    """Build HTML table of concordance metrics per tier."""
    if not concordance_data:
        return "<p>No concordance data available.</p>"

    rows = []
    for entry in concordance_data:
        rho = entry.get('Spearman_Rho', 0)
        ci_lo = entry.get('Spearman_CI_Lower', 0)
        ci_hi = entry.get('Spearman_CI_Upper', 0)
        jaccard = entry.get('Jaccard_Similarity', 0)
        n_shared = entry.get('N_Shared_Genes', 0)

        rho_cls = 'concordant' if rho > 0.7 else ('divergent' if rho < 0.3 else '')
        rho_style = f' style="color:var(--{rho_cls});"' if rho_cls else ''

        rows.append(
            f'<tr>'
            f'<td><strong>Tier {entry["Tier"]}</strong></td>'
            f'<td{rho_style}>{rho:.3f}</td>'
            f'<td>{ci_lo:.3f} \u2013 {ci_hi:.3f}</td>'
            f'<td>{jaccard:.3f}</td>'
            f'<td>{n_shared:,}</td>'
            f'</tr>'
        )

    return f"""
    <table>
        <thead><tr>
            <th>Tier</th><th>Spearman \u03C1</th><th>95% CI</th>
            <th>Jaccard (top genes)</th><th>Shared Genes</th>
        </tr></thead>
        <tbody>{''.join(rows)}</tbody>
    </table>
    <p style="font-size:0.82em;color:var(--text-light);">
        Spearman \u03C1 measures rank correlation of fold-change profiles across categories.
        Jaccard index measures overlap of top-k gene sets between conditions.
        CI = Fisher z-transform 95% confidence interval.
    </p>
    """


def _build_differential_table(diff_df, tiers: List[int],
                               label_a: str, label_b: str) -> str:
    """Build HTML table of differential enrichment results."""
    if diff_df is None or len(diff_df) == 0:
        return "<p>No differential enrichment data available.</p>"

    categories = sorted(diff_df['Category'].unique())

    header = f'<th>Category</th>'
    for t in tiers:
        header += (f'<th>Tier {t}<br>'
                   f'<small>{label_a} | {label_b} | \u0394 (FDR)</small></th>')

    rows = []
    for cat in categories:
        cells = f'<td><strong>{cat}</strong></td>'
        any_sig = False

        for t in tiers:
            row = diff_df[(diff_df['Category'] == cat) & (diff_df['Tier'] == t)]
            if len(row) == 0:
                cells += '<td style="color:var(--text-light)">\u2014</td>'
                continue

            r = row.iloc[0]
            fc_a = r['FC_A']
            fc_b = r['FC_B']
            fc_diff = r['FC_Diff']
            fdr = r['Perm_FDR']
            sig = fdr < 0.05

            if sig:
                any_sig = True

            # Direction arrow
            if fc_diff > 0.1:
                arrow = '<span class="arrow-up">\u25B2</span>'
            elif fc_diff < -0.1:
                arrow = '<span class="arrow-down">\u25BC</span>'
            else:
                arrow = '\u2022'

            sig_marker = ' *' if sig else ''
            fdr_str = _format_pval(fdr)

            cells += (f'<td>'
                      f'<span class="enriched-a">{fc_a:.2f}x</span> | '
                      f'<span class="enriched-b">{fc_b:.2f}x</span><br>'
                      f'<small>{arrow} {fc_diff:+.2f} ({fdr_str}){sig_marker}</small>'
                      f'</td>')

        row_cls = ' class="sig"' if any_sig else ''
        rows.append(f'<tr{row_cls}>{cells}</tr>')

    return f"""
    <table>
        <thead><tr>{header}</tr></thead>
        <tbody>{''.join(rows)}</tbody>
    </table>
    <p style="font-size:0.82em;color:var(--text-light);">
        * Significant (BH-adjusted permutation p &lt; 0.05).
        \u0394 = fold-change difference ({label_b} &minus; {label_a}).
        \u25B2 = enriched in {label_b}, \u25BC = enriched in {label_a}.
    </p>
    """


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
) -> Optional[str]:
    """
    Generate self-contained HTML comparison report.

    Args:
        diff_df: DataFrame with differential enrichment results
        concordance_data: List of per-tier concordance metric dicts
        label_a, label_b: Condition labels
        n_genes_a, n_genes_b: Gene counts per condition
        n_shared: Shared genes between conditions
        n_permutations: Number of permutations used
        tiers: List of tier sizes
        output_prefix: Output path prefix
        figures_dir: Directory containing figure files

    Returns:
        Path to generated HTML file, or None on error.
    """
    html_path = f"{output_prefix}_comparison_report.html"

    # Count significant results
    n_sig = int((diff_df['Perm_FDR'] < 0.05).sum()) if len(diff_df) > 0 else 0
    n_categories = len(diff_df['Category'].unique()) if len(diff_df) > 0 else 0

    # Find figures
    fig_pfx = os.path.join(figures_dir, os.path.basename(output_prefix)) if figures_dir else output_prefix

    images = {}
    for name in ['radar_overlay', 'differential_heatmap', 'fc_barplot',
                  'volcano', 'gradient_overlay']:
        b64 = _img_to_base64(f"{fig_pfx}_{name}.png")
        if b64:
            images[name] = b64

    # Build figure sections
    fig_num = 1

    radar_html = ""
    if 'radar_overlay' in images:
        radar_html = _build_figure_html(
            images['radar_overlay'],
            f"Side-by-side functional profiles: {label_a} vs {label_b}", fig_num)
        fig_num += 1

    heatmap_html = ""
    if 'differential_heatmap' in images:
        heatmap_html = _build_figure_html(
            images['differential_heatmap'],
            "Differential enrichment heatmap (blue = enriched in "
            f"{label_a}, red = enriched in {label_b})", fig_num)
        fig_num += 1

    barplot_html = ""
    if 'fc_barplot' in images:
        barplot_html = _build_figure_html(
            images['fc_barplot'],
            "Grouped fold-change comparison by category and tier", fig_num)
        fig_num += 1

    volcano_html = ""
    if 'volcano' in images:
        volcano_html = _build_figure_html(
            images['volcano'],
            "Volcano plot of fold-change differences vs permutation significance",
            fig_num)
        fig_num += 1

    gradient_html = ""
    if 'gradient_overlay' in images:
        gradient_html = _build_figure_html(
            images['gradient_overlay'],
            "Expression gradient comparison for significantly different categories",
            fig_num)
        fig_num += 1

    concordance_table = _build_concordance_table(concordance_data)
    differential_table = _build_differential_table(diff_df, tiers, label_a, label_b)

    try:
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SCEPTR Comparison: {label_a} vs {label_b}</title>
    <style>{_CSS}</style>
</head>
<body>

<div class="header">
    <h1>SCEPTR Cross-Sample Comparison</h1>
    <div class="subtitle">{label_a} vs {label_b}</div>
    <div class="meta">Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} &middot;
    {n_permutations:,} permutations &middot; SCEPTR Compare</div>
</div>

<div class="section">
    <h2>Comparison Summary</h2>
    <div class="stats-grid">
        <div class="stat-card cond-a">
            <div class="label">{label_a}</div>
            <div class="value">{n_genes_a:,}</div>
            <div class="detail">genes</div>
        </div>
        <div class="stat-card cond-b">
            <div class="label">{label_b}</div>
            <div class="value">{n_genes_b:,}</div>
            <div class="detail">genes</div>
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

<div class="section">
    <h2>Concordance Overview</h2>
    <p>How similar are the two conditions' functional enrichment profiles?
    High Spearman \u03C1 indicates the overall pattern of functional investment
    is conserved; low \u03C1 indicates divergent functional programmes.</p>
    {concordance_table}
</div>

<div class="section">
    <h2>Functional Profile Comparison</h2>
    <p>Side-by-side radar plots showing how each condition allocates transcriptional
    resources across functional categories at each expression tier.</p>
    {radar_html}
</div>

<div class="section">
    <h2>Differential Enrichment Analysis</h2>
    <p>Gene-label permutation test ({n_permutations:,} permutations) with per-tier
    Benjamini-Hochberg FDR correction. Tests whether fold-change differences exceed
    those expected by chance when gene labels are randomly reassigned between conditions.</p>
    {differential_table}
    {heatmap_html}
</div>

<div class="section">
    <h2>Expression Gradient Comparison</h2>
    <p>How enrichment changes across the expression gradient for categories showing
    significant differences between conditions.</p>
    {gradient_html}
</div>

{'<div class="section"><h2>Volcano Plot</h2><p>Overview of all category \u00D7 tier combinations, highlighting significant differences.</p>' + volcano_html + '</div>' if volcano_html else ''}

{'<div class="section"><h2>Grouped Fold-Change Comparison</h2>' + barplot_html + '</div>' if barplot_html else ''}

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

        <p><strong>Categorisation:</strong> Genes are assigned to functional categories
        using SCEPTR's dual-method approach (keyword matching + GO ID hierarchy overlap).
        Category assignments are computed on the shared gene set to ensure comparability.</p>

        <p><strong>Limitations:</strong> This comparison evaluates per-tier differences
        in category enrichment. It does not formally test for differences in enrichment
        gradient shape across tiers. Results from single-sample comparisons without
        replicates should be interpreted as exploratory and validated independently
        where possible.</p>
    </div>
</div>

<div class="section">
    <h2>Interpretation Guide</h2>
    <div class="methods-box">
        <p><strong>Significant positive \u0394:</strong> Category is more enriched in
        {label_b} than {label_a} at that expression tier &mdash; {label_b} invests
        more transcriptional resources in this function among its most highly
        expressed genes.</p>

        <p><strong>Significant negative \u0394:</strong> Category is more enriched in
        {label_a} than {label_b}.</p>

        <p><strong>High concordance (Spearman \u03C1 &gt; 0.7):</strong> The two
        conditions share a broadly similar functional investment pattern, differing
        in specific categories rather than globally.</p>

        <p><strong>Low concordance (Spearman \u03C1 &lt; 0.3):</strong> The two
        conditions have fundamentally different functional investment profiles,
        suggesting a major biological shift.</p>
    </div>
</div>

<div class="footer">
    SCEPTR v1.0.0 &middot; Cross-Sample Comparison Module &middot;
    {datetime.now().strftime('%Y-%m-%d')}
</div>

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
