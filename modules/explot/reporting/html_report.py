#!/usr/bin/env python3
"""
HTML report generation for SCEPTR ExPlot.

Generates self-contained dashboard-style HTML reports with base64-embedded
figures, stat cards, styled enrichment tables, and methods documentation.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import os
import glob
import base64
import logging
from datetime import datetime
from typing import Dict, Any, Optional

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.reporting.html_report')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def _find_images(prefix: str, chart_type: str,
                 figures_dir: Optional[str] = None) -> Dict[str, str]:
    """Find and base64-encode all generated images for a report.

    Args:
        prefix: Output file prefix (e.g. 'sceptr')
        chart_type: 'BP_MF' or 'CC'
        figures_dir: Directory containing figure files. If provided, images
            are looked up as {figures_dir}/{prefix}_{chart_type}_*.png.
            If None, uses {prefix}_{chart_type}_*.png (legacy flat layout).
    """
    images = {}

    # Build figure path prefix
    if figures_dir:
        fig_pfx = os.path.join(figures_dir, os.path.basename(prefix))
    else:
        fig_pfx = prefix

    # Main radar chart
    for pattern in [f"{fig_pfx}_{chart_type}_radar.png",
                    f"{fig_pfx}_{chart_type}_publication_radar.png"]:
        b64 = _img_to_base64(pattern)
        if b64:
            images['radar'] = b64
            break

    # Grouped bar chart
    b64 = _img_to_base64(f"{fig_pfx}_{chart_type}_grouped_bar.png")
    if b64:
        images['grouped_bar'] = b64

    # Multi-tier enrichment chart
    b64 = _img_to_base64(f"{fig_pfx}_{chart_type}_multi_tier_enrichment.png")
    if b64:
        images['multi_tier_enrichment'] = b64

    # Per-tier bar charts
    fig_base = os.path.basename(fig_pfx)
    for png in sorted(glob.glob(f"{fig_pfx}_{chart_type}_top_*_bar.png")):
        tier = os.path.basename(png).replace(f"{fig_base}_{chart_type}_", '').replace('_bar.png', '')
        b64 = _img_to_base64(png)
        if b64:
            images[f'bar_{tier}'] = b64

    # Per-tier enrichment charts
    for png in sorted(glob.glob(f"{fig_pfx}_{chart_type}_top_*_enrichment.png")):
        tier = os.path.basename(png).replace(f"{fig_base}_{chart_type}_", '').replace('_enrichment.png', '')
        b64 = _img_to_base64(png)
        if b64:
            images[f'enrich_{tier}'] = b64

    logger.info(f"Embedded {len(images)} images in report")
    return images


def _format_pval(p: float) -> str:
    """Format p-value for display."""
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


# ---------------------------------------------------------------------------
# CSS
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
    --enriched: #27ae60;
    --depleted: #c0392b;
    --sig-bg: #fef9e7;
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
.enriched { color: var(--enriched); font-weight: 700; }
.depleted { color: var(--depleted); font-weight: 700; }

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
"""


# ---------------------------------------------------------------------------
# Report builders
# ---------------------------------------------------------------------------

def _build_enrichment_table(results: Dict, tier_names: list) -> str:
    """Build HTML enrichment table across all tiers."""
    cats = []
    if tier_names and tier_names[0] in results:
        cats = list(results[tier_names[0]].get('enrichment', {}).keys())

    if not cats:
        return "<p>No enrichment data available.</p>"

    # Check if any tier has core specificity data with non-zero values
    has_core_data = any(
        results.get(tn, {}).get('core_pct', {}).get(cat, 0) > 0
        for tn in tier_names for cat in cats
    )

    header_cells = "<th>Category</th>"
    for tn in tier_names:
        label = tn.replace('top_', 'Top ').title()
        header_cells += f"<th>{label}<br><small>FC (adj. p)</small></th>"

    rows = []
    for cat in cats:
        cells = f"<td><strong>{cat}</strong></td>"
        any_sig = False
        for tn in tier_names:
            stats = results.get(tn, {}).get('enrichment', {}).get(cat, {})
            fc = stats.get('fold_change', 0)
            adj_p = stats.get('adjusted_p_value', stats.get('p_value', 1.0))
            sig = stats.get('significant', False)
            ci_lo = stats.get('ci_lower', 0)
            ci_hi = stats.get('ci_upper', 0)

            if sig:
                any_sig = True

            if fc == 0:
                cells += "<td style='color:var(--text-light)'>—</td>"
            else:
                cls = 'enriched' if fc > 1.5 and sig else ('depleted' if fc < 0.67 and sig else '')
                sig_marker = ' *' if sig else ''
                ci_str = f" ({ci_lo:.1f}–{ci_hi:.1f})" if 0 < ci_hi < 100 else ""
                # Core specificity badge (neutral, no colour judgment)
                core_badge = ''
                if has_core_data:
                    cpct = results.get(tn, {}).get('core_pct', {}).get(cat, 0)
                    if cpct > 0:
                        core_badge = (
                            f'<br><span style="font-size:0.75em;color:var(--text-light);">'
                            f'Core: {cpct:.0f}%</span>')
                cells += (f'<td class="{cls}">'
                          f'{fc:.2f}x{sig_marker}<br>'
                          f'<small>{_format_pval(adj_p)}{ci_str}</small>'
                          f'{core_badge}</td>')

        row_cls = ' class="sig"' if any_sig else ''
        rows.append(f"<tr{row_cls}>{cells}</tr>")

    core_note = ''
    if has_core_data:
        core_note = (' Core % = proportion of genes matched by diagnostic '
                     '(high-specificity) keywords.')

    return f"""
    <table>
        <thead><tr>{header_cells}</tr></thead>
        <tbody>{''.join(rows)}</tbody>
    </table>
    <p style="font-size:0.82em;color:var(--text-light);">
        * Significant (BH-adjusted p &lt; 0.05). FC = fold change.
        CI = 95% exact binomial confidence interval.{core_note}
    </p>
    """


def _build_figure_html(b64: Optional[str], caption: str, fig_num: int) -> str:
    """Build an embedded figure block, or empty if no image."""
    if not b64:
        return ""
    return f"""
    <div class="fig-container">
        <img src="data:image/png;base64,{b64}" alt="{caption}">
        <p class="fig-caption">Figure {fig_num}: {caption}</p>
    </div>
    """


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_functional_report(
    results: Dict[str, Dict[str, Any]],
    all_results: Dict[str, Any],
    output_prefix: str,
    total_genes: int,
    figures_dir: Optional[str] = None,
) -> Optional[str]:
    """Generate dashboard-style HTML report for functional profiling."""
    return _generate_report(results, all_results, output_prefix, total_genes,
                            chart_type="BP_MF", report_title="Functional Profiling",
                            description="biological process and molecular function",
                            figures_dir=figures_dir)


def generate_cellular_report(
    results: Dict[str, Dict[str, Any]],
    all_results: Dict[str, Any],
    output_prefix: str,
    total_genes: int,
    figures_dir: Optional[str] = None,
) -> Optional[str]:
    """Generate dashboard-style HTML report for cellular component profiling."""
    return _generate_report(results, all_results, output_prefix, total_genes,
                            chart_type="CC", report_title="Cellular Component Profiling",
                            description="cellular compartment localisation",
                            figures_dir=figures_dir)


def _build_method_summary_html(method_summary: Dict[str, Dict[str, int]],
                               core_summary: Optional[Dict[str, Dict[str, int]]] = None) -> str:
    """Build HTML table showing keyword vs GO ID assignment breakdown per category."""
    if not method_summary:
        return ""

    # Check if any GO-based assignments exist
    has_go = any(
        v.get('go_id_only', 0) + v.get('both', 0) > 0
        for v in method_summary.values()
    )
    if not has_go:
        return ""

    # Check if core specificity data exists with non-zero core counts
    has_core = (core_summary is not None and
                any(cs.get('core', 0) > 0 for cs in core_summary.values()))

    rows = []
    total_kw = total_go = total_both = total_core = total_ext = 0
    for cat, methods in sorted(method_summary.items()):
        kw = methods.get('keyword_only', 0)
        go = methods.get('go_id_only', 0)
        both = methods.get('both', 0)
        total = kw + go + both
        if total == 0:
            continue
        total_kw += kw
        total_go += go
        total_both += both

        core_cells = ''
        if has_core:
            cs = core_summary.get(cat, {})
            c = cs.get('core', 0)
            e = cs.get('extended', 0)
            total_core += c
            total_ext += e
            core_cells = f"<td>{c}</td><td>{e}</td>"

        rows.append(
            f"<tr><td><strong>{cat}</strong></td>"
            f"<td>{kw}</td><td>{go}</td><td>{both}</td>"
            f"{core_cells}<td>{total}</td></tr>"
        )

    if not rows:
        return ""

    core_total_cells = ''
    if has_core:
        core_total_cells = f"<td>{total_core}</td><td>{total_ext}</td>"

    rows.append(
        f'<tr style="font-weight:700;border-top:2px solid var(--primary);">'
        f"<td>Total</td><td>{total_kw}</td><td>{total_go}</td>"
        f"<td>{total_both}</td>{core_total_cells}"
        f"<td>{total_kw + total_go + total_both}</td></tr>"
    )

    core_headers = ''
    if has_core:
        core_headers = '<th>Core</th><th>Extended</th>'

    core_legend = ''
    if has_core:
        core_legend = (' Core = matched by a high-specificity diagnostic keyword.'
                       ' Extended = matched by a broader keyword or GO ID only.')

    return f"""
    <h3>Assignment Method Breakdown</h3>
    <p>Genes are assigned to categories by keyword matching against annotations
    and/or GO ID overlap with curated anchor terms. This table shows how each
    method contributes.</p>
    <table>
        <thead><tr>
            <th>Category</th><th>Keyword Only</th><th>GO ID Only</th>
            <th>Both</th>{core_headers}<th>Total</th>
        </tr></thead>
        <tbody>{''.join(rows)}</tbody>
    </table>
    <p style="font-size:0.82em;color:var(--text-light);">
        Keyword Only = matched by annotation text only.
        GO ID Only = matched by GO hierarchy overlap only.
        Both = matched by both methods (agreement).{core_legend}
    </p>
    """


def _generate_report(
    results: Dict, all_results: Dict, output_prefix: str,
    total_genes: int, chart_type: str, report_title: str, description: str,
    figures_dir: Optional[str] = None,
) -> Optional[str]:
    """Core report generator."""

    html_path = f"{output_prefix}_{chart_type}_report.html"
    tier_names = sorted(results.keys(), key=lambda x: int(x.split('_')[-1]))
    images = _find_images(output_prefix, chart_type, figures_dir=figures_dir)

    total_cat = all_results.get('total_categorised', 0)
    uncat = all_results.get('uncategorised', 0)
    multi = all_results.get('multi_category_genes', 0)
    pct_cat = round(total_cat / total_genes * 100, 1) if total_genes > 0 else 0
    multi_pct = round(multi / total_cat * 100, 1) if total_cat > 0 else 0

    sig_count = sum(
        1 for tn in tier_names
        for stats in results.get(tn, {}).get('enrichment', {}).values()
        if stats.get('significant', False)
    )

    # Build figure sections
    fig_num = 1

    radar_html = ""
    if 'radar' in images:
        radar_html = _build_figure_html(
            images['radar'],
            f"{report_title.replace(' Profiling', '')} distribution across expression tiers",
            fig_num)
        fig_num += 1

    grouped_html = ""
    if 'grouped_bar' in images:
        grouped_html = _build_figure_html(
            images['grouped_bar'],
            "Grouped comparison of category percentages across tiers",
            fig_num)
        fig_num += 1

    # Per-tier figures
    tier_figs_html = ""
    for tn in tier_names:
        bar_b64 = images.get(f'bar_{tn}')
        enrich_b64 = images.get(f'enrich_{tn}')
        if bar_b64 or enrich_b64:
            label = tn.replace('top_', 'Top ').title()
            tier_figs_html += f'<h3>{label}</h3><div class="two-col">'
            if bar_b64:
                tier_figs_html += _build_figure_html(bar_b64, f"Category distribution — {label}", fig_num)
                fig_num += 1
            if enrich_b64:
                tier_figs_html += _build_figure_html(enrich_b64, f"Fold enrichment — {label}", fig_num)
                fig_num += 1
            tier_figs_html += '</div>'

    enrichment_table = _build_enrichment_table(results, tier_names)
    method_summary_html = _build_method_summary_html(
        all_results.get('method_summary', {}),
        core_summary=all_results.get('core_summary'))

    try:
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SCEPTR {report_title} Report</title>
    <style>{_CSS}</style>
</head>
<body>

<div class="header">
    <h1>SCEPTR {report_title} Report</h1>
    <div class="subtitle">Multi-resolution {description} analysis</div>
    <div class="meta">Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} &middot; SCEPTR</div>
</div>

<div class="section">
    <h2>Dataset Summary</h2>
    <div class="stats-grid">
        <div class="stat-card">
            <div class="label">Total Genes</div>
            <div class="value">{total_genes:,}</div>
        </div>
        <div class="stat-card">
            <div class="label">Categorised</div>
            <div class="value">{total_cat:,}</div>
            <div class="detail">{pct_cat}% of total</div>
        </div>
        <div class="stat-card accent">
            <div class="label">Multi-Category</div>
            <div class="value">{multi:,}</div>
            <div class="detail">{multi_pct}% of categorised</div>
        </div>
        <div class="stat-card accent">
            <div class="label">Significant Enrichments</div>
            <div class="value">{sig_count}</div>
            <div class="detail">across all tiers</div>
        </div>
    </div>
</div>

<div class="section">
    <h2>Expression Tier Overview</h2>
    <p>Category distribution across expression tiers, showing how {description}
    patterns shift from the most highly expressed genes to the broader transcriptome.</p>
    {radar_html}
    {grouped_html}
</div>

<div class="section">
    <h2>Enrichment Analysis</h2>
    <p>Fold enrichment of each category in expression tiers relative to the full dataset.
    Statistical significance assessed by Fisher's exact test with Benjamini-Hochberg
    FDR correction.</p>
    {enrichment_table}
</div>

{'<div class="section"><h2>Assignment Methods</h2>' + method_summary_html + '</div>' if method_summary_html else ''}

<div class="section">
    <h2>Tier Detail</h2>
    <p>Per-tier category distributions and enrichment patterns.</p>
    {tier_figs_html if tier_figs_html else '<p style="color:var(--text-light)">No per-tier figures generated.</p>'}
</div>

<div class="section">
    <h2>Methods</h2>
    <div class="methods-box">
        <p><strong>Categorisation:</strong> Genes were assigned to categories using
        a dual-method approach: (1) keyword matching against protein annotations
        using word-boundary regular expressions, and (2) Gene Ontology (GO) ID
        overlap between gene annotations and curated anchor GO terms expanded
        through the GO hierarchy. A gene is assigned to a category if either
        method matches. Multi-category assignment is enabled by default,
        reflecting the multi-functional nature of many proteins.</p>

        <p><strong>Enrichment:</strong> For each expression tier, category enrichment
        was calculated as the fold change of the observed proportion relative to
        the background (full dataset) proportion. Statistical significance was
        assessed using Fisher's exact test (one-sided, greater) with
        Benjamini-Hochberg FDR correction across all categories.</p>

        <p><strong>Confidence intervals:</strong> 95% exact binomial
        (Clopper-Pearson) confidence intervals are reported for fold changes.</p>

        <p><strong>Multi-category note:</strong> When genes belong to multiple
        categories, category counts may sum to more than the total number of genes.
        Fisher's exact test assumes independent draws, which is approximated but
        not strictly met. Enrichment statistics should be interpreted as approximate
        in multi-category mode.</p>

        <p><strong>Expression tiers:</strong> Genes are ranked by TPM (Transcripts
        Per Million) from Salmon quantification. Each tier represents the top N
        most highly expressed genes.</p>

        <p><strong>Core specificity:</strong> Categories may optionally define
        high-confidence &ldquo;core&rdquo; keywords that unambiguously identify
        genes belonging to that function, alongside broader &ldquo;extended&rdquo;
        keywords and GO ancestry terms. The core percentage (Core %) indicates
        what proportion of genes in each enriched category were matched by
        diagnostic core terms. Enrichment statistics are computed identically
        regardless of core/extended status; the distinction is purely a
        specificity diagnostic.</p>
    </div>
</div>

<div class="section">
    <h2>Interpretation Guide</h2>
    <div class="methods-box">
        <p><strong>Fold change &gt; 1:</strong> Category is enriched in the expression
        tier relative to the background &mdash; these functions are disproportionately
        represented among highly expressed genes.</p>

        <p><strong>Fold change &lt; 1:</strong> Category is depleted &mdash; under-represented
        among highly expressed genes.</p>

        <p><strong>Limitations:</strong> This is a single-sample analysis. Results
        characterise expression-dependent functional enrichment patterns but require
        independent validation. Annotation bias may affect results &mdash; well-studied
        pathways may appear more enriched due to more complete annotation.</p>
    </div>
</div>

<div class="footer">
    SCEPTR v1.0.0 &middot; Multi-Resolution Enrichment Profiling
    &middot; {datetime.now().strftime('%Y-%m-%d')}
</div>

</body>
</html>"""

        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html)

        logger.info(f"Generated {report_title} report: {html_path}")
        return html_path

    except Exception as e:
        logger.error(f"Error generating report: {e}")
        import traceback
        traceback.print_exc()
        return None
