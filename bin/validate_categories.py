#!/usr/bin/env python3
"""
Validate SCEPTR ExPlot category definitions.

Performs structural checks, keyword analysis, GO ID validation, and
optional coverage analysis against annotation data. Designed to be run
standalone or as part of the pipeline.

Usage:
    python validate_categories.py category.json \
        [--annotation_tsv file.tsv] \
        [--go_obo go-basic.obo] \
        [--output_json report.json] \
        [--category_type functional|cellular]

Author: James McCabe
Module: SCEPTR ExPlot
"""

import os
import sys
import re
import json
import argparse
import logging
from collections import defaultdict
from typing import Dict, List, Set, Any, Optional

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.validate_categories')

# Add explot module to path
_script_dir = os.path.dirname(os.path.abspath(__file__))
_explot_dir = os.path.join(os.path.dirname(_script_dir), 'modules', 'explot')
if os.path.isdir(_explot_dir):
    sys.path.insert(0, _explot_dir)

GO_ID_PATTERN = re.compile(r'GO:\d{7}')


# ---------------------------------------------------------------------------
# Structural validation
# ---------------------------------------------------------------------------

def validate_structure(categories: dict) -> List[str]:
    """Check JSON structure. Returns list of issues."""
    issues = []

    if not isinstance(categories, dict):
        issues.append("Root element is not a dict")
        return issues

    if not categories:
        issues.append("Empty category file")
        return issues

    for cat_name, cat_data in categories.items():
        if isinstance(cat_data, list):
            issues.append(f"'{cat_name}': legacy format (list). Expected v2 dict.")
            continue

        if not isinstance(cat_data, dict):
            issues.append(f"'{cat_name}': unexpected type {type(cat_data).__name__}")
            continue

        if 'keywords' not in cat_data:
            issues.append(f"'{cat_name}': missing 'keywords' key")

        if 'anchor_go_ids' not in cat_data:
            issues.append(f"'{cat_name}': missing 'anchor_go_ids' key")

        kws = cat_data.get('keywords', [])
        if not isinstance(kws, list):
            issues.append(f"'{cat_name}': keywords is not a list")
        elif len(kws) == 0:
            issues.append(f"'{cat_name}': no keywords defined")

        anchors = cat_data.get('anchor_go_ids', [])
        if not isinstance(anchors, list):
            issues.append(f"'{cat_name}': anchor_go_ids is not a list")

        for aid in anchors:
            if not GO_ID_PATTERN.fullmatch(aid):
                issues.append(f"'{cat_name}': malformed GO ID '{aid}'")

    return issues


# ---------------------------------------------------------------------------
# Keyword analysis
# ---------------------------------------------------------------------------

def analyse_keywords(categories: dict) -> Dict[str, Any]:
    """Analyse keyword quality and cross-category overlaps."""
    all_keywords: Dict[str, List[str]] = defaultdict(list)
    short_keywords = []
    single_word_keywords = []
    per_category = {}

    for cat_name, cat_data in categories.items():
        kws = cat_data.get('keywords', []) if isinstance(cat_data, dict) else cat_data
        per_category[cat_name] = len(kws)

        for kw in kws:
            kw_lower = kw.lower().strip()
            all_keywords[kw_lower].append(cat_name)

            if len(kw_lower) < 4:
                short_keywords.append((kw, cat_name))

            if ' ' not in kw_lower and len(kw_lower) > 0:
                single_word_keywords.append((kw, cat_name))

    # Cross-category duplicates
    duplicates = {kw: cats for kw, cats in all_keywords.items() if len(cats) > 1}

    return {
        'total_keywords': sum(per_category.values()),
        'unique_keywords': len(all_keywords),
        'per_category_counts': per_category,
        'cross_category_duplicates': duplicates,
        'short_keywords': short_keywords,
        'single_word_keywords': single_word_keywords,
    }


# ---------------------------------------------------------------------------
# GO ID validation
# ---------------------------------------------------------------------------

def validate_go_ids(categories: dict, go_dag) -> Dict[str, Any]:
    """Validate anchor GO IDs against the GO hierarchy."""
    results = {}
    total_anchors = 0
    total_valid = 0
    total_descendants = 0
    missing_ids = []
    namespace_mismatches = []

    for cat_name, cat_data in categories.items():
        if not isinstance(cat_data, dict):
            continue
        anchors = cat_data.get('anchor_go_ids', [])
        cat_info = {'anchors': [], 'total_descendants': 0}

        for go_id in anchors:
            total_anchors += 1
            if go_id not in go_dag:
                missing_ids.append((go_id, cat_name))
                cat_info['anchors'].append({
                    'id': go_id, 'status': 'MISSING', 'name': None,
                    'namespace': None, 'descendants': 0
                })
                continue

            total_valid += 1
            term = go_dag[go_id]
            desc_count = len(term.get_all_children()) if hasattr(term, 'get_all_children') else 0
            total_descendants += desc_count

            cat_info['anchors'].append({
                'id': go_id,
                'status': 'OK',
                'name': term.name,
                'namespace': term.namespace,
                'descendants': desc_count,
            })
            cat_info['total_descendants'] += desc_count

        results[cat_name] = cat_info

    return {
        'total_anchors': total_anchors,
        'valid_anchors': total_valid,
        'missing_ids': missing_ids,
        'total_descendants': total_descendants,
        'per_category': results,
    }


# ---------------------------------------------------------------------------
# Coverage analysis (requires annotation data)
# ---------------------------------------------------------------------------

def analyse_coverage(categories: dict, annotation_tsv: str, go_dag=None) -> Dict[str, Any]:
    """Run categorisation on annotation data and report coverage."""
    try:
        import categorisation
        import go_utils
    except ImportError:
        logger.error("Cannot import categorisation module. "
                     "Ensure modules/explot is on PYTHONPATH.")
        return {}

    import pandas as pd

    try:
        df = pd.read_csv(annotation_tsv, sep='\t')
    except Exception:
        try:
            df = pd.read_csv(annotation_tsv, sep=',')
        except Exception as e:
            logger.error(f"Cannot read annotation file: {e}")
            return {}

    logger.info(f"Loaded {len(df)} genes from {annotation_tsv}")

    # Prepare GO expansion if anchors exist
    keyword_map, anchor_map = categorisation.normalize_categories(categories)
    expanded_go_sets = None

    has_anchors = any(len(ids) > 0 for ids in anchor_map.values())
    if has_anchors and go_dag is not None:
        expanded_go_sets = go_utils.expand_anchor_go_ids(anchor_map, go_dag)

    # Run categorisation
    results = categorisation.categorise_genes(
        df, categories, allow_multiple=True,
        expanded_go_sets=expanded_go_sets, verbose=False)

    # Build coverage report
    total = len(df)
    categorised = results['total_categorised']
    uncategorised = results['uncategorised']

    method_summary = results.get('method_summary', {})

    return {
        'total_genes': total,
        'categorised': categorised,
        'uncategorised': uncategorised,
        'coverage_pct': round(categorised / total * 100, 1) if total > 0 else 0,
        'multi_category_genes': results['multi_category_genes'],
        'per_category_counts': results['category_counts'],
        'method_summary': method_summary,
    }


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def print_report(structure_issues, keyword_analysis, go_validation, coverage):
    """Print formatted validation report to stdout."""
    print("\n" + "=" * 70)
    print("SCEPTR CATEGORY VALIDATION REPORT")
    print("=" * 70)

    # Structure
    print("\n--- STRUCTURAL VALIDATION ---")
    if structure_issues:
        for issue in structure_issues:
            print(f"  WARNING: {issue}")
    else:
        print("  All checks passed.")

    # Keywords
    print("\n--- KEYWORD ANALYSIS ---")
    ka = keyword_analysis
    print(f"  Total keywords: {ka['total_keywords']}")
    print(f"  Unique keywords: {ka['unique_keywords']}")
    print(f"  Categories: {len(ka['per_category_counts'])}")
    print()
    for cat, count in sorted(ka['per_category_counts'].items()):
        print(f"    {cat}: {count} keywords")

    if ka['cross_category_duplicates']:
        print(f"\n  Cross-category duplicates ({len(ka['cross_category_duplicates'])}):")
        for kw, cats in sorted(ka['cross_category_duplicates'].items()):
            print(f"    \"{kw}\" → {', '.join(cats)}")

    if ka['single_word_keywords']:
        print(f"\n  Single-word keywords ({len(ka['single_word_keywords'])}):")
        for kw, cat in ka['single_word_keywords'][:20]:
            print(f"    \"{kw}\" in {cat}")
        if len(ka['single_word_keywords']) > 20:
            print(f"    ... and {len(ka['single_word_keywords']) - 20} more")

    # GO validation
    if go_validation:
        print("\n--- GO ID VALIDATION ---")
        gv = go_validation
        print(f"  Total anchors: {gv['total_anchors']}")
        print(f"  Valid: {gv['valid_anchors']}")
        print(f"  Missing: {len(gv['missing_ids'])}")
        print(f"  Total descendants: {gv['total_descendants']}")

        if gv['missing_ids']:
            print("\n  Missing GO IDs:")
            for go_id, cat in gv['missing_ids']:
                print(f"    {go_id} in '{cat}'")

        print("\n  Per-category expansion:")
        for cat, info in gv['per_category'].items():
            if info['anchors']:
                desc = info['total_descendants']
                anchors = len(info['anchors'])
                print(f"    {cat}: {anchors} anchors → {desc} descendants")
                for a in info['anchors']:
                    status = a['status']
                    name = a['name'] or '???'
                    ns = a['namespace'] or '???'
                    d = a['descendants']
                    print(f"      {a['id']} ({status}) {ns}: {name} [{d} desc]")

    # Coverage
    if coverage:
        print("\n--- COVERAGE ANALYSIS ---")
        print(f"  Total genes: {coverage['total_genes']}")
        print(f"  Categorised: {coverage['categorised']} ({coverage['coverage_pct']}%)")
        print(f"  Uncategorised: {coverage['uncategorised']}")
        print(f"  Multi-category genes: {coverage['multi_category_genes']}")

        print("\n  Per-category counts:")
        for cat, count in sorted(coverage['per_category_counts'].items(),
                                  key=lambda x: x[1], reverse=True):
            print(f"    {cat}: {count}")

        ms = coverage.get('method_summary', {})
        if ms and any(v.get('go_id_only', 0) + v.get('both', 0) > 0
                      for v in ms.values()):
            print("\n  Assignment method breakdown:")
            print(f"    {'Category':<45} {'KW':>5} {'GO':>5} {'Both':>5}")
            print(f"    {'-'*45} {'-'*5} {'-'*5} {'-'*5}")
            for cat, methods in sorted(ms.items()):
                kw = methods.get('keyword_only', 0)
                go = methods.get('go_id_only', 0)
                both = methods.get('both', 0)
                if kw + go + both > 0:
                    print(f"    {cat:<45} {kw:>5} {go:>5} {both:>5}")

    print("\n" + "=" * 70)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Validate SCEPTR ExPlot category definitions')
    parser.add_argument('category_json', help='Path to category JSON file')
    parser.add_argument('--annotation_tsv', default=None,
                        help='Annotation TSV for coverage analysis')
    parser.add_argument('--go_obo', default=None,
                        help='Path to GO OBO file for ID validation')
    parser.add_argument('--output_json', default=None,
                        help='Save report as JSON')
    args = parser.parse_args()

    # Load categories
    with open(args.category_json) as f:
        categories = json.load(f)
    logger.info(f"Loaded {len(categories)} categories from {args.category_json}")

    # 1. Structural validation
    structure_issues = validate_structure(categories)

    # 2. Keyword analysis
    keyword_analysis = analyse_keywords(categories)

    # 3. GO ID validation (if OBO provided)
    go_validation = None
    go_dag = None
    if args.go_obo:
        try:
            from goatools.obo_parser import GODag
            go_dag = GODag(args.go_obo)
            go_validation = validate_go_ids(categories, go_dag)
        except ImportError:
            logger.warning("goatools not available. Skipping GO validation.")
        except Exception as e:
            logger.error(f"Error loading GO OBO: {e}")

    # 4. Coverage analysis (if annotation TSV provided)
    coverage = None
    if args.annotation_tsv:
        coverage = analyse_coverage(categories, args.annotation_tsv, go_dag)

    # Print report
    print_report(structure_issues, keyword_analysis, go_validation, coverage)

    # Save JSON report
    if args.output_json:
        report = {
            'category_file': args.category_json,
            'structure_issues': structure_issues,
            'keyword_analysis': {
                k: v for k, v in keyword_analysis.items()
                if k not in ('cross_category_duplicates',)
            },
        }
        # Convert non-serializable items
        report['keyword_analysis']['cross_category_duplicates'] = {
            k: v for k, v in keyword_analysis['cross_category_duplicates'].items()
        }
        if go_validation:
            report['go_validation'] = {
                k: v for k, v in go_validation.items()
                if k != 'per_category'
            }
        if coverage:
            report['coverage'] = coverage

        with open(args.output_json, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        logger.info(f"Saved report to {args.output_json}")


if __name__ == '__main__':
    main()
