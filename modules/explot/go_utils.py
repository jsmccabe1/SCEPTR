#!/usr/bin/env python3
"""
GO Term Utilities for SCEPTR ExPlot.

Provides functionality for:
- Loading and parsing Gene Ontology (GO) hierarchies from OBO files
- Expanding functional categories with GO terms (conservative, word-boundary)
- Reporting expansion statistics for transparency

Author: James McCabe
Module: SCEPTR ExPlot
"""

import os
import re
import sys
import logging
from typing import List, Dict, Optional, Set

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.go_utils')

try:
    from goatools.obo_parser import GODag
except ImportError as e:
    logger.error(f"Missing required package: {e}")
    logger.error("Install goatools: pip install goatools")
    sys.exit(1)


# ---------------------------------------------------------------------------
# GO DAG loading
# ---------------------------------------------------------------------------

def load_go_dag(ontology: Optional[str] = None,
                paths: Optional[List[str]] = None) -> GODag:
    """Load Gene Ontology DAG from a local OBO file."""
    if paths is None:
        paths = [
            os.environ.get('GO_OBO_PATH', '/data/go/go-basic.obo'),
            '/data/go/go-basic.obo',
            os.path.join(os.getcwd(), 'go-basic.obo'),
            '/SCEPTR/data/go/go-basic.obo',
            os.path.join(os.path.dirname(__file__), 'cli', 'go-basic.obo'),
        ]

    for path in paths:
        if os.path.exists(path):
            logger.info(f"Found GO OBO file at: {path}")
            try:
                go_dag = GODag(path)
                logger.info(f"GO hierarchy loaded ({len(go_dag)} terms)")
                if ontology:
                    ns_map = {'BP': 'biological_process',
                              'MF': 'molecular_function',
                              'CC': 'cellular_component'}
                    ns = ns_map.get(ontology, ontology)
                    n = sum(1 for t in go_dag.values()
                            if hasattr(t, 'namespace') and t.namespace == ns)
                    logger.info(f"  {ontology} ({ns}): {n} terms")
                return go_dag
            except Exception as e:
                logger.warning(f"Error loading GO OBO from {path}: {e}")

    raise FileNotFoundError(
        f"GO OBO file not found in: {', '.join(paths)}. "
        f"Set GO_OBO_PATH or place go-basic.obo in the working directory.")


# ---------------------------------------------------------------------------
# Keyword expansion
# ---------------------------------------------------------------------------

def _get_direct_children(go_dag, go_id: str) -> Set[str]:
    """Get immediate children of a GO term (one level down)."""
    children = set()
    term = go_dag.get(go_id)
    if term is None:
        return children
    # goatools stores children via the 'children' attribute
    if hasattr(term, 'children'):
        for child in term.children:
            if hasattr(child, 'id'):
                children.add(child.id)
    return children


def _collect_descendants_depth_limited(go_dag, go_id: str,
                                        max_depth: int) -> Set[str]:
    """Collect descendant GO IDs up to max_depth levels."""
    result = set()
    frontier = {go_id}
    for _ in range(max_depth):
        next_frontier = set()
        for fid in frontier:
            for child_id in _get_direct_children(go_dag, fid):
                if child_id not in result:
                    result.add(child_id)
                    next_frontier.add(child_id)
        frontier = next_frontier
        if not frontier:
            break
    return result


def expand_keywords_with_go_terms(
    keywords: List[str],
    go_dag: GODag,
    include_descendants: bool = True,
    max_descendant_depth: int = 2,
) -> List[str]:
    """
    Expand keywords with related GO term names.

    Matching rules (conservative):
    - Only matches on GO term *names*, NOT definitions
    - Uses word-boundary regex to prevent substring false positives
    - Descendant traversal is depth-limited (default: 2 levels)

    Args:
        keywords: Seed keyword strings
        go_dag: GO hierarchy object
        include_descendants: Whether to add descendant term names
        max_descendant_depth: Maximum descendant levels (default: 2)

    Returns:
        Expanded keyword list including matched GO term names
    """
    if not go_dag or len(go_dag) == 0:
        logger.warning("Empty GO DAG, no expansion possible")
        return keywords

    expanded = set(keywords)
    matched_go_ids = set()

    for keyword in keywords:
        kw_lower = keyword.lower().strip()
        if not kw_lower:
            continue

        pattern = re.compile(r'\b' + re.escape(kw_lower) + r'\b', re.IGNORECASE)

        for go_id, go_term in go_dag.items():
            if not hasattr(go_term, 'name') or not go_term.name:
                continue
            if pattern.search(go_term.name):
                expanded.add(go_term.name)
                matched_go_ids.add(go_id)

    if include_descendants and matched_go_ids:
        prev_size = len(expanded)
        for seed_id in list(matched_go_ids):
            for desc_id in _collect_descendants_depth_limited(
                    go_dag, seed_id, max_descendant_depth):
                if desc_id in go_dag and hasattr(go_dag[desc_id], 'name'):
                    expanded.add(go_dag[desc_id].name)
        added = len(expanded) - prev_size
        if added:
            logger.info(f"  +{added} terms from descendants (depth ≤{max_descendant_depth})")

    logger.info(f"  Expansion: {len(keywords)} → {len(expanded)} keywords "
                f"({len(matched_go_ids)} GO terms matched)")
    return list(expanded)


def expand_category_keywords(
    categories: Dict[str, List[str]],
    go_dag: GODag,
    max_descendant_depth: int = 2,
) -> Dict[str, List[str]]:
    """Expand all categories with GO terms."""
    expanded = {}
    for category, keywords in categories.items():
        logger.info(f"Expanding category: {category}")
        expanded[category] = expand_keywords_with_go_terms(
            keywords, go_dag,
            include_descendants=True,
            max_descendant_depth=max_descendant_depth)
    return expanded


# ---------------------------------------------------------------------------
# Simple helpers
# ---------------------------------------------------------------------------

def get_go_term_ancestors(go_id: str, go_dag: GODag) -> Set[str]:
    """Get all ancestor term IDs."""
    return set(go_dag[go_id].get_all_parents()) if go_id in go_dag else set()


def get_go_term_descendants(go_id: str, go_dag: GODag) -> Set[str]:
    """Get all descendant term IDs."""
    return set(go_dag[go_id].get_all_children()) if go_id in go_dag else set()


# ---------------------------------------------------------------------------
# Anchor GO ID expansion (for v2 category format)
# ---------------------------------------------------------------------------

def expand_anchor_go_ids(
    anchor_ids_by_category: Dict[str, List[str]],
    go_dag: GODag,
) -> Dict[str, Set[str]]:
    """
    Expand anchor GO IDs to include all descendants for each category.

    Unlike keyword-based expansion (which matches keyword text against GO
    term names), this performs direct GO hierarchy traversal from hand-curated
    anchor IDs. Uses unlimited depth since anchors are chosen at the
    appropriate level of specificity.

    Args:
        anchor_ids_by_category: {category_name: [GO:XXXXXXX, ...]}
        go_dag: GO hierarchy object

    Returns:
        {category_name: set of GO IDs (anchors + all descendants)}
    """
    expanded = {}
    total_anchors = 0
    total_expanded = 0

    for category, anchor_ids in anchor_ids_by_category.items():
        if not anchor_ids:
            expanded[category] = set()
            continue

        cat_ids = set()
        valid_anchors = 0
        for go_id in anchor_ids:
            if go_id not in go_dag:
                logger.warning(f"  {category}: anchor {go_id} not found in GO DAG - skipping")
                continue
            valid_anchors += 1
            cat_ids.add(go_id)
            descendants = get_go_term_descendants(go_id, go_dag)
            cat_ids.update(descendants)

        expanded[category] = cat_ids
        total_anchors += valid_anchors
        total_expanded += len(cat_ids)

        if valid_anchors > 0:
            logger.info(f"  {category}: {valid_anchors} anchors → {len(cat_ids)} GO IDs")

    logger.info(f"Anchor GO expansion: {total_anchors} anchors → {total_expanded} total GO IDs "
                f"across {sum(1 for s in expanded.values() if s)} categories")
    return expanded


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='SCEPTR GO Term Utilities')
    parser.add_argument('--ontology', choices=['BP', 'MF', 'CC'])
    parser.add_argument('--keywords', nargs='+', help='Keywords to expand')
    parser.add_argument('--max_depth', type=int, default=2)
    args = parser.parse_args()

    dag = load_go_dag(ontology=args.ontology)
    if args.keywords:
        result = expand_keywords_with_go_terms(
            args.keywords, dag, max_descendant_depth=args.max_depth)
        print(f"Original: {args.keywords}")
        print(f"Expanded: {len(result)} keywords")
