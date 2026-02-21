#!/usr/bin/env python3
"""
Gene categorisation for SCEPTR ExPlot.

Assigns genes to functional or cellular categories using a dual-method
approach: word-boundary keyword matching against annotation text AND/OR
GO ID overlap with anchor GO terms expanded via the GO hierarchy.

Supports both legacy format ({"Cat": ["kw1", ...]}) and v2 format
({"Cat": {"keywords": [...], "anchor_go_ids": [...]}}).

Author: James McCabe
Module: SCEPTR ExPlot
"""

import re
import logging
from typing import Dict, List, Set, Tuple, Any, Optional
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('sceptr.explot.categorisation')


# ---------------------------------------------------------------------------
# GO ID columns used for structural GO-based assignment
# ---------------------------------------------------------------------------

_DEFAULT_GO_COLUMNS = [
    'GO_Biological_Process', 'GO_Molecular_Function', 'GO_Cellular_Component',
    'GO_terms', 'GO', 'go_terms', 'GO_annotation',
]

_GO_ID_PATTERN = re.compile(r'GO:\d{7}')


# ---------------------------------------------------------------------------
# Category format detection and normalisation
# ---------------------------------------------------------------------------

def detect_category_format(categories: Dict[str, Any]) -> str:
    """Detect whether categories use legacy (list) or v2 (dict) format."""
    if not categories:
        return 'legacy'
    first_val = next(iter(categories.values()))
    if isinstance(first_val, dict) and 'keywords' in first_val:
        return 'v2'
    return 'legacy'


def normalize_categories(
    categories: Dict[str, Any],
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, frozenset]]:
    """
    Split categories into keyword map, anchor GO ID map, and core keyword set.

    Returns:
        (keyword_map, anchor_id_map, core_keyword_sets) where keyword_map and
        anchor_id_map are {category: list} and core_keyword_sets is
        {category: frozenset of lowercased core keywords}.
    """
    fmt = detect_category_format(categories)
    if fmt == 'v2':
        keyword_map = {cat: data['keywords'] for cat, data in categories.items()}
        anchor_map = {cat: data.get('anchor_go_ids', []) for cat, data in categories.items()}
        core_sets = {
            cat: frozenset(k.lower().strip() for k in data.get('core_keywords', []) if k.strip())
            for cat, data in categories.items()
        }
    else:
        keyword_map = categories
        anchor_map = {cat: [] for cat in categories}
        core_sets = {cat: frozenset() for cat in categories}
    return keyword_map, anchor_map, core_sets


# ---------------------------------------------------------------------------
# GO ID extraction from annotation rows
# ---------------------------------------------------------------------------

def extract_go_ids_from_row(
    row: pd.Series,
    go_columns: Optional[List[str]] = None,
) -> Set[str]:
    """Extract GO IDs from a gene's annotation columns.

    Parses entries like 'glucose metabolic process [GO:0006006]; ...'
    and returns a set of GO ID strings.
    """
    columns = go_columns or _DEFAULT_GO_COLUMNS
    go_ids = set()
    for col in columns:
        val = row.get(col)
        if pd.notna(val):
            go_ids.update(_GO_ID_PATTERN.findall(str(val)))
    return go_ids


# ---------------------------------------------------------------------------
# Annotation extraction
# ---------------------------------------------------------------------------

# Standard annotation fields to check (order matters: more specific first)
_ANNOTATION_FIELDS = [
    'Protein_Name', 'protein_name', 'protein_description',
    'annotation', 'description', 'product',
    'GO_Biological_Process', 'GO_Molecular_Function', 'GO_Cellular_Component',
    'GO_terms', 'GO', 'go_terms', 'GO_annotation',
    'gene_name', 'Name', 'name', 'function', 'Function',
]


def identify_available_annotation_fields(df: pd.DataFrame) -> List[str]:
    """Identify which annotation fields are present in the DataFrame."""
    available = [f for f in _ANNOTATION_FIELDS if f in df.columns]
    if not available:
        logger.warning("No standard annotation fields found!")
        logger.info(f"Available columns: {df.columns.tolist()}")
    else:
        logger.info(f"Annotation fields: {available}")
    return available


def get_annotation_text(row: pd.Series, fields: List[str]) -> str:
    """Extract and combine annotation text from a gene row."""
    parts = []
    for field in fields:
        val = row.get(field)
        if pd.notna(val):
            parts.append(str(val))
    return ' '.join(parts).lower()


# ---------------------------------------------------------------------------
# Keyword matching (word-boundary)
# ---------------------------------------------------------------------------

def match_keywords_to_annotation(
    annotation_text: str,
    keywords: List[str],
    require_word_boundary: bool = True,
) -> Optional[str]:
    """
    Check if any keyword matches the annotation text.

    Args:
        annotation_text: Lowercased text to search in
        keywords: Keywords to search for
        require_word_boundary: Use \\b regex (default True)

    Returns:
        First matching keyword, or None
    """
    for keyword in keywords:
        kw = keyword.lower().strip()
        if not kw:
            continue

        if require_word_boundary:
            if re.search(r'\b' + re.escape(kw) + r'\b', annotation_text):
                return keyword
        else:
            if kw in annotation_text:
                return keyword

    return None


# ---------------------------------------------------------------------------
# Core categorisation (dual-method: keyword + GO ID)
# ---------------------------------------------------------------------------

def categorise_genes(
    df: pd.DataFrame,
    categories: Dict[str, Any],
    allow_multiple: bool = True,
    require_word_boundary: bool = True,
    max_example_matches: int = 3,
    verbose: bool = True,
    expanded_go_sets: Optional[Dict[str, Set[str]]] = None,
    go_columns: Optional[List[str]] = None,
    core_keyword_sets: Optional[Dict[str, frozenset]] = None,
) -> Dict[str, Any]:
    """
    Categorise genes using keyword matching and/or GO ID overlap.

    Args:
        df: DataFrame with gene annotations
        categories: Category definitions (legacy or v2 format)
        allow_multiple: Allow genes in multiple categories
        require_word_boundary: Word-boundary matching (default True)
        max_example_matches: Example matches to log per category
        verbose: Detailed logging
        expanded_go_sets: {category: set of expanded GO IDs} for GO-based
            assignment. If None, only keyword matching is used.
        go_columns: Columns to extract gene GO IDs from. Defaults to
            standard GO annotation columns.
        core_keyword_sets: {category: frozenset of lowercased core keywords}.
            If provided, tracks whether each gene-category match is via a
            core (diagnostic) keyword or an extended (broad) keyword/GO.

    Returns:
        dict with keys: category_counts, total_categorised, uncategorised,
                        genes_in_category, multi_category_genes, category_overlaps,
                        gene_categories, method_used, method_summary,
                        assignment_methods, core_match, core_summary
    """
    keyword_map, _, _core_sets = normalize_categories(categories)
    # Use caller-provided core sets if available, else fall back to parsed ones
    if core_keyword_sets is None:
        core_keyword_sets = _core_sets

    category_counts = {cat: 0 for cat in keyword_map}
    genes_in_category = {cat: [] for cat in keyword_map}
    categorised_genes: Set[str] = set()
    gene_categories: Dict[str, List[str]] = {}
    multi_category_genes = 0

    # Method tracking
    use_go = expanded_go_sets is not None and any(
        len(s) > 0 for s in expanded_go_sets.values())
    assignment_methods: Dict[str, Dict[str, str]] = {}
    method_summary: Dict[str, Dict[str, int]] = {
        cat: {'keyword_only': 0, 'go_id_only': 0, 'both': 0}
        for cat in keyword_map
    }

    # Core specificity tracking
    has_core = core_keyword_sets is not None and any(
        len(s) > 0 for s in core_keyword_sets.values())
    core_match: Dict[str, Dict[str, bool]] = {}
    core_summary: Dict[str, Dict[str, int]] = {
        cat: {'core': 0, 'extended': 0} for cat in keyword_map
    }

    available_fields = identify_available_annotation_fields(df)

    if not available_fields and len(df.columns) > 0:
        available_fields = [df.columns[0]]
        logger.warning(f"Using {df.columns[0]} as fallback annotation field")

    # Identify available GO columns for extraction
    available_go_cols = None
    if use_go:
        cols = go_columns or _DEFAULT_GO_COLUMNS
        available_go_cols = [c for c in cols if c in df.columns]
        if not available_go_cols:
            logger.warning("No GO ID columns found in data — GO-based assignment disabled")
            use_go = False

    for idx, row in df.iterrows():
        gene_id = str(row.get('gene_id', idx))
        annotation = get_annotation_text(row, available_fields)

        if not annotation.strip() and len(df.columns) > 0:
            first_col = df.columns[0]
            if pd.notna(row.get(first_col)):
                annotation += ' ' + str(row[first_col]).lower()

        # Extract gene GO IDs once per row
        gene_go_ids = None
        if use_go:
            gene_go_ids = extract_go_ids_from_row(row, available_go_cols)

        matched = []
        gene_methods: Dict[str, str] = {}

        for category, keywords in keyword_map.items():
            kw_hit = match_keywords_to_annotation(
                annotation, keywords, require_word_boundary)

            go_hit = False
            if use_go and gene_go_ids and expanded_go_sets.get(category):
                go_hit = bool(gene_go_ids & expanded_go_sets[category])

            if kw_hit or go_hit:
                if category not in matched:
                    category_counts[category] += 1
                    genes_in_category[category].append(gene_id)
                    matched.append(category)

                    # Determine method
                    if kw_hit and go_hit:
                        method = 'both'
                    elif kw_hit:
                        method = 'keyword_only'
                    else:
                        method = 'go_id_only'
                    gene_methods[category] = method
                    method_summary[category][method] += 1

                    # Core specificity: keyword match against core set
                    is_core = (
                        kw_hit is not None
                        and has_core
                        and kw_hit.lower().strip() in core_keyword_sets.get(category, frozenset())
                    )
                    if gene_id not in core_match:
                        core_match[gene_id] = {}
                    core_match[gene_id][category] = is_core
                    if is_core:
                        core_summary[category]['core'] += 1
                    else:
                        core_summary[category]['extended'] += 1

                    if verbose and len(genes_in_category[category]) <= max_example_matches:
                        via = f"'{kw_hit}'" if kw_hit else 'GO ID'
                        core_tag = ' [core]' if is_core else ''
                        logger.info(f"Match: {gene_id} → '{category}' via {via} [{method}]{core_tag}")

                    if not allow_multiple:
                        break

        if matched:
            categorised_genes.add(gene_id)
            gene_categories[gene_id] = matched
            assignment_methods[gene_id] = gene_methods
            if len(matched) > 1:
                multi_category_genes += 1

    total_categorised = len(categorised_genes)
    uncategorised = len(df) - total_categorised

    # Category overlaps
    category_overlaps: Dict[Tuple[str, str], int] = {}
    if allow_multiple and multi_category_genes > 0:
        for cats in gene_categories.values():
            if len(cats) > 1:
                for i, c1 in enumerate(cats):
                    for c2 in cats[i + 1:]:
                        pair = tuple(sorted([c1, c2]))
                        category_overlaps[pair] = category_overlaps.get(pair, 0) + 1

    if verbose:
        logger.info(f"Categorisation summary:")
        logger.info(f"  Total genes: {len(df)}")
        pct = round(total_categorised / len(df) * 100, 1) if len(df) > 0 else 0
        logger.info(f"  Categorised: {total_categorised} ({pct}%)")
        logger.info(f"  Uncategorised: {uncategorised}")
        if total_categorised > 0:
            mc_pct = round(multi_category_genes / total_categorised * 100, 1)
            logger.info(f"  Multi-category: {multi_category_genes} ({mc_pct}% of categorised)")
        for cat, count in category_counts.items():
            ms = method_summary[cat]
            parts = [f"{cat}: {count}"]
            if use_go and count > 0:
                parts.append(f"(kw={ms['keyword_only']}, go={ms['go_id_only']}, both={ms['both']})")
            logger.info(f"  {' '.join(parts)}")

    return {
        'category_counts': category_counts,
        'total_categorised': total_categorised,
        'uncategorised': uncategorised,
        'genes_in_category': genes_in_category,
        'multi_category_genes': multi_category_genes,
        'category_overlaps': category_overlaps,
        'gene_categories': gene_categories,
        'method_used': 'dual' if use_go else 'standard',
        'method_summary': method_summary,
        'assignment_methods': assignment_methods,
        'core_match': core_match,
        'core_summary': core_summary,
    }


# ---------------------------------------------------------------------------
# Percentage calculation
# ---------------------------------------------------------------------------

def calculate_percentages(
    category_counts: Dict[str, int],
    total_genes: int,
    method: str = 'multi',
) -> Dict[str, float]:
    """
    Calculate percentages for each category.

    Args:
        category_counts: {category: count}
        total_genes: Total number of genes
        method: 'multi' (% of all assignments), 'single' (% of total genes),
                or 'relative' (% of categorised genes)
    """
    if method == 'single':
        if total_genes == 0:
            return {cat: 0.0 for cat in category_counts}
        return {cat: count / total_genes * 100 for cat, count in category_counts.items()}

    elif method == 'multi':
        total_assignments = sum(category_counts.values())
        if total_assignments == 0:
            return {cat: 0.0 for cat in category_counts}
        return {cat: count / total_assignments * 100 for cat, count in category_counts.items()}

    elif method == 'relative':
        categorised = sum(1 for c in category_counts.values() if c > 0)
        if categorised == 0:
            return {cat: 0.0 for cat in category_counts}
        return {cat: count / categorised * 100 for cat, count in category_counts.items()}

    else:
        raise ValueError(f"Unknown percentage method: {method}")


# ---------------------------------------------------------------------------
# Fallback categorisation (uses word boundaries, searches all text columns)
# ---------------------------------------------------------------------------

def fallback_categorisation(
    df: pd.DataFrame,
    categories: Dict[str, Any],
    allow_multiple: bool = True,
    expanded_go_sets: Optional[Dict[str, Set[str]]] = None,
    go_columns: Optional[List[str]] = None,
    core_keyword_sets: Optional[Dict[str, frozenset]] = None,
) -> Dict[str, Any]:
    """
    Fallback: search ALL text columns with word-boundary matching.

    Used when standard categorisation finds <5% of genes categorised.
    Unlike the standard approach, this checks every text column rather
    than just known annotation fields.
    """
    logger.warning("Using fallback categorisation (all text columns, word-boundary)...")

    keyword_map, _, _core_sets = normalize_categories(categories)
    if core_keyword_sets is None:
        core_keyword_sets = _core_sets

    category_counts = {cat: 0 for cat in keyword_map}
    genes_in_category = {cat: [] for cat in keyword_map}
    categorised_genes: Set[str] = set()
    gene_categories: Dict[str, List[str]] = {}
    multi_category_genes = 0
    method_summary: Dict[str, Dict[str, int]] = {
        cat: {'keyword_only': 0, 'go_id_only': 0, 'both': 0}
        for cat in keyword_map
    }
    assignment_methods: Dict[str, Dict[str, str]] = {}

    # Core specificity tracking
    has_core = core_keyword_sets is not None and any(
        len(s) > 0 for s in core_keyword_sets.values())
    core_match: Dict[str, Dict[str, bool]] = {}
    core_summary: Dict[str, Dict[str, int]] = {
        cat: {'core': 0, 'extended': 0} for cat in keyword_map
    }

    use_go = expanded_go_sets is not None and any(
        len(s) > 0 for s in expanded_go_sets.values())
    available_go_cols = None
    if use_go:
        cols = go_columns or _DEFAULT_GO_COLUMNS
        available_go_cols = [c for c in cols if c in df.columns]
        if not available_go_cols:
            use_go = False

    for idx, row in df.iterrows():
        gene_id = str(row.get('gene_id', idx))
        all_text = ' '.join(str(v).lower() for v in row.values if isinstance(v, str))

        gene_go_ids = None
        if use_go:
            gene_go_ids = extract_go_ids_from_row(row, available_go_cols)

        matched = []
        gene_methods: Dict[str, str] = {}

        for category, keywords in keyword_map.items():
            kw_hit = match_keywords_to_annotation(
                all_text, keywords, require_word_boundary=True)

            go_hit = False
            if use_go and gene_go_ids and expanded_go_sets.get(category):
                go_hit = bool(gene_go_ids & expanded_go_sets[category])

            if (kw_hit or go_hit) and category not in matched:
                category_counts[category] += 1
                genes_in_category[category].append(gene_id)
                matched.append(category)

                if kw_hit and go_hit:
                    method = 'both'
                elif kw_hit:
                    method = 'keyword_only'
                else:
                    method = 'go_id_only'
                gene_methods[category] = method
                method_summary[category][method] += 1

                # Core specificity
                is_core = (
                    kw_hit is not None
                    and has_core
                    and kw_hit.lower().strip() in core_keyword_sets.get(category, frozenset())
                )
                if gene_id not in core_match:
                    core_match[gene_id] = {}
                core_match[gene_id][category] = is_core
                if is_core:
                    core_summary[category]['core'] += 1
                else:
                    core_summary[category]['extended'] += 1

                if not allow_multiple:
                    break

        if matched:
            categorised_genes.add(gene_id)
            gene_categories[gene_id] = matched
            assignment_methods[gene_id] = gene_methods
            if len(matched) > 1:
                multi_category_genes += 1

    total_categorised = len(categorised_genes)
    uncategorised = len(df) - total_categorised

    category_overlaps: Dict[Tuple[str, str], int] = {}
    if allow_multiple and multi_category_genes > 0:
        for cats in gene_categories.values():
            if len(cats) > 1:
                for i, c1 in enumerate(cats):
                    for c2 in cats[i + 1:]:
                        pair = tuple(sorted([c1, c2]))
                        category_overlaps[pair] = category_overlaps.get(pair, 0) + 1

    pct = round(total_categorised / len(df) * 100, 1) if len(df) > 0 else 0
    logger.info(f"Fallback: {total_categorised} categorised ({pct}%)")

    return {
        'category_counts': category_counts,
        'total_categorised': total_categorised,
        'uncategorised': uncategorised,
        'genes_in_category': genes_in_category,
        'multi_category_genes': multi_category_genes,
        'category_overlaps': category_overlaps,
        'gene_categories': gene_categories,
        'method_used': 'fallback',
        'method_summary': method_summary,
        'assignment_methods': assignment_methods,
        'core_match': core_match,
        'core_summary': core_summary,
    }


def categorise_and_fallback(
    df: pd.DataFrame,
    categories: Dict[str, Any],
    allow_multiple: bool = True,
    expanded_go_sets: Optional[Dict[str, Set[str]]] = None,
    go_columns: Optional[List[str]] = None,
    core_keyword_sets: Optional[Dict[str, frozenset]] = None,
) -> Dict[str, Any]:
    """Try standard categorisation, fall back to broader search if <5% hit."""
    results = categorise_genes(
        df, categories, allow_multiple,
        expanded_go_sets=expanded_go_sets, go_columns=go_columns,
        core_keyword_sets=core_keyword_sets)

    if results['total_categorised'] < 0.05 * len(df):
        logger.warning(f"Only {results['total_categorised']} genes categorised "
                       f"({results['total_categorised'] / max(len(df), 1) * 100:.1f}%), "
                       f"trying fallback...")
        fallback = fallback_categorisation(
            df, categories, allow_multiple,
            expanded_go_sets=expanded_go_sets, go_columns=go_columns,
            core_keyword_sets=core_keyword_sets)
        if fallback['total_categorised'] > results['total_categorised']:
            logger.info(f"Fallback improved: {fallback['total_categorised']} genes")
            return fallback
        else:
            logger.info("Fallback did not improve, keeping standard results")

    return results


# ---------------------------------------------------------------------------
# Category relationship analysis
# ---------------------------------------------------------------------------

def analyse_category_relationships(
    category_overlaps: Dict[Tuple[str, str], int],
    category_counts: Dict[str, int],
) -> Dict[str, Any]:
    """Analyse overlap/similarity between categories (Jaccard index)."""
    cats = sorted(category_counts.keys())
    overlap_matrix = pd.DataFrame(0, index=cats, columns=cats)
    normalised = pd.DataFrame(0.0, index=cats, columns=cats)

    for (c1, c2), count in category_overlaps.items():
        overlap_matrix.loc[c1, c2] = count
        overlap_matrix.loc[c2, c1] = count

    for c in cats:
        overlap_matrix.loc[c, c] = category_counts[c]

    for c1 in cats:
        for c2 in cats:
            if c1 == c2:
                normalised.loc[c1, c2] = 1.0
            else:
                union = category_counts[c1] + category_counts[c2] - overlap_matrix.loc[c1, c2]
                if union > 0:
                    normalised.loc[c1, c2] = overlap_matrix.loc[c1, c2] / union

    return {
        'overlap_matrix': overlap_matrix,
        'normalised_overlaps': normalised,
        'similarity_scores': normalised.copy(),
    }
