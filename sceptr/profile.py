"""
SCEPTR enrichment profiling orchestrator.

Runs the complete continuous enrichment profiling workflow:
1. Load and validate expression data
2. Load or auto-assign functional categories
3. Compute enrichment at discrete tiers
4. Compute continuous enrichment E_C(k) at every gene rank
5. Compute D_KL functional specialisation gradient
6. Run permutation-based global profile test
7. Generate interactive HTML report + TSV outputs
"""

import json
import logging
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from sceptr import categorisation, enrichment, continuous, io

logger = logging.getLogger("sceptr.profile")


def _resolve_category_json(category_set: str, analysis_type: str = "functional"):
    """Find the bundled category JSON for a given category set."""
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    cat_dir = os.path.join(pkg_dir, "categories", analysis_type)

    # Try standard naming patterns
    candidates = [
        os.path.join(cat_dir, f"{category_set}_{analysis_type}_categories.json"),
        os.path.join(cat_dir, f"{category_set}_categories.json"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path

    available = [f.replace(f"_{analysis_type}_categories.json", "")
                 for f in os.listdir(cat_dir) if f.endswith(".json")]
    raise FileNotFoundError(
        f"No {analysis_type} category set '{category_set}'. "
        f"Available: {', '.join(sorted(available))}"
    )


def _load_categories(category_set: Optional[str] = None,
                     custom_json: Optional[str] = None,
                     analysis_type: str = "functional"):
    """Load category definitions from JSON."""
    if custom_json:
        with open(custom_json) as f:
            return json.load(f)

    path = _resolve_category_json(category_set or "general", analysis_type)
    with open(path) as f:
        cats = json.load(f)
    logger.info(f"Loaded {len(cats)} categories from {category_set} "
                f"({analysis_type})")
    return cats


def run(expression_file: str,
        category_set: str = "general",
        custom_categories: Optional[str] = None,
        category_mapping_file: Optional[str] = None,
        output_dir: str = "sceptr_results",
        prefix: str = "sceptr",
        tiers: Optional[List[int]] = None,
        permutations: int = 1000,
        with_go: bool = False,
        analysis_type: str = "functional",
        quiet: bool = False) -> Dict:
    """Run SCEPTR enrichment profiling.

    Parameters
    ----------
    expression_file : str
        Path to expression data TSV/CSV. Must contain gene ID and expression
        value columns. May also contain protein_name, GO columns, or a
        categories column.
    category_set : str
        Name of built-in category set (e.g. 'bacteria', 'human_host').
    custom_categories : str, optional
        Path to custom category definitions JSON.
    category_mapping_file : str, optional
        Path to pre-computed gene-to-category mapping TSV.
    output_dir : str
        Output directory.
    prefix : str
        Output file prefix.
    tiers : list of int, optional
        Expression tiers for discrete enrichment (default: [50,100,250,500]).
    permutations : int
        Number of permutations for profile significance test.
    with_go : bool
        Enable GO hierarchy expansion (requires goatools + OBO file).
    analysis_type : str
        'functional' or 'cellular'.
    quiet : bool
        Suppress progress logging.

    Returns
    -------
    dict
        Results dictionary with keys: 'tier_results', 'continuous',
        'all_results', 'df_sorted', 'categories'.
    """
    if not quiet:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    if tiers is None:
        tiers = [50, 100, 250, 500]

    # Create output directories
    analysis_dir = os.path.join(output_dir, analysis_type)
    fig_dir = os.path.join(analysis_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    data_prefix = os.path.join(analysis_dir, prefix)
    fig_prefix = os.path.join(fig_dir, prefix)

    type_label = "BP_MF" if analysis_type == "functional" else "CC"

    # --- Step 1: Load expression data ---
    logger.info("Step 1: Loading expression data...")
    df = io.load_expression(expression_file)
    input_mode = io.detect_input_mode(df)
    logger.info(f"  Input mode: {input_mode} ({len(df)} genes)")

    # --- Step 2: Load categories and categorise genes ---
    logger.info("Step 2: Loading categories and assigning genes...")

    # If user provided a category mapping file, use that directly
    if category_mapping_file:
        ext_mapping = io.load_category_mapping(category_mapping_file)
        # Build category definitions from the mapping
        all_cats = set()
        for cats in ext_mapping.values():
            all_cats.update(cats)
        # Create a simple keyword map (each category maps to itself)
        keyword_map = {cat: [cat] for cat in sorted(all_cats)}
        anchor_map = {cat: [] for cat in keyword_map}
        core_keyword_sets = {cat: set() for cat in keyword_map}
        expanded_go_sets = None

        # Build categorisation results from external mapping
        category_counts = {cat: 0 for cat in keyword_map}
        genes_in_category = {cat: [] for cat in keyword_map}
        total_categorised = 0
        for _, row in df.iterrows():
            gid = str(row["sequence_id"])
            if gid in ext_mapping:
                total_categorised += 1
                for cat in ext_mapping[gid]:
                    if cat in category_counts:
                        category_counts[cat] += 1
                        genes_in_category[cat].append(gid)

        all_results = {
            "category_counts": category_counts,
            "total_categorised": total_categorised,
            "uncategorised": len(df) - total_categorised,
            "genes_in_category": genes_in_category,
            "multi_category_genes": 0,
        }
        use_external_mapping = True

    elif input_mode == "categorised":
        # Categories embedded in the expression file
        ext_mapping = io.extract_categories_from_column(df)
        all_cats = set()
        for cats in ext_mapping.values():
            all_cats.update(cats)
        keyword_map = {cat: [cat] for cat in sorted(all_cats)}
        anchor_map = {cat: [] for cat in keyword_map}
        core_keyword_sets = {cat: set() for cat in keyword_map}
        expanded_go_sets = None

        category_counts = {cat: 0 for cat in keyword_map}
        genes_in_category = {cat: [] for cat in keyword_map}
        total_categorised = 0
        for _, row in df.iterrows():
            gid = str(row["sequence_id"])
            if gid in ext_mapping:
                total_categorised += 1
                for cat in ext_mapping[gid]:
                    if cat in category_counts:
                        category_counts[cat] += 1
                        genes_in_category[cat].append(gid)

        all_results = {
            "category_counts": category_counts,
            "total_categorised": total_categorised,
            "uncategorised": len(df) - total_categorised,
            "genes_in_category": genes_in_category,
            "multi_category_genes": 0,
        }
        use_external_mapping = True

    else:
        # Use SCEPTR's built-in categorisation (keyword + optional GO)
        functional_categories = _load_categories(
            category_set, custom_categories, analysis_type)

        cat_format = categorisation.detect_category_format(functional_categories)
        keyword_map, anchor_map, core_keyword_sets = \
            categorisation.normalize_categories(functional_categories)
        logger.info(f"  Category format: {cat_format}")

        # GO expansion (optional)
        expanded_go_sets = None
        has_anchors = any(len(ids) > 0 for ids in anchor_map.values())

        if has_anchors and with_go:
            try:
                from sceptr import go_utils
                go_dag = go_utils.load_go_dag()
                expanded_go_sets = go_utils.expand_anchor_go_ids(anchor_map, go_dag)
                logger.info("  GO hierarchy expansion enabled")
            except (ImportError, Exception) as e:
                logger.warning(f"  GO expansion unavailable: {e}")
                logger.warning("  Falling back to keyword-only assignment")

        # Categorise all genes
        all_results = categorisation.categorise_genes(
            df, functional_categories, allow_multiple=True,
            expanded_go_sets=expanded_go_sets,
            core_keyword_sets=core_keyword_sets)
        use_external_mapping = False
        ext_mapping = None

    n_categorised = all_results["total_categorised"]
    n_cats = sum(1 for v in all_results["category_counts"].values() if v > 0)
    logger.info(f"  {n_categorised}/{len(df)} genes categorised into "
                f"{n_cats} categories")

    # --- Step 3: Discrete tier enrichment ---
    logger.info("Step 3: Computing discrete tier enrichment...")
    gene_subsets = {f"top_{t}": df.head(t) for t in tiers if t <= len(df)}

    tier_results = {}
    for tier_name, subset_df in gene_subsets.items():
        if use_external_mapping:
            # Count categories in this tier from external mapping
            tier_counts = {cat: 0 for cat in keyword_map}
            tier_genes = {cat: [] for cat in keyword_map}
            tier_categorised = 0
            mapping = ext_mapping
            for _, row in subset_df.iterrows():
                gid = str(row["sequence_id"])
                if gid in mapping:
                    tier_categorised += 1
                    for cat in mapping[gid]:
                        if cat in tier_counts:
                            tier_counts[cat] += 1
                            tier_genes[cat].append(gid)

            tier_cat_results = {
                "category_counts": tier_counts,
                "total_categorised": tier_categorised,
                "uncategorised": len(subset_df) - tier_categorised,
                "genes_in_category": tier_genes,
                "multi_category_genes": 0,
            }
        else:
            tier_cat_results = categorisation.categorise_genes(
                subset_df, functional_categories, allow_multiple=True,
                expanded_go_sets=expanded_go_sets,
                core_keyword_sets=core_keyword_sets)

        enrichment_stats = enrichment.calculate_enrichment(
            tier_cat_results["category_counts"],
            len(subset_df), len(df),
            all_results["category_counts"])

        enrichment_stats = enrichment.adjust_p_values(enrichment_stats)

        tier_results[tier_name] = {
            "counts": tier_cat_results["category_counts"],
            "categorised": tier_cat_results["total_categorised"],
            "uncategorised": tier_cat_results["uncategorised"],
            "enrichment": enrichment_stats,
            "genes_in_category": tier_cat_results["genes_in_category"],
            "multi_category_genes": tier_cat_results.get("multi_category_genes", 0),
            "core_pct": {},
            "core_summary": tier_cat_results.get("core_summary", {}),
        }

        sig = enrichment.get_significant_categories(
            enrichment_stats, p_threshold=0.05, use_adjusted_p=True)
        logger.info(f"  {tier_name}: {len(sig)} significant categories")

    # --- Step 4: Continuous enrichment ---
    logger.info("Step 4: Computing continuous enrichment profiles...")

    if use_external_mapping:
        # Build membership matrix from external mapping
        cat_names_sorted = sorted(keyword_map.keys())
        cat_names_sorted = [c for c in cat_names_sorted
                            if "ncharacteris" not in c.lower()]
        N = len(df)
        membership = np.zeros((N, len(cat_names_sorted)), dtype=bool)
        gene_ids = df["sequence_id"].astype(str).tolist()
        cat_to_idx = {c: i for i, c in enumerate(cat_names_sorted)}
        bg_counts = {c: 0 for c in cat_names_sorted}

        for row_idx, gid in enumerate(gene_ids):
            if gid in ext_mapping:
                for cat in ext_mapping[gid]:
                    if cat in cat_to_idx:
                        membership[row_idx, cat_to_idx[cat]] = True
                        bg_counts[cat] += 1

        cont_cat_names = cat_names_sorted
        cont_bg_counts = bg_counts
    else:
        membership, cont_cat_names, cont_bg_counts, _ = \
            continuous.build_membership_matrix(
                df, functional_categories, categorisation.categorise_genes,
                expanded_go_sets=expanded_go_sets,
                core_keyword_sets=core_keyword_sets)

    N = len(df)
    k_max = N // 2

    k_values, enrichment_matrix = continuous.compute_continuous_enrichment(
        membership, cont_cat_names, cont_bg_counts, N,
        k_min=10, k_max=k_max, step=5)
    logger.info(f"  Computed enrichment at {len(k_values)} points "
                f"(k={k_values[0]} to {k_values[-1]})")

    _, dkl_values = continuous.compute_continuous_dkl(
        membership, cont_cat_names, cont_bg_counts, N,
        k_min=10, k_max=k_max, step=5)

    shape_stats = continuous.classify_profile_shapes(
        k_values, enrichment_matrix, cont_cat_names)

    logger.info(f"  Running global profile test ({permutations} permutations)...")
    if not quiet:
        print(f"  Running permutation test ({permutations} permutations)...", flush=True)
    _, profile_stats = continuous.permutation_global_test(
        membership, cont_cat_names, cont_bg_counts, N,
        k_min=10, k_max=k_max, step=5,
        n_permutations=permutations)

    n_sig = sum(1 for s in profile_stats.values() if s["supremum_p"] < 0.05)
    logger.info(f"  {n_sig}/{len(cont_cat_names)} categories with "
                f"significant profiles")

    cont_results = {
        "k_values": k_values,
        "enrichment_matrix": enrichment_matrix,
        "dkl_values": dkl_values,
        "profile_stats": profile_stats,
        "shape_stats": shape_stats,
        "cat_names": cont_cat_names,
    }

    # --- Step 5: Save outputs ---
    logger.info("Step 5: Saving results...")

    # Save enrichment TSV
    enrichment.save_enrichment_tsv(
        tier_results, list(gene_subsets.keys()), data_prefix, type_label)

    # Save continuous enrichment TSV
    continuous.save_continuous_enrichment_tsv(
        k_values, enrichment_matrix, cont_cat_names, data_prefix, type_label)

    # Save profile test TSV
    continuous.save_profile_test_tsv(
        profile_stats, cont_cat_names, data_prefix, type_label)

    # Save DKL TSV
    continuous.save_continuous_dkl_tsv(
        k_values, dkl_values, data_prefix, type_label)

    # Save shape classification
    continuous.save_shape_classification_tsv(
        shape_stats, data_prefix, type_label)

    # --- Step 6: Generate interactive HTML report ---
    logger.info("Step 6: Generating interactive report...")
    try:
        from sceptr import report as interactive_report
        interactive_report.generate_interactive_report(
            tier_results, all_results, data_prefix, len(df),
            cont_results=cont_results,
            chart_type=type_label,
            report_title="Functional Profiling" if analysis_type == "functional"
                         else "Cellular Component Profiling",
            description="biological process and molecular function"
                        if analysis_type == "functional"
                        else "cellular component localisation",
            figures_dir=fig_dir,
            df_sorted=df)
        logger.info(f"  Report: {data_prefix}_{type_label}_report.html")
    except Exception as e:
        logger.error(f"Report generation failed: {e}")
        import traceback
        traceback.print_exc()

    # --- Step 7: Static figures (optional, requires matplotlib) ---
    try:
        from sceptr._visualisation import create_figures
        create_figures(k_values, enrichment_matrix, dkl_values,
                       cont_cat_names, profile_stats, tier_results,
                       fig_prefix, type_label, tiers)
    except ImportError:
        pass  # matplotlib not installed, skip static figures
    except Exception as e:
        logger.warning(f"Static figure generation failed: {e}")

    logger.info(f"Done! Results saved to {output_dir}/")

    return {
        "tier_results": tier_results,
        "continuous": cont_results,
        "all_results": all_results,
        "df_sorted": df,
        "categories": keyword_map,
    }
