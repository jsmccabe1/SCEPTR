"""
Flexible input parsing for SCEPTR standalone mode.

Supports multiple input formats:
1. SCEPTR integrated_annotations_expression.tsv (full format)
2. Simple expression + categories (gene_id, TPM, categories)
3. Expression + protein descriptions (auto-categorisation)
4. Separate expression and category mapping files
"""

import logging
import os
from typing import Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger("sceptr.io")

# Column name aliases for flexible matching
_TPM_ALIASES = ["TPM", "tpm", "expression", "Expression", "FPKM", "fpkm",
                "CPM", "cpm", "counts", "Counts", "value", "Value",
                "abundance", "Abundance"]
_GENE_ID_ALIASES = ["sequence_id", "gene_id", "gene", "Gene", "transcript_id",
                    "transcript", "Transcript", "ID", "id", "Name", "name",
                    "feature_id"]
_PROTEIN_ALIASES = ["protein_name", "Protein_Name", "description",
                    "Description", "annotation", "Annotation",
                    "protein_description"]
_CATEGORY_ALIASES = ["categories", "Categories", "category", "Category",
                     "functional_categories"]


def _find_column(df: pd.DataFrame, aliases: List[str],
                 required: bool = True, label: str = "") -> Optional[str]:
    """Find a column matching any of the given aliases."""
    for alias in aliases:
        if alias in df.columns:
            return alias
    if required:
        raise ValueError(
            f"Could not find {label} column. Tried: {aliases[:5]}... "
            f"Available columns: {list(df.columns)}"
        )
    return None


def load_expression(path: str) -> pd.DataFrame:
    """Load and validate an expression table.

    Accepts TSV or CSV. Requires at minimum a gene ID column and an
    expression value column. Returns a DataFrame sorted by expression
    (descending) with standardised column names.
    """
    sep = "\t" if path.endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(path, sep=sep)

    if df.shape[1] == 1 and "\t" in open(path).readline():
        df = pd.read_csv(path, sep="\t")

    gene_col = _find_column(df, _GENE_ID_ALIASES, label="gene ID")
    tpm_col = _find_column(df, _TPM_ALIASES, label="expression value")

    if gene_col != "sequence_id":
        df = df.rename(columns={gene_col: "sequence_id"})
    if tpm_col != "TPM":
        df = df.rename(columns={tpm_col: "TPM"})

    df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce")
    df = df.dropna(subset=["TPM"])
    df = df.sort_values("TPM", ascending=False).reset_index(drop=True)

    logger.info(f"Loaded {len(df)} genes from {os.path.basename(path)}")
    return df


def load_category_mapping(path: str) -> Dict[str, List[str]]:
    """Load a category mapping file.

    Expected format: TSV/CSV with gene_id and category columns.
    One row per gene-category assignment, OR a single row per gene with
    semicolon-delimited categories.

    Returns: {gene_id: [category1, category2, ...]}
    """
    sep = "\t" if path.endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(path, sep=sep)

    gene_col = _find_column(df, _GENE_ID_ALIASES, label="gene ID")
    cat_col = _find_column(df, _CATEGORY_ALIASES + ["category_name"],
                           label="category")

    mapping = {}
    for _, row in df.iterrows():
        gid = str(row[gene_col])
        cats_raw = str(row[cat_col])
        cats = [c.strip() for c in cats_raw.split(";") if c.strip()]
        if gid in mapping:
            mapping[gid].extend(cats)
        else:
            mapping[gid] = cats

    n_genes = len(mapping)
    n_cats = len(set(c for cats in mapping.values() for c in cats))
    logger.info(f"Loaded category mapping: {n_genes} genes, {n_cats} categories")
    return mapping


def detect_input_mode(df: pd.DataFrame) -> str:
    """Detect the input mode based on available columns.

    Returns one of:
        'full'        - SCEPTR integrated format (has protein_name + GO columns)
        'categorised' - Has a categories column (pre-assigned)
        'annotated'   - Has protein descriptions (needs auto-categorisation)
        'minimal'     - Only gene ID + expression (needs external categories)
    """
    has_protein = _find_column(df, _PROTEIN_ALIASES, required=False) is not None
    has_categories = _find_column(df, _CATEGORY_ALIASES, required=False) is not None
    has_go = any(c for c in df.columns if "GO_" in c or "go_" in c)

    if has_protein and has_go:
        return "full"
    if has_categories:
        return "categorised"
    if has_protein:
        return "annotated"
    return "minimal"


def apply_category_mapping(df: pd.DataFrame,
                           mapping: Dict[str, List[str]]) -> pd.DataFrame:
    """Add category assignments from an external mapping to the DataFrame."""
    df["_mapped_categories"] = df["sequence_id"].map(
        lambda x: ";".join(mapping.get(str(x), []))
    )
    n_mapped = (df["_mapped_categories"] != "").sum()
    logger.info(f"Mapped categories to {n_mapped}/{len(df)} genes")
    return df


def extract_categories_from_column(df: pd.DataFrame) -> Dict[str, List[str]]:
    """Extract category assignments from a 'categories' column in the DataFrame."""
    cat_col = _find_column(df, _CATEGORY_ALIASES, required=True,
                           label="category")
    mapping = {}
    for _, row in df.iterrows():
        gid = str(row["sequence_id"])
        raw = str(row[cat_col]) if pd.notna(row[cat_col]) else ""
        cats = [c.strip() for c in raw.split(";") if c.strip()]
        if cats:
            mapping[gid] = cats
    return mapping
