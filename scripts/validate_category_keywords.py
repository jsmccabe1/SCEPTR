#!/usr/bin/env python3
"""
Multi-layer provenance audit for SCEPTR category keywords.

Validates every keyword in every category JSON file against four independent
authoritative sources, producing a per-keyword provenance tag and summary
statistics suitable for manuscript methods and supplementary tables.

Layers:
  1. GO Ontology    - keyword matches GO term name or synonym (anchor + descendants)
  2. UniProt FASTA  - keyword appears in Swiss-Prot protein descriptions
  3. GO Slim        - anchor GO IDs overlap with expert-curated organism GO slims
  4. Provenance tag - explicit source label per keyword for supplementary table

Usage:
    python3 scripts/validate_category_keywords.py \
        --obo data/go/go-basic.obo \
        --categories modules/explot/categories \
        --uniprot-fasta data/uniprot/uniprot_sprot.fasta \
        --depth 3 \
        --output scripts/category_keyword_audit.tsv

The --uniprot-fasta flag enables Layer 2 (empirical hit-rate validation).
Without it, only Layers 1 and 3 run.
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


# ──────────────────────────────────────────────────────────────────────
# Layer 1: GO OBO parsing
# ──────────────────────────────────────────────────────────────────────

def parse_obo(obo_path):
    """Parse GO OBO file into term records, child graph, and slim sets."""
    terms = {}
    current = None
    slims = defaultdict(set)  # slim_name -> set of GO IDs

    with open(obo_path) as f:
        for line in f:
            line = line.rstrip()
            if line == "[Term]":
                current = {
                    "name": "",
                    "synonyms": [],
                    "is_obsolete": False,
                    "namespace": "",
                    "parents": [],
                    "subsets": [],
                }
            elif line == "[Typedef]":
                current = None
            elif current is not None:
                if line.startswith("id: GO:"):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("namespace: "):
                    current["namespace"] = line[11:]
                elif line.startswith("synonym: "):
                    m = re.match(r'synonym:\s+"(.+?)"\s+\w+', line)
                    if m:
                        current["synonyms"].append(m.group(1))
                elif line.startswith("is_a: GO:"):
                    parent_id = line.split("!")[0].strip().replace("is_a: ", "")
                    current["parents"].append(parent_id)
                elif line.startswith("subset: goslim_"):
                    slim_name = line.split()[1]
                    current["subsets"].append(slim_name)
                elif line == "is_obsolete: true":
                    current["is_obsolete"] = True
                elif line == "":
                    if "id" in current:
                        terms[current["id"]] = current
                        for s in current["subsets"]:
                            slims[s].add(current["id"])
                    current = None

    if current and "id" in current:
        terms[current["id"]] = current
        for s in current.get("subsets", []):
            slims[s].add(current["id"])

    children = defaultdict(list)
    for go_id, term in terms.items():
        for parent in term["parents"]:
            children[parent].append(go_id)

    return terms, children, slims


def get_descendants(go_id, children_map, terms, max_depth):
    """BFS to collect descendant GO IDs up to max_depth."""
    visited = set()
    queue = [(go_id, 0)]
    descendants = []
    while queue:
        current_id, depth = queue.pop(0)
        if current_id in visited:
            continue
        visited.add(current_id)
        if current_id != go_id:
            descendants.append(current_id)
        if depth < max_depth:
            for child_id in children_map.get(current_id, []):
                if child_id not in visited and child_id in terms:
                    queue.append((child_id, depth + 1))
    return descendants


def get_go_vocabulary(go_ids, terms, children_map, max_depth):
    """Collect all term names and synonyms from anchors and descendants."""
    vocabulary = set()
    vocab_sources = {}

    for go_id in go_ids:
        if go_id not in terms:
            continue
        all_ids = [go_id] + get_descendants(go_id, children_map, terms, max_depth)
        for gid in all_ids:
            if gid not in terms or terms[gid]["is_obsolete"]:
                continue
            term = terms[gid]
            name_lower = term["name"].lower()
            vocabulary.add(name_lower)
            if name_lower not in vocab_sources:
                vocab_sources[name_lower] = (gid, "GO:name")
            for syn in term["synonyms"]:
                syn_lower = syn.lower()
                vocabulary.add(syn_lower)
                if syn_lower not in vocab_sources:
                    vocab_sources[syn_lower] = (gid, "GO:synonym")

    return vocabulary, vocab_sources


def keyword_matches_vocabulary(keyword, vocabulary):
    """Multi-strategy keyword-to-GO matching."""
    kw_lower = keyword.lower()
    if kw_lower in vocabulary:
        return "exact", kw_lower
    for term in vocabulary:
        if kw_lower in term:
            return "substring_of_go", term
    for term in vocabulary:
        if len(term) >= 4 and term in kw_lower:
            return "go_substring_of_kw", term
    return None, None


# ──────────────────────────────────────────────────────────────────────
# Layer 2: UniProt FASTA empirical hit-rate
# ──────────────────────────────────────────────────────────────────────

def load_uniprot_descriptions(fasta_path):
    """Extract all protein descriptions from Swiss-Prot FASTA headers."""
    descriptions = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                # Format: >sp|ACCESSION|ENTRY_NAME Description OS=... OX=... GN=... PE=... SV=...
                # Extract description between entry name and OS=
                header = line[1:].strip()
                # Remove the sp|...|ENTRY_NAME prefix
                parts = header.split(" ", 1)
                if len(parts) > 1:
                    desc_part = parts[1]
                    # Cut at OS= if present
                    os_idx = desc_part.find(" OS=")
                    if os_idx > 0:
                        desc_part = desc_part[:os_idx]
                    descriptions.append(desc_part.lower())
    return descriptions


def count_keyword_hits(keyword, descriptions):
    """Count how many UniProt descriptions contain this keyword."""
    kw_lower = keyword.lower()
    # Use word boundary matching for short keywords to avoid false positives
    if len(kw_lower) <= 3:
        pattern = re.compile(r'\b' + re.escape(kw_lower) + r'\b')
        return sum(1 for d in descriptions if pattern.search(d))
    else:
        return sum(1 for d in descriptions if kw_lower in d)


def batch_keyword_hits(keywords_set, descriptions):
    """Count hits for all keywords at once (more efficient)."""
    results = {}
    # Pre-compile patterns for short keywords
    patterns = {}
    for kw in keywords_set:
        kw_lower = kw.lower()
        if len(kw_lower) <= 3:
            patterns[kw] = re.compile(r'\b' + re.escape(kw_lower) + r'\b')
        else:
            patterns[kw] = None  # Use simple 'in' check

    for kw in keywords_set:
        kw_lower = kw.lower()
        pat = patterns[kw]
        if pat is not None:
            results[kw] = sum(1 for d in descriptions if pat.search(d))
        else:
            results[kw] = sum(1 for d in descriptions if kw_lower in d)
    return results


# ──────────────────────────────────────────────────────────────────────
# Layer 3: GO Slim alignment
# ──────────────────────────────────────────────────────────────────────

# Mapping from SCEPTR category file prefixes to relevant GO slims
SLIM_MAP = {
    "fungi": ["goslim_candida", "goslim_yeast", "goslim_pombe"],
    "cancer": ["goslim_generic", "goslim_chembl"],
    "insect": ["goslim_drosophila", "goslim_flybase_ribbon"],
    "general": ["goslim_generic"],
    "parasite_protozoan": ["goslim_generic"],
    "parasite_metazoan": ["goslim_generic"],
    "protist_dinoflagellate": ["goslim_generic"],
    "model_organism": ["goslim_generic", "goslim_mouse", "goslim_yeast", "goslim_drosophila"],
    "plant": ["goslim_plant"],
    "bacteria": ["goslim_prokaryote"],
    "bacteria_gram_negative": ["goslim_prokaryote"],
    "bacteria_gram_positive": ["goslim_prokaryote"],
    "virus": ["goslim_virus"],
    "vertebrate_host": ["goslim_generic", "goslim_chembl"],
    "vertebrate_host_hallmark": ["goslim_generic", "goslim_chembl"],
}


def get_slim_set_for_file(json_filename, slims):
    """Get the union of relevant GO slim term IDs for a category file."""
    # Extract prefix: fungi_functional_categories.json -> fungi
    basename = os.path.basename(json_filename)
    for prefix in sorted(SLIM_MAP.keys(), key=len, reverse=True):
        if basename.startswith(prefix):
            slim_names = SLIM_MAP[prefix]
            combined = set()
            for s in slim_names:
                combined |= slims.get(s, set())
            return combined, slim_names
    return set(), []


def check_slim_coverage(anchor_go_ids, slim_go_ids, terms, children_map):
    """
    Check if anchor GO IDs (or their ancestors/descendants) overlap with
    the GO slim. Returns overlap details.
    """
    if not slim_go_ids or not anchor_go_ids:
        return {"direct": [], "via_ancestor": [], "via_descendant": [], "none": list(anchor_go_ids)}

    direct = []
    via_ancestor = []
    via_descendant = []
    no_match = []

    for gid in anchor_go_ids:
        if gid not in terms:
            no_match.append(gid)
            continue

        # Direct overlap
        if gid in slim_go_ids:
            direct.append(gid)
            continue

        # Check if any ancestor is in the slim
        found_ancestor = False
        visited = set()
        queue = [gid]
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            if current != gid and current in slim_go_ids:
                via_ancestor.append((gid, current))
                found_ancestor = True
                break
            for parent in terms.get(current, {}).get("parents", []):
                if parent not in visited:
                    queue.append(parent)

        if found_ancestor:
            continue

        # Check if any descendant is in the slim
        descendants = get_descendants(gid, children_map, terms, 5)
        slim_descendants = [d for d in descendants if d in slim_go_ids]
        if slim_descendants:
            via_descendant.append((gid, slim_descendants[0]))
        else:
            no_match.append(gid)

    return {"direct": direct, "via_ancestor": via_ancestor,
            "via_descendant": via_descendant, "none": no_match}


# ──────────────────────────────────────────────────────────────────────
# Layer 4: Provenance tagging
# ──────────────────────────────────────────────────────────────────────

def assign_provenance(keyword, go_match_type, go_matched_term, go_source_id,
                      uniprot_hits, vocab_sources):
    """
    Assign a provenance tag to a keyword based on all available evidence.

    Tags (in priority order):
      GO:name        - exact match to GO term name
      GO:synonym     - exact match to GO synonym
      GO:derived     - substring match to GO vocabulary
      UniProt:N      - found in N Swiss-Prot protein descriptions
      literature     - needs manual citation (no automated source found)
    """
    tags = []

    # GO provenance
    if go_match_type == "exact":
        source_type = vocab_sources.get(go_matched_term, ("", ""))[1]
        if source_type:
            tags.append(source_type)
        else:
            tags.append("GO:name")
    elif go_match_type in ("substring_of_go", "go_substring_of_kw"):
        tags.append("GO:derived")

    # UniProt provenance
    if uniprot_hits is not None and uniprot_hits > 0:
        tags.append(f"UniProt:{uniprot_hits}")

    if not tags:
        tags.append("literature")

    return "; ".join(tags)


# ──────────────────────────────────────────────────────────────────────
# Audit engine
# ──────────────────────────────────────────────────────────────────────

def audit_category_file(json_path, terms, children_map, slims, max_depth,
                        uniprot_descriptions=None):
    """Full multi-layer audit of a single category JSON file."""
    with open(json_path) as f:
        categories = json.load(f)

    # Layer 3 prep: get relevant slim
    slim_go_ids, slim_names = get_slim_set_for_file(str(json_path), slims)

    # Collect all unique keywords for batch UniProt lookup
    all_keywords = set()
    for cat_data in categories.values():
        all_keywords.update(cat_data.get("keywords", []))
        all_keywords.update(cat_data.get("core_keywords", []))

    # Layer 2: batch UniProt hit counting
    uniprot_hits_map = None
    if uniprot_descriptions is not None:
        uniprot_hits_map = batch_keyword_hits(all_keywords, uniprot_descriptions)

    results = []
    file_total = 0
    file_go_backed = 0
    file_uniprot_backed = 0
    file_any_backed = 0

    for cat_name, cat_data in categories.items():
        keywords = cat_data.get("keywords", [])
        anchor_ids = cat_data.get("anchor_go_ids", [])
        core_keywords = cat_data.get("core_keywords", [])

        # Layer 1: GO vocabulary
        obsolete_anchors = [g for g in anchor_ids if g in terms and terms[g]["is_obsolete"]]
        missing_anchors = [g for g in anchor_ids if g not in terms]
        vocabulary, vocab_sources = get_go_vocabulary(anchor_ids, terms, children_map, max_depth)

        # Layer 3: slim coverage
        slim_coverage = check_slim_coverage(anchor_ids, slim_go_ids, terms, children_map)

        # Audit each keyword
        keyword_results = []
        go_backed = 0
        uniprot_backed = 0
        any_backed = 0

        for kw in keywords:
            # Layer 1
            match_type, matched_term = keyword_matches_vocabulary(kw, vocabulary)
            source_go_id = ""
            if match_type and matched_term in vocab_sources:
                source_go_id = vocab_sources[matched_term][0]

            is_go_backed = match_type is not None

            # Layer 2
            hits = uniprot_hits_map.get(kw) if uniprot_hits_map else None
            is_uniprot_backed = hits is not None and hits > 0

            # Provenance tag
            provenance = assign_provenance(
                kw, match_type, matched_term, source_go_id,
                hits, vocab_sources
            )

            if is_go_backed:
                go_backed += 1
            if is_uniprot_backed:
                uniprot_backed += 1
            if is_go_backed or is_uniprot_backed:
                any_backed += 1

            keyword_results.append({
                "keyword": kw,
                "go_status": "GO-backed" if is_go_backed else "not-in-GO",
                "go_match_type": match_type or "",
                "go_matched_term": matched_term or "",
                "go_source_id": source_go_id,
                "uniprot_hits": hits,
                "provenance": provenance,
                "is_core": kw in core_keywords,
            })

        file_total += len(keywords)
        file_go_backed += go_backed
        file_uniprot_backed += uniprot_backed
        file_any_backed += any_backed

        results.append({
            "category": cat_name,
            "n_keywords": len(keywords),
            "n_go_backed": go_backed,
            "n_uniprot_backed": uniprot_backed,
            "n_any_backed": any_backed,
            "n_literature_only": len(keywords) - any_backed,
            "go_coverage_pct": (go_backed / len(keywords) * 100) if keywords else 0,
            "uniprot_coverage_pct": (uniprot_backed / len(keywords) * 100) if keywords else 0,
            "combined_coverage_pct": (any_backed / len(keywords) * 100) if keywords else 0,
            "anchor_go_ids": anchor_ids,
            "obsolete_anchors": obsolete_anchors,
            "missing_anchors": missing_anchors,
            "go_vocab_size": len(vocabulary),
            "slim_names": slim_names,
            "slim_coverage": slim_coverage,
            "keyword_results": keyword_results,
        })

    return results, file_total, file_go_backed, file_uniprot_backed, file_any_backed


# ──────────────────────────────────────────────────────────────────────
# Reporting
# ──────────────────────────────────────────────────────────────────────

def print_report(json_file, results, file_total, file_go_backed,
                 file_uniprot_backed, file_any_backed, verbose=False):
    """Print human-readable report for one file."""
    go_pct = (file_go_backed / file_total * 100) if file_total > 0 else 0
    up_pct = (file_uniprot_backed / file_total * 100) if file_total > 0 else 0
    any_pct = (file_any_backed / file_total * 100) if file_total > 0 else 0
    lit_only = file_total - file_any_backed

    print(f"{'='*80}")
    print(f"FILE: {json_file}")
    print(f"  Layer 1 (GO):      {file_go_backed}/{file_total} ({go_pct:.1f}%)")
    if file_uniprot_backed > 0 or any(r["keyword_results"][0].get("uniprot_hits") is not None
                                       for r in results if r["keyword_results"]):
        print(f"  Layer 2 (UniProt): {file_uniprot_backed}/{file_total} ({up_pct:.1f}%)")
    print(f"  Combined:          {file_any_backed}/{file_total} ({any_pct:.1f}%)")
    if lit_only > 0:
        print(f"  Literature-only:   {lit_only}")
    print()

    for cat in results:
        # Status icon
        if cat["combined_coverage_pct"] >= 80:
            icon = "OK"
        elif cat["combined_coverage_pct"] >= 50:
            icon = "ok"
        else:
            icon = "REVIEW"

        parts = [f"GO:{cat['n_go_backed']}"]
        if cat["n_uniprot_backed"] > 0:
            parts.append(f"UP:{cat['n_uniprot_backed']}")
        parts.append(f"combined:{cat['n_any_backed']}/{cat['n_keywords']}")
        parts.append(f"({cat['combined_coverage_pct']:.0f}%)")

        print(f"  [{icon}] {cat['category']}: {' | '.join(parts)}")

        # Warnings
        if cat["obsolete_anchors"]:
            print(f"    WARNING: Obsolete anchor GO IDs: {', '.join(cat['obsolete_anchors'])}")
        if cat["missing_anchors"]:
            print(f"    WARNING: Missing anchor GO IDs: {', '.join(cat['missing_anchors'])}")

        # Slim coverage
        sc = cat["slim_coverage"]
        if cat["slim_names"]:
            n_anchors = len(cat["anchor_go_ids"])
            n_slim_covered = len(sc["direct"]) + len(sc["via_ancestor"]) + len(sc["via_descendant"])
            if n_anchors > 0:
                print(f"    GO slim ({', '.join(cat['slim_names'])}): "
                      f"{n_slim_covered}/{n_anchors} anchors covered "
                      f"(direct:{len(sc['direct'])}, ancestor:{len(sc['via_ancestor'])}, "
                      f"descendant:{len(sc['via_descendant'])})")

        # Literature-only keywords
        lit_only = [kr for kr in cat["keyword_results"]
                    if kr["provenance"] == "literature"]
        if lit_only:
            print(f"    Literature-only keywords ({len(lit_only)}):")
            for kr in lit_only:
                print(f"      - \"{kr['keyword']}\"")

        if verbose:
            # Show all keywords with provenance
            for kr in cat["keyword_results"]:
                core_tag = " [CORE]" if kr["is_core"] else ""
                hits_tag = f" (UP:{kr['uniprot_hits']})" if kr["uniprot_hits"] else ""
                print(f"      {kr['provenance']:30s} | \"{kr['keyword']}\"{core_tag}{hits_tag}")

    print()


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Multi-layer provenance audit for SCEPTR category keywords"
    )
    parser.add_argument("--obo", required=True, help="Path to go-basic.obo")
    parser.add_argument("--categories", required=True,
                        help="Path to categories directory (or single JSON)")
    parser.add_argument("--uniprot-fasta", default=None,
                        help="Path to Swiss-Prot FASTA for empirical hit-rate (Layer 2)")
    parser.add_argument("--depth", type=int, default=3,
                        help="Max GO descendant depth (default: 3)")
    parser.add_argument("--output", default=None,
                        help="Output TSV path for per-keyword audit table")
    parser.add_argument("--verbose", action="store_true",
                        help="Print per-keyword provenance details")
    parser.add_argument("--files", nargs="*", default=None,
                        help="Specific JSON files to audit (default: all)")
    args = parser.parse_args()

    # Parse OBO
    print(f"Parsing GO OBO: {args.obo}", file=sys.stderr)
    terms, children_map, slims = parse_obo(args.obo)
    n_active = sum(1 for t in terms.values() if not t["is_obsolete"])
    n_obsolete = sum(1 for t in terms.values() if t["is_obsolete"])
    n_slims = len(slims)
    print(f"  {n_active} active terms, {n_obsolete} obsolete, {n_slims} GO slims loaded",
          file=sys.stderr)

    # Load UniProt descriptions (Layer 2)
    uniprot_descriptions = None
    if args.uniprot_fasta:
        print(f"Loading UniProt FASTA: {args.uniprot_fasta}", file=sys.stderr)
        uniprot_descriptions = load_uniprot_descriptions(args.uniprot_fasta)
        print(f"  {len(uniprot_descriptions)} protein descriptions loaded", file=sys.stderr)

    # Find category files
    cat_path = Path(args.categories)
    if cat_path.is_file():
        json_files = [cat_path]
    elif args.files:
        json_files = [Path(f) for f in args.files]
    else:
        json_files = sorted(cat_path.rglob("*.json"))

    if not json_files:
        print("No category JSON files found.", file=sys.stderr)
        sys.exit(1)

    print(f"Auditing {len(json_files)} files (depth={args.depth})\n", file=sys.stderr)

    # TSV output
    tsv_rows = []
    grand_total = 0
    grand_go = 0
    grand_uniprot = 0
    grand_any = 0

    for json_file in json_files:
        rel_path = json_file.relative_to(cat_path) if cat_path.is_dir() else json_file.name

        results, ft, fg, fu, fa = audit_category_file(
            json_file, terms, children_map, slims, args.depth,
            uniprot_descriptions
        )
        grand_total += ft
        grand_go += fg
        grand_uniprot += fu
        grand_any += fa

        print_report(str(rel_path), results, ft, fg, fu, fa, verbose=args.verbose)

        # Build TSV rows
        for cat in results:
            sc = cat["slim_coverage"]
            n_slim = len(sc["direct"]) + len(sc["via_ancestor"]) + len(sc["via_descendant"])
            for kr in cat["keyword_results"]:
                tsv_rows.append({
                    "file": str(rel_path),
                    "category": cat["category"],
                    "keyword": kr["keyword"],
                    "is_core_keyword": kr["is_core"],
                    "go_status": kr["go_status"],
                    "go_match_type": kr["go_match_type"],
                    "go_matched_term": kr["go_matched_term"],
                    "go_source_id": kr["go_source_id"],
                    "uniprot_hits": kr["uniprot_hits"] if kr["uniprot_hits"] is not None else "",
                    "provenance": kr["provenance"],
                    "slim_names": "; ".join(cat["slim_names"]),
                    "anchor_slim_coverage": f"{n_slim}/{len(cat['anchor_go_ids'])}",
                })

    # Grand summary
    go_pct = (grand_go / grand_total * 100) if grand_total else 0
    up_pct = (grand_uniprot / grand_total * 100) if grand_total else 0
    any_pct = (grand_any / grand_total * 100) if grand_total else 0

    print(f"{'='*80}")
    print(f"GRAND TOTAL: {grand_total} keywords across {len(json_files)} files")
    print(f"  Layer 1 (GO ontology):     {grand_go}/{grand_total} ({go_pct:.1f}%)")
    if uniprot_descriptions:
        print(f"  Layer 2 (UniProt empirical): {grand_uniprot}/{grand_total} ({up_pct:.1f}%)")
    print(f"  Combined (any source):     {grand_any}/{grand_total} ({any_pct:.1f}%)")
    print(f"  Literature-only:           {grand_total - grand_any}")
    print()

    # Write TSV
    if args.output:
        header = ["file", "category", "keyword", "is_core_keyword",
                  "go_status", "go_match_type", "go_matched_term", "go_source_id",
                  "uniprot_hits", "provenance", "slim_names", "anchor_slim_coverage"]
        with open(args.output, "w") as f:
            f.write("\t".join(header) + "\n")
            for row in tsv_rows:
                f.write("\t".join(str(row[h]) for h in header) + "\n")
        print(f"Audit table written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
