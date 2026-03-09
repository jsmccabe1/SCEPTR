"""
SCEPTR command-line interface.

Usage:
    sceptr profile --expression data.tsv --category-set bacteria -o results/
    sceptr profile --expression data.tsv --categories mapping.tsv -o results/
    sceptr categories --list
"""

import argparse
import json
import os
import sys

from sceptr import __version__


def _list_category_sets():
    """List all available built-in category sets."""
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    cat_dir = os.path.join(pkg_dir, "categories", "functional")
    sets = []
    for f in sorted(os.listdir(cat_dir)):
        if f.endswith("_functional_categories.json"):
            name = f.replace("_functional_categories.json", "")
            path = os.path.join(cat_dir, f)
            with open(path) as fh:
                cats = json.load(fh)
            sets.append((name, len(cats)))
    return sets


def cmd_profile(args):
    """Run enrichment profiling."""
    from sceptr.profile import run

    run(
        expression_file=args.expression,
        category_set=args.category_set,
        custom_categories=args.custom_categories,
        category_mapping_file=args.categories,
        output_dir=args.output,
        prefix=args.prefix,
        tiers=[int(t) for t in args.tiers.split(",")],
        permutations=args.permutations,
        with_go=args.with_go,
        analysis_type=args.analysis_type,
        quiet=args.quiet,
    )


def cmd_categories(args):
    """List or inspect category sets."""
    if args.list:
        print(f"\nAvailable category sets ({len(_list_category_sets())}):\n")
        for name, n_cats in _list_category_sets():
            print(f"  {name:<30s} {n_cats} categories")
        print(f"\nUsage: sceptr profile --expression data.tsv "
              f"--category-set <name>\n")
    elif args.show:
        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(pkg_dir, "categories", "functional",
                            f"{args.show}_functional_categories.json")
        if not os.path.exists(path):
            print(f"Category set '{args.show}' not found.")
            sys.exit(1)
        with open(path) as f:
            cats = json.load(f)
        print(f"\n{args.show} ({len(cats)} categories):\n")
        for name, defn in cats.items():
            if isinstance(defn, dict):
                n_kw = len(defn.get("keywords", []))
                n_go = len(defn.get("anchor_go_ids", []))
                print(f"  {name:<40s} {n_kw} keywords, {n_go} GO anchors")
            else:
                print(f"  {name:<40s} {len(defn)} keywords")
        print()


def main():
    parser = argparse.ArgumentParser(
        prog="sceptr",
        description="SCEPTR: Statistical Characterisation of Expression "
                    "Profiles in Transcriptomes",
    )
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # --- profile ---
    p_profile = subparsers.add_parser(
        "profile",
        help="Run continuous enrichment profiling",
        description="Compute continuous enrichment profiles, permutation "
                    "significance, and DKL gradient from expression data.",
    )
    p_profile.add_argument(
        "--expression", "-e", required=True,
        help="Path to expression data (TSV/CSV). Must have gene ID + "
             "expression columns. May also contain protein_name, GO terms, "
             "or a categories column.")
    p_profile.add_argument(
        "--category-set", "-s", default="general",
        help="Built-in category set name (default: general). "
             "Run 'sceptr categories --list' to see options.")
    p_profile.add_argument(
        "--categories", "-c",
        help="Path to external category mapping file (TSV/CSV). "
             "Columns: gene_id, category. Overrides --category-set.")
    p_profile.add_argument(
        "--custom-categories",
        help="Path to custom category definitions JSON (SCEPTR format).")
    p_profile.add_argument(
        "--output", "-o", default="sceptr_results",
        help="Output directory (default: sceptr_results)")
    p_profile.add_argument(
        "--prefix", "-p", default="sceptr",
        help="Output file prefix (default: sceptr)")
    p_profile.add_argument(
        "--tiers", default="50,100,250,500",
        help="Comma-separated expression tiers (default: 50,100,250,500)")
    p_profile.add_argument(
        "--permutations", type=int, default=1000,
        help="Permutations for profile significance test (default: 1000)")
    p_profile.add_argument(
        "--with-go", action="store_true",
        help="Enable GO hierarchy expansion (requires: pip install sceptr[go])")
    p_profile.add_argument(
        "--analysis-type", choices=["functional", "cellular"],
        default="functional",
        help="Analysis type (default: functional)")
    p_profile.add_argument(
        "--quiet", "-q", action="store_true",
        help="Suppress progress logging")
    p_profile.set_defaults(func=cmd_profile)

    # --- categories ---
    p_cats = subparsers.add_parser(
        "categories",
        help="List or inspect built-in category sets",
    )
    p_cats.add_argument("--list", "-l", action="store_true",
                        help="List all available category sets")
    p_cats.add_argument("--show",
                        help="Show categories in a specific set")
    p_cats.set_defaults(func=cmd_categories)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    args.func(args)


if __name__ == "__main__":
    main()
