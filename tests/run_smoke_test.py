#!/usr/bin/env python3
"""
Self-test for SCEPTR. Verifies that the Python package is installed and
that a complete enrichment run on a tiny bundled dataset succeeds.

Designed to be runnable on a fresh clone immediately after either
``pip install sceptr-profiling`` or ``pip install -e .``. Does not
require Docker, Nextflow, or any external databases.

Exit codes:
    0  smoke test passed
    1  smoke test failed (any step)

The smoke test is deliberately small:
    - 300 synthetic genes
    - the built-in "general" category set (11 categories)
    - 100 permutations
    - typical wall-clock under 30 seconds on a laptop

A FAIL message points users to the issue tracker so problems on real
user machines surface as actionable bug reports rather than silent
breakage.
"""

import logging
import os
import sys
import time
import traceback
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_TSV = REPO_ROOT / "tests" / "data" / "smoke_expression.tsv"
ISSUE_URL = "https://github.com/jsmccabe1/SCEPTR/issues/new"


def _print_box(title: str, body: str, char: str = "=") -> None:
    bar = char * 70
    print(bar)
    print(f"  {title}")
    print(bar)
    if body:
        print(body)
    print(bar)


def _fail(step: str, detail: str) -> int:
    body = (
        f"Step that failed: {step}\n"
        f"Detail:\n{detail}\n\n"
        f"Please report this so it can be fixed for other users:\n"
        f"  {ISSUE_URL}\n\n"
        f"Including the output above is enough; no personal data is logged."
    )
    _print_box("SCEPTR smoke test: FAIL", body, "!")
    return 1


def _check_input() -> int:
    if not INPUT_TSV.exists():
        return _fail(
            "locate bundled test data",
            f"Expected to find:\n  {INPUT_TSV}\nThis file ships with the "
            f"repository. Try a fresh clone or run:\n"
            f"  python3 tests/data/generate_smoke_dataset.py",
        )
    return 0


def _check_imports() -> int:
    try:
        import sceptr  # noqa: F401
        from sceptr import profile, continuous, enrichment, categorisation  # noqa: F401
    except Exception as e:
        return _fail(
            "import sceptr Python package",
            f"{type(e).__name__}: {e}\n\n"
            f"Try:  pip install sceptr-profiling   (or, from a clone)  pip install .",
        )
    return 0


def _run_profile(outdir: Path) -> int:
    try:
        from sceptr.profile import run as profile_run
    except Exception as e:
        return _fail("import sceptr.profile.run", f"{type(e).__name__}: {e}")

    try:
        profile_run(
            expression_file=str(INPUT_TSV),
            category_set="general",
            custom_categories=None,
            category_mapping_file=None,
            output_dir=str(outdir),
            prefix="smoke",
            tiers=[10, 25, 50, 100],
            permutations=100,
            with_go=False,
            analysis_type="functional",
            quiet=True,
        )
    except Exception:
        return _fail(
            "run sceptr.profile",
            traceback.format_exc(),
        )
    return 0


def _check_outputs(outdir: Path) -> int:
    produced = list(outdir.rglob("smoke*"))
    if not produced:
        return _fail(
            "verify output files",
            f"No files matching 'smoke*' produced anywhere under {outdir}.\n"
            f"This usually means profile() exited silently without writing.",
        )
    return 0


def main() -> int:
    import tempfile

    # Quiet down library logging so the smoke test output is clean.
    logging.getLogger("sceptr").setLevel(logging.WARNING)
    logging.getLogger("sceptr.explot").setLevel(logging.WARNING)
    logging.getLogger("sceptr.explot.reporting.interactive").setLevel(logging.WARNING)
    logging.getLogger("sceptr.profile").setLevel(logging.WARNING)

    start = time.time()
    print("SCEPTR smoke test")
    print(f"  input:    {INPUT_TSV.relative_to(REPO_ROOT)}")
    print(f"  category: general (11 built-in categories)")
    print(f"  permutations: 100")
    print()

    rc = _check_input()
    if rc:
        return rc
    print("[ ok ] bundled test data present")

    rc = _check_imports()
    if rc:
        return rc
    print("[ ok ] sceptr Python package imports")

    with tempfile.TemporaryDirectory(prefix="sceptr_smoke_") as tmp:
        outdir = Path(tmp)
        rc = _run_profile(outdir)
        if rc:
            return rc
        print(f"[ ok ] profile run completed in {time.time() - start:.1f}s")

        rc = _check_outputs(outdir)
        if rc:
            return rc
        n_out = len(list(outdir.rglob("smoke*")))
        print(f"[ ok ] {n_out} output files produced")

    elapsed = time.time() - start
    body = (
        f"SCEPTR is installed and producing valid output on this machine.\n"
        f"Total time: {elapsed:.1f} seconds.\n\n"
        f"You can now try a real dataset, e.g.:\n"
        f"  sceptr profile --expression your_data.tsv --category-set general -o results/\n\n"
        f"If anything during a real run is unclear or unexpected, a one-line\n"
        f"issue at {ISSUE_URL} is the only way the developer hears about it."
    )
    _print_box("SCEPTR smoke test: PASS", body)
    return 0


if __name__ == "__main__":
    sys.exit(main())
