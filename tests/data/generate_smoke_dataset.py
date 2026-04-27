#!/usr/bin/env python3
"""
Generate the bundled smoke-test expression file used by ``--smoke-test``.

Produces a deterministic 300-gene synthetic transcriptome with realistic
TPM distribution (Pareto-like apex, long low-expression tail) and GO
annotations spanning the eleven "general" SCEPTR categories. The result
is checked into the repo so that running the smoke test does not require
re-generating the dataset.

Run this script only when intentionally regenerating the bundled file:
    python3 tests/data/generate_smoke_dataset.py
"""

import os
import random

OUT_PATH = os.path.join(os.path.dirname(__file__), "smoke_expression.tsv")
N_GENES = 300
SEED = 42

CATEGORY_GO = {
    "Central Metabolism":              ["GO:0006096", "GO:0006099", "GO:0006629"],
    "Translation & Protein Synthesis": ["GO:0006412", "GO:0042254", "GO:0003735"],
    "Protein Folding & Quality Control": ["GO:0006457", "GO:0006508", "GO:0051082"],
    "DNA Replication & Repair":        ["GO:0006260", "GO:0006281", "GO:0006259"],
    "Transcriptional Regulation":      ["GO:0006355", "GO:0016070"],
    "Cell Division & Proliferation":   ["GO:0007049", "GO:0051301", "GO:0008283"],
    "Signal Transduction":             ["GO:0007165", "GO:0006468", "GO:0007166"],
    "Transport & Membrane":            ["GO:0055085", "GO:0006811", "GO:0006810"],
    "RNA Processing":                  ["GO:0006396", "GO:0016070"],
    "Cytoskeletal Organization":       ["GO:0030048", "GO:0048870", "GO:0016043"],
}

CATEGORY_PROTEIN = {
    "Central Metabolism":              "Glyceraldehyde-3-phosphate dehydrogenase",
    "Translation & Protein Synthesis": "60S ribosomal protein",
    "Protein Folding & Quality Control": "Heat shock protein 70",
    "DNA Replication & Repair":        "DNA polymerase delta",
    "Transcriptional Regulation":      "Transcription factor IIA",
    "Cell Division & Proliferation":   "Cyclin-dependent kinase",
    "Signal Transduction":             "Serine/threonine protein kinase",
    "Transport & Membrane":            "ABC transporter",
    "RNA Processing":                  "mRNA splicing factor",
    "Cytoskeletal Organization":       "Actin",
}


def main() -> None:
    rng = random.Random(SEED)
    categories = list(CATEGORY_GO.keys())

    # Plant Translation genes near the apex so the smoke test exercises a
    # genuinely apex-concentrated category. Distribute the rest with mild
    # bias so other categories are non-trivial too.
    apex_translation_count = 18
    rows = []

    for rank in range(N_GENES):
        # Pareto-like TPM: top genes have ~10000 TPM, tail has near-zero.
        tpm = max(0.01, 12000.0 * (1.0 / (rank + 1)) ** 0.95 + rng.uniform(-2.0, 2.0))

        if rank < apex_translation_count:
            cat = "Translation & Protein Synthesis"
        else:
            cat = rng.choice(categories)

        go_ids = CATEGORY_GO[cat]
        # Format matches what SCEPTR's GO parser expects in the
        # GO_Biological_Process column: 'name [GO:NNNNNNN]; name [...]'.
        go_field = "; ".join(
            f"smoke term [{gid}]" for gid in go_ids
        )

        rows.append({
            "sequence_id": f"smoke_gene_{rank+1:04d}",
            "TPM": f"{tpm:.4f}",
            "uniprot_id": f"P{(rank+1):05d}",
            "protein_name": CATEGORY_PROTEIN[cat],
            "gene_names": f"smoke{rank+1}",
            "organism": "smoke test (synthetic)",
            "length": str(rng.randint(100, 800)),
            "GO_Biological_Process": go_field,
            "GO_Cellular_Component": "",
            "GO_Molecular_Function": "",
        })

    cols = list(rows[0].keys())
    with open(OUT_PATH, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(r[c] for c in cols) + "\n")

    print(f"Wrote {len(rows)} genes to {OUT_PATH}")


if __name__ == "__main__":
    main()
