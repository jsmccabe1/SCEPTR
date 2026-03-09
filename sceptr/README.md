# SCEPTR

**Statistical Characterisation of Expression Profiles in Transcriptomes**

A statistical framework for continuous enrichment profiling of ranked transcriptomes. SCEPTR computes enrichment functions E_C(k) at every gene rank, permutation-based significance testing, and D_KL functional specialisation gradients - all from a single sample with no replicates required.

## Installation

```bash
pip install sceptr
```

For GO hierarchy support:
```bash
pip install sceptr[go]
```

## Quick Start

```bash
# List available category sets
sceptr categories --list

# Run enrichment profiling
sceptr profile --expression my_data.tsv --category-set bacteria -o results/
```

### Input formats

SCEPTR accepts several input formats:

**1. Annotated expression table** (recommended) - a TSV with gene IDs, expression values, and protein descriptions:
```
sequence_id    TPM    protein_name    GO_Biological_Process
gene001        2500   Ribosomal protein L3    translation [GO:0006412]
gene002        1800   ATP synthase subunit    ATP synthesis [GO:0015986]
```

**2. Expression + external category mapping** - when you have your own category assignments:
```bash
sceptr profile --expression expr.tsv --categories mapping.tsv -o results/
```
Where `mapping.tsv` is:
```
gene_id    category
gene001    Translation & Ribosome
gene001    Protein Folding
gene002    Central Metabolism
```

**3. Pre-categorised expression table** - a TSV with a `categories` column:
```
sequence_id    TPM    categories
gene001        2500   Translation & Ribosome;Protein Folding
gene002        1800   Central Metabolism
```

### Custom categories

You can define your own category sets as JSON:
```bash
sceptr profile --expression data.tsv --custom-categories my_categories.json -o results/
```

## What SCEPTR computes

- **Continuous enrichment profiles** E_C(k) at every integer gene rank
- **Discrete tier enrichment** with Fisher's exact test and FDR correction
- **Permutation-based global profile test** (supremum and integral statistics)
- **D_KL functional specialisation gradient** quantifying transcriptome organisation
- **Profile shape classification** (apex-concentrated, distributed, flat)
- **Interactive HTML report** with all results in a single self-contained file

## Python API

```python
from sceptr.profile import run

results = run(
    expression_file="data.tsv",
    category_set="bacteria",
    output_dir="results/",
    permutations=1000,
)

# Access results programmatically
tier_results = results["tier_results"]
continuous = results["continuous"]
k_values = continuous["k_values"]
enrichment_matrix = continuous["enrichment_matrix"]
```

## Citation

McCabe, J.S. and Janouskovek, J. (2026). SCEPTR: continuous enrichment profiling reveals functional architecture across the expression gradient.

## Full framework

For end-to-end analysis from raw reads (QC, quantification, annotation, and profiling), see the [full SCEPTR framework](https://github.com/jsmccabe1/SCEPTR) which uses Nextflow and Docker.
