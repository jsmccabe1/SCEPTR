"""
SCEPTR: Statistical Characterisation of Expression Profiles in Transcriptomes.

A statistical framework for continuous enrichment profiling of ranked
transcriptomes. Computes enrichment functions E_C(k) at every gene rank,
permutation-based significance testing, and D_KL functional specialisation
gradients from single samples.

Usage:
    # Command line
    sceptr profile --expression data.tsv --category-set bacteria -o results/

    # Python API
    from sceptr import profile
    results = profile.run(expression_file="data.tsv", category_set="bacteria")
"""

__version__ = "1.0.0"
__author__ = "James S. McCabe"
