#!/usr/bin/env nextflow

/*
 * SCEPTR Cross-Sample Comparison Module
 *
 * Compares functional enrichment profiles between two conditions using a
 * gene-label permutation test with concordance metrics. Takes two pre-existing
 * SCEPTR integrated_annotations_expression.tsv files as input.
 *
 * Author: James McCabe
 * Module: SCEPTR
 */

process CompareConditions {
    label 'process_high'
    publishDir "${params.outdir}/comparison", mode: 'copy'

    input:
    path condition_a
    path condition_b
    val label_a
    val label_b
    val category_set
    path functional_categories
    val tiers
    val n_permutations
    val output_prefix
    val expand_go
    val go_expansion_depth
    val seed

    output:
    path "*.html",                 emit: report
    path "*.tsv",                  emit: results
    path "figures/*.png",          emit: figures_png, optional: true
    path "figures/*.svg",          emit: figures_svg, optional: true

    script:
    def expand_flag = expand_go ? "--expand_go --go_expansion_depth ${go_expansion_depth}" : ""
    def cat_flag = functional_categories.name != 'NO_CUSTOM_CATEGORIES' ? "--functional_categories ${functional_categories}" : "--category_set ${category_set}"
    """
    python3 ${projectDir}/bin/sceptr_compare.py \
        --condition_a ${condition_a} \
        --condition_b ${condition_b} \
        --label_a "${label_a}" \
        --label_b "${label_b}" \
        ${cat_flag} \
        --tiers "${tiers}" \
        --n_permutations ${n_permutations} \
        --output_dir . \
        --output_prefix ${output_prefix} \
        --seed ${seed} \
        ${expand_flag}
    """
}
