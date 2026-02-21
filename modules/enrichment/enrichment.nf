#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process PrepareGOEnrichmentData {
    tag "prepare_go_enrichment"
    publishDir "${params.outdir}/go_enrichment/data", mode: 'copy'
    errorStrategy 'terminate'
    
    input:
    path uniprot_file
    
    output:
    path "go_terms_long_format.csv", emit: go_terms
    path "go_terms_analysis.csv", emit: go_analysis
    
    script:
    """
    cp ${params.enrichment_bin}/prepare_go_enrichment.py ./
    python3 prepare_go_enrichment.py ${uniprot_file} \
        --output go_terms_long_format.csv \
        --analysis go_terms_analysis.csv
    """
}

process MergeExpressionAnnotations {
    tag "merge_expression_annotations"
    publishDir "${params.outdir}/integrated_data", mode: 'copy'
    errorStrategy 'terminate'
    
    input:
    path annotations_file
    path quant_dir
    
    output:
    path "integrated_annotations_expression.tsv", emit: merged_annotations
    
    script:
    """
    cp ${params.enrichment_bin}/merge_annotations_expression.py ./
    python3 merge_annotations_expression.py \
        --annotations ${annotations_file} \
        --quant ${quant_dir} \
        --output integrated_annotations_expression.tsv
    """
}

process TPMGOEnrichment {
    tag "tpm_go_enrichment"
    publishDir "${params.outdir}/go_enrichment/data", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/go_enrichment/reports", mode: 'copy', pattern: '*.html'
    publishDir "${params.outdir}/go_enrichment/figures", mode: 'copy', pattern: '*.{png,svg}'
    errorStrategy 'terminate'
    time '4h'
    memory '8 GB'

    input:
    path go_data_file
    path expression_file

    output:
    path "${params.output_prefix ?: 'sceptr'}_go_enrichment_results.tsv", emit: go_enrichment_results
    path "${params.output_prefix ?: 'sceptr'}_tpm_enrichment.html", emit: go_enrichment_html
    path "${params.output_prefix ?: 'sceptr'}_*.png", optional: true, emit: go_enrichment_plots
    path "${params.output_prefix ?: 'sceptr'}_*.svg", optional: true, emit: go_enrichment_plots_svg
    
    script:
    def prefix = params.output_prefix ?: 'sceptr'
    def tpm_threshold = params.tpm_threshold ?: 1.0
    def expression_tiers = params.expression_tiers ?: '50,100,250,500'
    def pvalue_cutoff = params.go_pvalue_cutoff ?: 0.05
    def bootstrap_n = params.bootstrap_n ?: 50
    
    """
    cp ${params.enrichment_bin}/tpm_go_enrichment_analysis.R ./
    
    Rscript tpm_go_enrichment_analysis.R \\
        --go_data ${go_data_file} \\
        --expression ${expression_file} \\
        --output_prefix ${prefix} \\
        --tpm_threshold ${tpm_threshold} \\
        --tiers "${expression_tiers}" \\
        --pvalue ${pvalue_cutoff} \\
        --bootstrap_n ${bootstrap_n}
    """
}
