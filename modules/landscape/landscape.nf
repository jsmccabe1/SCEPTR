#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process TranscriptomeLandscape {
    tag "landscape"
    publishDir "${params.outdir}/landscape", mode: 'copy'
    errorStrategy 'ignore'

    input:
        path integrated_results
        path categories_json

    output:
        path "${params.output_prefix ?: 'sceptr'}_landscape_report.html", emit: landscape_report
        path "${params.output_prefix ?: 'sceptr'}_landscape_stats.json", emit: landscape_stats, optional: true
        path "figures/${params.output_prefix ?: 'sceptr'}_expression_distribution.png", emit: expr_dist_png, optional: true
        path "figures/${params.output_prefix ?: 'sceptr'}_expression_distribution.svg", emit: expr_dist_svg, optional: true
        path "figures/${params.output_prefix ?: 'sceptr'}_annotation_quality.png", emit: ann_quality_png, optional: true
        path "figures/${params.output_prefix ?: 'sceptr'}_annotation_quality.svg", emit: ann_quality_svg, optional: true
        path "figures/${params.output_prefix ?: 'sceptr'}_taxonomic_distribution.png", emit: tax_dist_png, optional: true
        path "figures/${params.output_prefix ?: 'sceptr'}_taxonomic_distribution.svg", emit: tax_dist_svg, optional: true
        path "figures/${params.output_prefix ?: 'sceptr'}_functional_shifts.png", emit: func_shifts_png, optional: true
        path "figures/${params.output_prefix ?: 'sceptr'}_functional_shifts.svg", emit: func_shifts_svg, optional: true

    script:
        def categories_arg = categories_json.name != 'NO_CATEGORIES' ? "--categories ${categories_json}" : ""
        def top_genes = params.landscape_top_genes ?: 30
        def top_taxa = params.landscape_top_taxa ?: 20

        """
        python3 ${params.landscape_script ?: '/app/modules/landscape/transcriptome_landscape.py'} \
            ${integrated_results} \
            --prefix ${params.output_prefix ?: "sceptr"} \
            --outdir . \
            --top_genes ${top_genes} \
            --top_taxa ${top_taxa} \
            ${categories_arg}
        """
}

workflow Landscape {
    take:
        integrated_results

    main:
        log.info "Transcriptome Landscape analysis"

        // Resolve categories file - use host path for Nextflow staging
        // params.explot_functional_categories points to container path; derive host path
        def cat_set = params.category_set ?: 'general'
        def host_categories = "${projectDir}/modules/explot/categories/functional/${cat_set}_functional_categories.json"

        if (cat_set != 'custom' && file(host_categories).exists()) {
            categories_ch = Channel.fromPath(host_categories)
        } else if (params.custom_functional_categories && file(params.custom_functional_categories).exists()) {
            categories_ch = Channel.fromPath(params.custom_functional_categories)
        } else {
            // Create an empty sentinel file so the process can still run
            def sentinel = file("${workDir}/NO_CATEGORIES")
            sentinel.text = '{}'
            categories_ch = Channel.fromPath(sentinel)
        }

        TranscriptomeLandscape(integrated_results, categories_ch)

    emit:
        report = TranscriptomeLandscape.out.landscape_report
}
