#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// All parameter defaults are defined in nextflow.config (single source of truth).
// params.explot_functional_categories and params.explot_cellular_categories
// are set dynamically by resolveCategoryFiles() in main.nf at runtime.

process ExPlotFunctional {
    tag "explot_functional"
    publishDir "${params.outdir}/enrichment_profiles", mode: 'copy'
    errorStrategy 'terminate'

    input:
        path integrated_results

    output:
        path "functional/figures/*.png", emit: bp_mf_plots, optional: true
        path "functional/figures/*.svg", emit: bp_mf_plots_svg, optional: true
        path "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_report.html", emit: bp_mf_report
        path "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_enrichment_results.tsv", emit: bp_mf_enrichment_tsv, optional: true
        path "functional/${params.output_prefix ?: 'sceptr'}_expanded_categories.json", emit: expanded_categories, optional: true
        path "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_assignment_methods.json", emit: bp_mf_methods, optional: true
        path "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_core_specificity.json", emit: bp_mf_core, optional: true

    script:
        def cat_set = params.category_set ?: "general"
        def cat_path = "${params.modules_dir}/explot/categories/functional/${cat_set}_functional_categories.json"
        def custom_categories = "--custom_categories ${cat_path}"
        def tiers = params.explot_tiers ? "--tiers ${params.explot_tiers}" : ""
        def multi_category = params.explot_multi_category ? "--multi_category" : ""
        def expand_go = params.explot_expand_go ? "--expand_go --go_expansion_depth ${params.explot_go_expansion_depth}" : ""
        def debug = params.debug ? "--debug" : ""

        """
        if [ ! -f "${params.explot_cli}/functional_profiling_cli.py" ]; then
            echo "ERROR: Script not found at ${params.explot_cli}/functional_profiling_cli.py"
            exit 1
        fi

        python3 ${params.explot_cli}/functional_profiling_cli.py ${integrated_results} \\
            --prefix ${params.output_prefix ?: "sceptr"} \\
            --outdir . \\
            ${custom_categories} ${tiers} ${multi_category} ${expand_go} ${debug}

        if [ ! -f "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_report.html" ]; then
            echo "ERROR: Failed to generate report file"
            exit 1
        fi
        """
}

process ExPlotCellular {
    tag "explot_cellular"
    publishDir "${params.outdir}/enrichment_profiles", mode: 'copy'
    errorStrategy 'terminate'

    input:
        path integrated_results

    output:
        path "cellular/figures/*.png", emit: cc_plots, optional: true
        path "cellular/figures/*.svg", emit: cc_plots_svg, optional: true
        path "cellular/${params.output_prefix ?: 'sceptr'}_CC_report.html", emit: cc_report
        path "cellular/${params.output_prefix ?: 'sceptr'}_CC_enrichment_results.tsv", emit: cc_enrichment_tsv, optional: true
        path "cellular/${params.output_prefix ?: 'sceptr'}_expanded_cc_categories.json", emit: expanded_categories, optional: true
        path "cellular/${params.output_prefix ?: 'sceptr'}_CC_assignment_methods.json", emit: cc_methods, optional: true
        path "cellular/${params.output_prefix ?: 'sceptr'}_CC_core_specificity.json", emit: cc_core, optional: true

    script:
        def cat_set = params.category_set ?: "general"
        def cat_path = "${params.modules_dir}/explot/categories/cellular/${cat_set}_cellular_categories.json"
        def custom_categories = "--custom_categories ${cat_path}"
        def tiers = params.explot_tiers ? "--tiers ${params.explot_tiers}" : ""
        def multi_category = params.explot_multi_category ? "--multi_category" : ""
        def expand_go = params.explot_expand_go ? "--expand_go --go_expansion_depth ${params.explot_go_expansion_depth}" : ""
        def debug = params.debug ? "--debug" : ""

        """
        if [ ! -f "${params.explot_cli}/cellular_profiling_cli.py" ]; then
            echo "ERROR: Script not found at ${params.explot_cli}/cellular_profiling_cli.py"
            exit 1
        fi

        python3 ${params.explot_cli}/cellular_profiling_cli.py ${integrated_results} \\
            --prefix ${params.output_prefix ?: "sceptr"} \\
            --outdir . \\
            ${custom_categories} ${tiers} ${multi_category} ${expand_go} ${debug}

        if [ ! -f "cellular/${params.output_prefix ?: 'sceptr'}_CC_report.html" ]; then
            echo "ERROR: Failed to generate report file"
            exit 1
        fi
        """
}

workflow ExPlot {
    take:
        integrated_results

    main:
        log.info "Enrichment profiling workflow started"
        log.info "Using script directory: ${params.explot_cli}"
        log.info "GO expansion: ${params.explot_expand_go}"

        integrated_results
            .ifEmpty { error "No integrated results file found" }

        ExPlotFunctional(integrated_results)
        ExPlotCellular(integrated_results)

    emit:
        bp_mf_plots = ExPlotFunctional.out.bp_mf_plots
        bp_mf_report = ExPlotFunctional.out.bp_mf_report
        cc_plots = ExPlotCellular.out.cc_plots
        cc_report = ExPlotCellular.out.cc_report
}
