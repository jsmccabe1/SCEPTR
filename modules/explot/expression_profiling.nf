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
        path "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_continuous_enrichment.tsv", emit: bp_mf_continuous, optional: true
        path "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_profile_test.tsv", emit: bp_mf_profile_test, optional: true
        path "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_continuous_dkl.tsv", emit: bp_mf_continuous_dkl, optional: true
        path "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_profile_shapes.tsv", emit: bp_mf_shapes, optional: true
        path "functional/${params.output_prefix ?: 'sceptr'}_BP_MF_report_data.json", emit: bp_mf_report_data, optional: true

    script:
        def cat_set = params.category_set ?: "general"
        def cat_path = "${params.modules_dir}/explot/categories/functional/${cat_set}_functional_categories.json"
        def custom_categories = "--custom_categories ${cat_path}"
        def tiers = params.explot_tiers ? "--tiers ${params.explot_tiers}" : ""
        def multi_category = params.explot_multi_category ? "--multi_category" : ""
        def expand_go = params.explot_expand_go ? "--expand_go --go_expansion_depth ${params.explot_go_expansion_depth}" : ""
        def continuous = params.explot_continuous ? "" : "--no_continuous"
        def cont_step = params.explot_continuous_step ? "--continuous_step ${params.explot_continuous_step}" : ""
        def cont_k_min = params.explot_continuous_k_min ? "--continuous_k_min ${params.explot_continuous_k_min}" : ""
        def cont_k_max = params.explot_continuous_k_max ? "--continuous_k_max ${params.explot_continuous_k_max}" : ""
        def cont_perms = params.explot_profile_permutations ? "--profile_permutations ${params.explot_profile_permutations}" : ""
        def debug = params.debug ? "--debug" : ""

        """
        if [ ! -f "${params.explot_cli}/functional_profiling_cli.py" ]; then
            echo "ERROR: Script not found at ${params.explot_cli}/functional_profiling_cli.py"
            exit 1
        fi

        python3 ${params.explot_cli}/functional_profiling_cli.py ${integrated_results} \\
            --prefix ${params.output_prefix ?: "sceptr"} \\
            --outdir . \\
            ${custom_categories} ${tiers} ${multi_category} ${expand_go} \\
            ${continuous} ${cont_step} ${cont_k_min} ${cont_k_max} ${cont_perms} \\
            ${debug}

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
        path "cellular/${params.output_prefix ?: 'sceptr'}_CC_continuous_enrichment.tsv", emit: cc_continuous, optional: true
        path "cellular/${params.output_prefix ?: 'sceptr'}_CC_profile_test.tsv", emit: cc_profile_test, optional: true
        path "cellular/${params.output_prefix ?: 'sceptr'}_CC_continuous_dkl.tsv", emit: cc_continuous_dkl, optional: true
        path "cellular/${params.output_prefix ?: 'sceptr'}_CC_profile_shapes.tsv", emit: cc_shapes, optional: true
        path "cellular/${params.output_prefix ?: 'sceptr'}_CC_report_data.json", emit: cc_report_data, optional: true

    script:
        def cat_set = params.category_set ?: "general"
        def cat_path = "${params.modules_dir}/explot/categories/cellular/${cat_set}_cellular_categories.json"
        def custom_categories = "--custom_categories ${cat_path}"
        def tiers = params.explot_tiers ? "--tiers ${params.explot_tiers}" : ""
        def multi_category = params.explot_multi_category ? "--multi_category" : ""
        def expand_go = params.explot_expand_go ? "--expand_go --go_expansion_depth ${params.explot_go_expansion_depth}" : ""
        def continuous = params.explot_continuous ? "" : "--no_continuous"
        def cont_step = params.explot_continuous_step ? "--continuous_step ${params.explot_continuous_step}" : ""
        def cont_k_min = params.explot_continuous_k_min ? "--continuous_k_min ${params.explot_continuous_k_min}" : ""
        def cont_k_max = params.explot_continuous_k_max ? "--continuous_k_max ${params.explot_continuous_k_max}" : ""
        def cont_perms = params.explot_profile_permutations ? "--profile_permutations ${params.explot_profile_permutations}" : ""
        def debug = params.debug ? "--debug" : ""

        """
        if [ ! -f "${params.explot_cli}/cellular_profiling_cli.py" ]; then
            echo "ERROR: Script not found at ${params.explot_cli}/cellular_profiling_cli.py"
            exit 1
        fi

        python3 ${params.explot_cli}/cellular_profiling_cli.py ${integrated_results} \\
            --prefix ${params.output_prefix ?: "sceptr"} \\
            --outdir . \\
            ${custom_categories} ${tiers} ${multi_category} ${expand_go} \\
            ${continuous} ${cont_step} ${cont_k_min} ${cont_k_max} ${cont_perms} \\
            ${debug}

        if [ ! -f "cellular/${params.output_prefix ?: 'sceptr'}_CC_report.html" ]; then
            echo "ERROR: Failed to generate report file"
            exit 1
        fi
        """
}

process GenerateCombinedReport {
    tag "combined_report"
    publishDir "${params.outdir}/enrichment_profiles", mode: 'copy'
    errorStrategy 'ignore'

    input:
        path integrated_results
        path func_report_data
        path cell_report_data

    output:
        path "${params.output_prefix ?: 'sceptr'}_report.html", emit: combined_report

    script:
        """
        python3 ${params.explot_cli}/generate_report_cli.py ${integrated_results} \\
            --functional ${func_report_data} \\
            --cellular ${cell_report_data} \\
            --output ${params.output_prefix ?: "sceptr"}_report.html
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

        // Generate unified tabbed report if both analyses produce report data
        GenerateCombinedReport(
            integrated_results,
            ExPlotFunctional.out.bp_mf_report_data,
            ExPlotCellular.out.cc_report_data
        )

    emit:
        bp_mf_plots = ExPlotFunctional.out.bp_mf_plots
        bp_mf_report = ExPlotFunctional.out.bp_mf_report
        cc_plots = ExPlotCellular.out.cc_plots
        cc_report = ExPlotCellular.out.cc_report
        combined_report = GenerateCombinedReport.out.combined_report
}
