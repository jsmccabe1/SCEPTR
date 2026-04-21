#!/usr/bin/env nextflow

/*
 * Focused integration test for the new InterProScan + AugmentAnnotations
 * processes. Exercises them through Nextflow + Docker so we validate:
 *
 *   - Process module syntax under real Nextflow execution
 *   - Channel wiring between RunInterProScan and AugmentAnnotations
 *   - File staging into the container
 *   - Docker bind mount of /opt/interproscan
 *   - The AugmentAnnotations script finding /app/bin/interproscan/...
 *
 * Doesn't test main.nf's conditional `if (!params.skip_interproscan)` block;
 * that's covered by code review.
 *
 * Usage:
 *   nextflow run test_interproscan_integration.nf \
 *     --proteome    <fasta path> \
 *     --annotations <integrated_annotations_expression.tsv path> \
 *     --outdir      test_interproscan_output \
 *     -profile docker
 */

include { RunInterProScan; AugmentAnnotations } from "${projectDir}/modules/interproscan.nf"

workflow {
    if (!params.proteome || !params.annotations) {
        error """
        Missing required inputs. Provide:
          --proteome    <fasta>
          --annotations <integrated_annotations_expression.tsv>
        """
    }

    proteome_ch    = Channel.fromPath(params.proteome,    checkIfExists: true)
    annotations_ch = Channel.fromPath(params.annotations, checkIfExists: true)

    log.info "Test InterProScan integration"
    log.info "  proteome:    ${params.proteome}"
    log.info "  annotations: ${params.annotations}"

    ipr_results = RunInterProScan(proteome_ch)
    augmented   = AugmentAnnotations(
        annotations_ch,
        ipr_results.iprscan_tsv
    )

    augmented.augmented_annotations.view { f ->
        "TEST PASSED: augmented annotation produced: ${f}"
    }
}
