#!/usr/bin/env nextflow

/*
 * SCEPTR: InterProScan Annotation Module
 *
 * Runs InterProScan (Pfam by default) on the proteome and merges the
 * resulting Pfam/InterPro/GO annotations into SCEPTR's
 * integrated_annotations_expression.tsv. Designed as the second
 * annotation source after DIAMOND/UniProt, particularly important for
 * non-model organisms where UniProt sequence-similarity coverage is
 * sparse.
 *
 * Two processes:
 *   RunInterProScan         - executes interproscan.sh on a proteome
 *   AugmentAnnotations      - merges InterProScan TSV into the existing
 *                             integrated_annotations_expression.tsv
 *
 * Author: James McCabe
 * Module: SCEPTR
 */

process RunInterProScan {
    tag "interproscan_pfam"
    label 'process_high'
    publishDir "${params.outdir}/preprocessing/interproscan", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path proteome

    output:
    path "interproscan_results.tsv", emit: iprscan_tsv

    script:
    def applications = params.interproscan_applications ?: 'Pfam'
    def interproscan_dir = params.interproscan_dir ?: '/opt/interproscan'
    """
    # Strip stop codons that TransDecoder includes (* characters);
    # InterProScan rejects them.
    sed 's/\\*//g' ${proteome} > proteome_clean.fasta

    n_in=\$(grep -c '^>' proteome_clean.fasta || echo 0)
    echo "Running InterProScan on \$n_in proteins (applications: ${applications})"

    ${interproscan_dir}/interproscan.sh \\
        -i proteome_clean.fasta \\
        -f tsv \\
        -appl ${applications} \\
        -dp \\
        -goterms \\
        -iprlookup \\
        -cpu ${task.cpus} \\
        -d .

    if [ ! -f proteome_clean.fasta.tsv ]; then
        echo "ERROR: InterProScan did not produce an output TSV"
        exit 1
    fi

    mv proteome_clean.fasta.tsv interproscan_results.tsv
    n_hits=\$(wc -l < interproscan_results.tsv)
    n_proteins=\$(cut -f1 interproscan_results.tsv | sort -u | wc -l)
    echo "InterProScan complete: \$n_hits hits across \$n_proteins proteins"
    """
}

process AugmentAnnotations {
    tag "augment_annotations"
    label 'process_low'
    publishDir "${params.outdir}/integrated_data", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path annotations
    path iprscan_tsv

    output:
    path "integrated_annotations_expression.tsv", emit: augmented_annotations

    script:
    """
    # InterProScan-only rows (proteins without UniProt annotation) get
    # NaN TPM, which sorts them to the bottom of the rank list - which
    # is the correct behaviour because they have no expression evidence.
    python3 ${params.interproscan_bin}/augment_annotations.py \\
        --annotations ${annotations} \\
        --interproscan ${iprscan_tsv} \\
        --output integrated_annotations_expression.tsv
    """
}
