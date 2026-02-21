#!/usr/bin/env nextflow
// Salmon Module
// Contains processes for SalmonIndex and SalmonQuant
process SalmonIndex {
    tag "salmon_index"
    label 'process_medium'
    publishDir "${params.outdir}/quantification/index", mode: 'copy'
    
    input:
        path transcripts
    
    output:
        path "salmon_index", emit: index
    
    script:
    """
    salmon index -t ${transcripts} -i salmon_index --threads ${task.cpus}
    """
}

process SalmonQuant {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/quantification", mode: 'copy'
    errorStrategy 'terminate'
    
    input:
    tuple val(sample_id), path(reads)
    path index
    
    output:
    tuple val(sample_id), path("quant.sf"), emit: quant
    
    script:
    def is_single = params.single_end ?: false
    if (is_single) {
        """
        echo "Processing sample (single-end): ${sample_id}"
        echo "Read file: ${reads}"
        echo "Index directory: ${index}"

        mkdir -p logs

        salmon quant \\
            -i ${index} \\
            -l A \\
            -r ${reads} \\
            --fldMean ${params.fld_mean ?: 250} \\
            --fldSD ${params.fld_sd ?: 25} \\
            --validateMappings \\
            --threads ${task.cpus} \\
            -o . \\
            2> logs/salmon_quant.log

        if [ ! -f quant.sf ]; then
            echo "ERROR: Salmon quantification failed - no quant.sf file produced"
            exit 1
        fi
        """
    } else {
        def read1 = reads[0]
        def read2 = reads[1]
        """
        echo "Processing sample (paired-end): ${sample_id}"
        echo "Read files: ${read1}, ${read2}"
        echo "Index directory: ${index}"

        mkdir -p logs

        salmon quant \\
            -i ${index} \\
            -l A \\
            -1 ${read1} \\
            -2 ${read2} \\
            --validateMappings \\
            --threads ${task.cpus} \\
            -o . \\
            2> logs/salmon_quant.log

        if [ ! -f quant.sf ]; then
            echo "ERROR: Salmon quantification failed - no quant.sf file produced"
            exit 1
        fi
        """
    }
}
