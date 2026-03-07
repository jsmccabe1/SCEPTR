#!/usr/bin/env nextflow
// QC Module
// Contains processes for FastQC and MultiQC
process FastQC {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/preprocessing/qc/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip", emit: zip
    
    script:
    """
    # Debug info
    echo "Processing: ${sample_id}"
    echo "Reads file: ${reads}"
    echo "Current directory: \$(pwd)"
    echo "Directory contents:"
    ls -la
    
    # Run FastQC with default output naming
    fastqc --noextract ${reads}
    
    # Verify outputs
    echo "After FastQC - Directory contents:"
    ls -la
    """
}

process MultiQC {
    tag "multiqc_report"
    label 'process_low'
    publishDir "${params.outdir}/preprocessing/qc/multiqc", mode: 'copy'
    errorStrategy 'terminate'
    
    input:
    path('*')
    
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    
    script:
    """
    # Debug information
    echo "Running MultiQC"
    echo "Current directory contents:"
    ls -la
    
    # Run MultiQC
    multiqc --force .
    
    # Verify output exists
    if [ ! -f multiqc_report.html ] || [ ! -d multiqc_data ]; then
        echo "ERROR: MultiQC outputs not found!"
        echo "Current directory contents after MultiQC:"
        ls -la
        exit 1
    fi
    
    echo "MultiQC completed successfully."
    """
}
