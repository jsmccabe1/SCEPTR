#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ============================================================================
// OPTIONAL HOST SEQUENCE FILTERING (for parasite/pathogen transcriptomes)
// ============================================================================

process TranslateHostTranscriptome {
    label 'process_medium'
    
    publishDir "${params.outdir}/preprocessing/contamination/host_filter", mode: 'copy'
    
    input:
    path host_transcriptome
    
    output:
    path "host_proteins.fasta", emit: host_proteome
    path "host_translation.log", emit: log
    
    when:
    params.host_transcriptome && !params.skip_host_filter
    
    script:
    """
    echo "Translating host transcriptome to proteins..."
    echo "Input: ${host_transcriptome}" > host_translation.log
    echo "Minimum protein length: ${params.min_prot_len}" >> host_translation.log
    echo "Using direct CDS translation (standard genetic code)" >> host_translation.log
    
    python3 << 'PYEOF'
from Bio import SeqIO
from Bio.Seq import Seq
import sys

min_len = ${params.min_prot_len}
count = 0
skipped = 0

with open("host_proteins.fasta", "w") as out_fh:
    for record in SeqIO.parse("${host_transcriptome}", "fasta"):
        seq = str(record.seq).upper().replace('U', 'T')
        # Trim to multiple of 3
        trim_len = len(seq) - (len(seq) % 3)
        if trim_len < min_len * 3:
            skipped += 1
            continue
        try:
            protein = str(Seq(seq[:trim_len]).translate(table=1))
            # Remove trailing stop codon if present
            if protein.endswith('*'):
                protein = protein[:-1]
            if len(protein) >= min_len:
                out_fh.write(f">{record.id}\\n{protein}\\n")
                count += 1
            else:
                skipped += 1
        except Exception as e:
            skipped += 1

print(f"Translated {count} sequences, skipped {skipped}")
PYEOF

    protein_count=\$(grep -c "^>" host_proteins.fasta || echo "0")
    echo "Translated sequences: \$protein_count" >> host_translation.log
    
    if [ "\$protein_count" -eq 0 ]; then
        echo "ERROR: No proteins were predicted from host transcriptome"
        exit 1
    fi
    
    echo "Host translation complete" >> host_translation.log
    echo "Host transcriptome successfully translated to \$protein_count proteins"
    
    cat host_translation.log
    """
}

process BuildHostDatabase {
    label 'process_medium'
    
    publishDir "${params.outdir}/preprocessing/contamination/host_filter", mode: 'copy'
    
    input:
    path host_proteome
    
    output:
    path "host.dmnd", emit: host_db
    path "host_db_build.log", emit: log
    
    when:
    (params.host_transcriptome || params.host_proteome) && !params.skip_host_filter
    
    script:
    """
    echo "Building DIAMOND database from host proteome..."
    echo "Host proteome: ${host_proteome}"
    
    # Verify input is valid FASTA
    if ! grep -q "^>" "${host_proteome}"; then
        echo "ERROR: Host proteome does not appear to be in FASTA format"
        exit 1
    fi
    
    # Count sequences in host reference
    host_seq_count=\$(grep -c "^>" "${host_proteome}")
    echo "Host protein sequences: \$host_seq_count"
    
    if [ "\$host_seq_count" -eq 0 ]; then
        echo "ERROR: No sequences found in host proteome"
        exit 1
    fi
    
    # Build DIAMOND database from proteins (no --ignore-warnings needed)
    diamond makedb \\
        --in "${host_proteome}" \\
        --db host \\
        --threads ${task.cpus} \\
        --verbose
    
    # Log build info
    echo "Host DIAMOND database built successfully" > host_db_build.log
    echo "Input file: ${host_proteome}" >> host_db_build.log
    echo "Sequences indexed: \$host_seq_count" >> host_db_build.log
    echo "Database file: host.dmnd" >> host_db_build.log
    echo "Build date: \$(date)" >> host_db_build.log
    
    # Verify database was created
    if [ ! -f "host.dmnd" ]; then
        echo "ERROR: Failed to create DIAMOND database"
        exit 1
    fi
    
    echo "Host database build complete"
    """
}

process FilterHostSequences {
    label 'process_medium'
    
    publishDir "${params.outdir}/preprocessing/contamination/host_filter", mode: 'copy'
    
    input:
    path proteome
    path host_db
    
    output:
    path "after_host_filter.fasta", emit: filtered_proteome
    path "host_sequences.fasta", emit: host_sequences, optional: true
    path "host_filter_report.txt", emit: report
    path "host_blast_results.tsv", emit: blast_results, optional: true
    
    when:
    (params.host_transcriptome || params.host_proteome) && !params.skip_host_filter
    
    script:
    """
    echo "SCEPTR Host Sequence Filtering"
    echo "=============================="
    echo "Input proteome: ${proteome}"
    echo "Host database: ${host_db}"
    echo ""
    echo "Parameters:"
    echo "  Identity threshold: ${params.host_identity}%"
    echo "  Coverage threshold: ${params.host_coverage}%"
    echo "  E-value threshold: ${params.host_evalue}"
    echo ""
    
    # Count input sequences
    input_count=\$(grep -c "^>" "${proteome}")
    echo "Input sequences: \$input_count"
    
    # Run DIAMOND BLAST against host database
    echo "Running DIAMOND BLAST against host sequences..."
    diamond blastp \\
        --query "${proteome}" \\
        --db "${host_db}" \\
        --out host_blast_results.tsv \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \\
        --id ${params.host_identity} \\
        --query-cover ${params.host_coverage} \\
        --evalue ${params.host_evalue} \\
        --max-target-seqs 5 \\
        --threads ${task.cpus} \\
        --sensitive
    
    # Check if BLAST produced results
    if [ ! -f "host_blast_results.tsv" ]; then
        echo "No BLAST results - creating empty file"
        touch host_blast_results.tsv
    fi
    
    # Create Python script for host sequence filtering
    cat << 'PYEOF' > filter_host_sequences.py
#!/usr/bin/env python3

import sys
from Bio import SeqIO
from collections import defaultdict

print("Processing host BLAST results...")

# Parse BLAST results
host_sequences = set()
host_details = []

try:
    with open('host_blast_results.tsv', 'r') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\\t')
            if len(parts) >= 13:
                qseqid = parts[0]
                sseqid = parts[1]
                pident = float(parts[2])
                evalue = float(parts[10])
                stitle = parts[12] if len(parts) > 12 else "Unknown"
                
                host_sequences.add(qseqid)
                host_details.append({
                    'query': qseqid,
                    'subject': sseqid,
                    'identity': pident,
                    'evalue': evalue,
                    'description': stitle
                })
except FileNotFoundError:
    print("No BLAST results file found")

print(f"Identified {len(host_sequences)} sequences matching host")

# Filter proteome
input_count = 0
output_count = 0
host_count = 0

with open('${proteome}', 'r') as infile, \\
     open('after_host_filter.fasta', 'w') as outfile, \\
     open('host_sequences.fasta', 'w') as hostfile:
    
    for record in SeqIO.parse(infile, 'fasta'):
        input_count += 1
        if record.id in host_sequences:
            SeqIO.write(record, hostfile, 'fasta')
            host_count += 1
        else:
            SeqIO.write(record, outfile, 'fasta')
            output_count += 1

# Generate report
print(f"\\nFiltering complete:")
print(f"  Input sequences: {input_count}")
print(f"  Host sequences: {host_count}")
print(f"  Retained sequences: {output_count}")

with open('host_filter_report.txt', 'w') as report:
    report.write("SCEPTR Host Sequence Removal Report\\n")
    report.write("===================================\\n\\n")
    report.write(f"Host Database: ${host_db}\\n\\n")
    report.write("Filtering Parameters:\\n")
    report.write(f"- Identity threshold: ${params.host_identity}%\\n")
    report.write(f"- Coverage threshold: ${params.host_coverage}%\\n")
    report.write(f"- E-value threshold: ${params.host_evalue}\\n\\n")
    report.write("Results Summary:\\n")
    report.write(f"- Input sequences: {input_count}\\n")
    host_pct = (host_count / input_count * 100) if input_count > 0 else 0.0
    report.write(f"- Host sequences identified: {host_count} ({host_pct:.1f}%)\\n")
    report.write(f"- Sequences after host filtering: {output_count}\\n\\n")
    
    if host_details:
        report.write("Top Host Matches (by identity):\\n")
        sorted_hits = sorted(host_details, key=lambda x: x['identity'], reverse=True)
        for i, hit in enumerate(sorted_hits[:10], 1):
            report.write(f"{i}. {hit['query']} → {hit['description'][:60]} ")
            report.write(f"({hit['identity']:.1f}% identity, E-value: {hit['evalue']:.2e})\\n")
    else:
        report.write("No host sequences detected above threshold.\\n")
    
    report.write("\\nHost sequences saved to: host_sequences.fasta\\n")
    report.write("Filtered proteome saved to: after_host_filter.fasta\\n")
    report.write("\\nProceeding to general contamination filtering...\\n")

PYEOF

    chmod +x filter_host_sequences.py
    python3 filter_host_sequences.py
    
    # Display summary
    cat host_filter_report.txt
    
    # Verify output
    final_count=\$(grep -c "^>" after_host_filter.fasta)
    echo ""
    echo "Host filtering complete: \$final_count sequences remaining"
    """
}

// ============================================================================
// GENERAL CONTAMINATION DETECTION (always runs, after optional host filtering)
// ============================================================================

process DiamondBlastContaminants {
    label 'process_medium'
    
    publishDir "${params.outdir}/preprocessing/contamination", mode: 'copy'
    
    input:
    path proteome
    path contaminants_db
    
    output:
    path "filtered_proteome.fasta", emit: filtered_proteome
    path "contaminant_details.csv", emit: details
    path "contaminant_report.txt", emit: report
    path "contaminants_hits.tsv", emit: blast_results, optional: true
    path "filtering_params.json", emit: params_json, optional: true
    path "contaminant_sequences.fasta", emit: contaminant_seqs, optional: true
    path "*_contaminants.fasta", emit: category_fastas, optional: true
    
    script:
    def input_source = (params.host_transcriptome || params.host_proteome) && !params.skip_host_filter ? 
        "after host filtering" : "TransDecoder"
    """
    echo "SCEPTR General Contamination Detection"
    echo "======================================"
    echo "Input proteome (${input_source}): ${proteome}"
    echo "Contaminants database: ${contaminants_db}"
    
    # Verify input files exist
    if [ ! -f "${proteome}" ]; then
        echo "ERROR: Proteome file not found: ${proteome}"
        exit 1
    fi
    
    if [ ! -f "${contaminants_db}" ]; then
        echo "ERROR: Contaminants database not found: ${contaminants_db}"
        exit 1
    fi
    
    # Count input sequences
    input_count=\$(grep -c "^>" "${proteome}" || echo "0")
    echo "Input sequences: \$input_count"
    
    # Run DIAMOND BLAST against contaminants database
    echo "Running DIAMOND BLAST against general contaminants database..."
    diamond blastp \\
        --query "${proteome}" \\
        --db "${contaminants_db}" \\
        --out contaminants_hits.tsv \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen stitle \\
        --evalue ${params.evalue_threshold} \\
        --max-target-seqs ${params.max_target_seqs} \\
        --query-cover ${params.coverage_threshold} \\
        --threads ${task.cpus} \\
        --sensitive
    
    # Check if BLAST produced results
    if [ ! -f "contaminants_hits.tsv" ]; then
        echo "No BLAST results file created"
        touch contaminants_hits.tsv
    fi
    
    # Run contamination filtering using external Python script
    echo "Running contamination filtering..."
    python3 ${params.contamination_bin}/filter_contaminants.py \\
        --identity ${params.identity_threshold} \\
        --coverage ${params.coverage_threshold} \\
        --evalue ${params.evalue_threshold} \\
        --proteome "${proteome}" \\
        --blast contaminants_hits.tsv
    
    # Verify outputs were created
    if [ ! -f "filtered_proteome.fasta" ]; then
        echo "ERROR: Output file filtered_proteome.fasta was not created"
        exit 1
    fi
    
    # Final verification
    output_count=\$(grep -c "^>" filtered_proteome.fasta || echo "0")
    echo ""
    echo "General contamination filtering complete: \$output_count clean sequences"
    """
}

process VisualiseContaminants {
    label 'process_low'
    
    publishDir "${params.outdir}/preprocessing/contamination/visualisation", mode: 'copy'
    
    input:
    path contaminant_details
    path params_json
    
    output:
    path "contaminant_analysis_report.html", emit: html_report, optional: true
    path "*.png", emit: plots, optional: true
    
    script:
    """
    echo "Starting contamination visualization..."
    
    # Check if we have contamination data
    if [ ! -f "${contaminant_details}" ] || [ ! -s "${contaminant_details}" ]; then
        echo "No contamination data found or file is empty"
        
        # Create minimal HTML report for no contamination case
        cat > contaminant_analysis_report.html << 'HTMLEOF'
<!DOCTYPE html>
<html>
<head>
    <title>SCEPTR Contamination Analysis</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        .success { background-color: #d4edda; padding: 20px; border-radius: 5px; border-left: 5px solid #28a745; }
    </style>
</head>
<body>
    <h1>SCEPTR Contamination Analysis Report</h1>
    <div class="success">
        <h2>✓ No Contamination Detected</h2>
        <p>No contaminating sequences were identified in this dataset.</p>
        <p>This indicates either:</p>
        <ul>
            <li>Very clean input data</li>
            <li>No matches to the contamination database</li>
            <li>Filtering parameters excluded all potential contaminants</li>
        </ul>
    </div>
</body>
</html>
HTMLEOF
        
        echo "Created minimal report for no contamination case"
        exit 0
    fi
    
    # Run the external contamination visualization script
    echo "Running contamination visualization with external script..."
    Rscript ${params.contamination_bin}/contaminant_viz.R "${contaminant_details}" "."
    
    # Verify output was created
    if [ ! -f "contaminant_analysis_report.html" ]; then
        echo "Warning: HTML report was not created by visualization script"
        
        # Create a basic fallback report
        cat > contaminant_analysis_report.html << 'HTMLEOF'
<!DOCTYPE html>
<html>
<head><title>SCEPTR Contamination Analysis</title></head>
<body>
    <h1>Contamination Analysis Report</h1>
    <p>Visualization script completed but report generation failed.</p>
    <p>Please check the contaminant_details.csv file for raw data.</p>
</body>
</html>
HTMLEOF
    fi
    
    echo "Contamination visualization completed"
    """
}
