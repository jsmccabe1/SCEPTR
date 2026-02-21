#!/usr/bin/env nextflow
// TransDecoder Module
// Contains the TransDecoder process for ORF prediction
// and DirectTranslate for pre-existing CDS inputs (bacteria, viruses)
process TransDecoder {
    tag "transdecoder"
    publishDir "${params.outdir}/proteome", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path transcripts

    output:
    path "${transcripts}.transdecoder.pep", emit: predicted_proteome
    path "${transcripts}.transdecoder.*", emit: all_outputs
    
    script:
    """
    # Debug info
    echo "Transcripts file: ${transcripts}"
    echo "Current directory: \$(pwd)"
    echo "Directory contents:"
    ls -la
    
    # Verify input file
    if [ ! -f "${transcripts}" ]; then
        echo "ERROR: Input file ${transcripts} not found"
        exit 1
    fi
    
    echo "File size: \$(du -h ${transcripts} | cut -f1)"
    echo "Line count: \$(wc -l ${transcripts})"
    echo "First few lines:"
    head -n 5 "${transcripts}"
    
    # Run TransDecoder
    echo "Running TransDecoder.LongOrfs..."
    TransDecoder.LongOrfs -m ${params.min_prot_len ?: 50} -t "${transcripts}"
    
    if [ \$? -ne 0 ]; then
        echo "ERROR: TransDecoder.LongOrfs failed"
        exit 1
    fi
    
    echo "Running TransDecoder.Predict..."
    TransDecoder.Predict -t "${transcripts}"
    
    if [ \$? -ne 0 ]; then
        echo "ERROR: TransDecoder.Predict failed"
        exit 1
    fi
    
    # Verify output
    if [ ! -f "${transcripts}.transdecoder.pep" ]; then
        echo "ERROR: Expected output file not found: ${transcripts}.transdecoder.pep"
        echo "Available files:"
        ls -la
        find . -name "*.pep" -type f
        exit 1
    fi
    
    echo "TransDecoder completed successfully"
    echo "Output size: \$(du -h ${transcripts}.transdecoder.pep | cut -f1)"
    """
}

// DirectTranslate: translates CDS nucleotide FASTA directly to protein
// Used when --skip_transdecoder is set (e.g. bacterial/viral CDS inputs)
// This avoids TransDecoder's ORF-finding which fails on pre-existing CDS
process DirectTranslate {
    tag "direct_translate"
    publishDir "${params.outdir}/proteome", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path transcripts

    output:
    path "translated_proteome.pep", emit: predicted_proteome
    path "translation_report.txt", emit: all_outputs
    
    script:
    """
    #!/usr/bin/env python3
    import sys
    from Bio import SeqIO
    from Bio.Seq import Seq

    input_file = "${transcripts}"
    output_file = "translated_proteome.pep"
    report_file = "translation_report.txt"

    total = 0
    translated = 0
    skipped_short = 0
    skipped_no_start = 0
    has_stop_internal = 0
    min_len = ${params.min_prot_len ?: 50}

    with open(output_file, 'w') as out_fh:
        for record in SeqIO.parse(input_file, 'fasta'):
            total += 1
            nuc_seq = str(record.seq).upper()

            # Translate in frame (CDS should already be in-frame)
            # Use standard bacterial/archaeal codon table (table 11)
            try:
                protein = str(Seq(nuc_seq).translate(table=11, to_stop=False))
            except Exception as e:
                print(f"Warning: could not translate {record.id}: {e}", file=sys.stderr)
                continue

            # Remove trailing stop codon if present
            if protein.endswith('*'):
                protein = protein[:-1]

            # Check for internal stop codons
            if '*' in protein:
                has_stop_internal += 1
                # Still include but trim at first internal stop
                protein = protein.split('*')[0]

            # Check minimum length
            if len(protein) < min_len:
                skipped_short += 1
                continue

            translated += 1
            # Write in TransDecoder-compatible format
            out_fh.write(f">{record.id}.p1 type:complete len:{len(protein)} {record.description}\\n")
            # Write protein sequence in 60-char lines
            for i in range(0, len(protein), 60):
                out_fh.write(protein[i:i+60] + "\\n")

    # Write report
    with open(report_file, 'w') as rpt:
        rpt.write("=== Direct CDS Translation Report ===\\n")
        rpt.write(f"Input sequences: {total}\\n")
        rpt.write(f"Successfully translated: {translated}\\n")
        rpt.write(f"Skipped (too short, <{min_len} aa): {skipped_short}\\n")
        rpt.write(f"Sequences with internal stop codons: {has_stop_internal}\\n")
        rpt.write(f"Translation table: 11 (Bacterial/Archaeal/Plant Plastid)\\n")
        rpt.write(f"Output file: {output_file}\\n")

    print(f"Direct translation complete: {translated}/{total} sequences translated")
    print(f"Skipped {skipped_short} sequences shorter than {min_len} aa")
    if has_stop_internal > 0:
        print(f"Warning: {has_stop_internal} sequences had internal stop codons (trimmed)")
    """
}
