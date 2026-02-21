#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Annotation Module
// Contains processes for functional annotation and GO term assignment

process DiamondBlastUniProt {
    tag "diamond_uniprot"
    label 'process_medium'
    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path proteome
    path uniprot_db

    output:
    path "uniprot_results.tsv", emit: blast_results

    script:
    """
    echo "Running DIAMOND BLASTp against UniProt database..."

    if [[ "${uniprot_db}" != *.dmnd ]]; then
        echo "Creating diamond database from FASTA..."
        diamond makedb --in ${uniprot_db} --db uniprot_db.dmnd
        DIAMOND_DB="uniprot_db.dmnd"
    else
        DIAMOND_DB="${uniprot_db}"
    fi

    diamond blastp -d \$DIAMOND_DB \\
        -q ${proteome} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
        --evalue ${params.uniprot_evalue} \\
        --max-target-seqs ${params.uniprot_max_target_seqs} \\
        --threads ${task.cpus} \\
        --out uniprot_results.tsv

    echo "UniProt BLAST search completed successfully."
    """
}

process ExtractUniProtIDs {
    tag "extract_uniprot"
    publishDir "${params.outdir}/annotation", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path blast_results

    output:
    path "sequence_uniprot_mapping.tsv", emit: uniprot_mapping

    script:
    """
    # Create inline Python script for UniProt ID extraction
    cat << 'EOF' > extract_uniprot_ids.py
#!/usr/bin/env python3

import sys
import pandas as pd

def extract_uniprot_ids(blast_results_file, output_file):
    print(f"Processing BLAST results: {blast_results_file}")
    
    # Read BLAST results
    try:
        columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        
        blast_df = pd.read_csv(blast_results_file, sep='\\t', names=columns, header=None)
        print(f"Loaded {len(blast_df)} BLAST results")
        
        if len(blast_df) == 0:
            print("Warning: No BLAST results found")
            # Create empty output file
            empty_df = pd.DataFrame(columns=['sequence_id', 'uniprot_id'])
            empty_df.to_csv(output_file, sep='\\t', index=False)
            return
        
        # Extract UniProt IDs from subject sequence IDs
        # Handle various UniProt ID formats: sp|P12345|NAME, tr|A0A123|NAME, or just P12345
        def extract_uniprot_id(sseqid):
            if '|' in sseqid:
                parts = sseqid.split('|')
                if len(parts) >= 2 and parts[0] in ['sp', 'tr']:
                    return parts[1]  # UniProt ID is the second part
                elif len(parts) >= 1:
                    return parts[0]  # Take first part
            return sseqid  # Return as-is if no | separator
        
        blast_df['uniprot_id'] = blast_df['sseqid'].apply(extract_uniprot_id)
        
        # Create mapping dataframe
        mapping_df = blast_df[['qseqid', 'uniprot_id']].copy()
        mapping_df.columns = ['sequence_id', 'uniprot_id']
        
        # Remove duplicates, keeping first occurrence
        mapping_df = mapping_df.drop_duplicates(subset=['sequence_id'])
        
        # Save mapping
        mapping_df.to_csv(output_file, sep='\\t', index=False)
        
        print(f"Extracted UniProt IDs for {len(mapping_df)} sequences")
        print(f"Output saved to {output_file}")
        
    except Exception as e:
        print(f"Error processing BLAST results: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 extract_uniprot_ids.py <blast_results.tsv> <output_mapping.tsv>")
        sys.exit(1)
    
    blast_results_file = sys.argv[1]
    output_file = sys.argv[2]
    
    extract_uniprot_ids(blast_results_file, output_file)
EOF

    # Make script executable and run it
    chmod +x extract_uniprot_ids.py
    python3 extract_uniprot_ids.py ${blast_results} sequence_uniprot_mapping.tsv
    """
}

process AnnotateProteins {
    tag "protein_annotation"
    label 'process_low'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    errorStrategy 'terminate'
    
    input:
    path uniprot_mapping
    
    output:
    path "protein_annotations.tsv", emit: protein_annotations
    
    script:
    """
    # Create inline Python script for protein annotation using batch UniProt API
    cat << 'PYEOF' > annotate_proteins.py
#!/usr/bin/env python3

import sys
import pandas as pd
import requests
import time

def batch_fetch_uniprot(uniprot_ids, batch_size=100, max_retries=3):
    \"\"\"Fetch protein info for multiple UniProt IDs using batch search API.
    
    Uses the UniProt REST API search endpoint with OR-joined accession queries.
    Much faster than individual lookups: ~2s per batch of 400 vs ~1.5s per single ID.
    \"\"\"
    
    all_results = {}
    total = len(uniprot_ids)
    
    for start in range(0, total, batch_size):
        batch = uniprot_ids[start:start + batch_size]
        batch_num = start // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size
        print(f"Fetching batch {batch_num}/{total_batches} ({len(batch)} IDs)...")
        
        # Build OR query: accession:P12345 OR accession:Q67890 ...
        query = " OR ".join([f"accession:{uid}" for uid in batch])
        
        params = {
            'query': query,
            'fields': 'accession,id,protein_name,gene_names,organism_name,length,go_p,go_c,go_f',
            'format': 'tsv',
            'size': str(len(batch))
        }
        
        for attempt in range(max_retries):
            try:
                response = requests.get(
                    "https://rest.uniprot.org/uniprotkb/search",
                    params=params,
                    timeout=120
                )
                
                if response.status_code == 200:
                    lines = response.text.strip().split('\\n')
                    if len(lines) >= 2:
                        header = lines[0].split('\\t')
                        for line in lines[1:]:
                            fields = line.split('\\t')
                            if len(fields) >= 6:
                                accession = fields[0]
                                all_results[accession] = {
                                    'accession': accession,
                                    'entry_name': fields[1] if len(fields) > 1 else '',
                                    'protein_name': fields[2] if len(fields) > 2 and fields[2] else 'Unknown protein',
                                    'gene_names': fields[3] if len(fields) > 3 else '',
                                    'organism': fields[4] if len(fields) > 4 else 'Unknown organism',
                                    'length': fields[5] if len(fields) > 5 else '0',
                                    'go_biological_process': fields[6] if len(fields) > 6 and fields[6] else '',
                                    'go_cellular_component': fields[7] if len(fields) > 7 and fields[7] else '',
                                    'go_molecular_function': fields[8] if len(fields) > 8 and fields[8] else ''
                                }
                    print(f"  Retrieved {len(all_results)} total annotations so far")
                    break  # Success, move to next batch
                    
                elif response.status_code == 429:
                    # Rate limited - wait and retry
                    retry_after = int(response.headers.get('Retry-After', 5))
                    print(f"  Rate limited, waiting {retry_after}s...")
                    time.sleep(retry_after)
                else:
                    print(f"  HTTP {response.status_code} on attempt {attempt + 1}")
                    if attempt < max_retries - 1:
                        time.sleep(2 ** attempt)
                        
            except requests.RequestException as e:
                print(f"  Request error on attempt {attempt + 1}: {e}")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
        
        # Brief pause between batches to respect rate limits
        if start + batch_size < total:
            time.sleep(1)
    
    return all_results

def annotate_proteins(mapping_file, output_file):
    print(f"Loading UniProt mapping from: {mapping_file}")
    
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\\t')
        print(f"Loaded {len(mapping_df)} sequence mappings")
        
        if len(mapping_df) == 0:
            print("Warning: No mappings found")
            empty_df = pd.DataFrame(columns=[
                'sequence_id', 'uniprot_id', 'protein_name', 'gene_names',
                'organism', 'length', 'go_biological_process', 
                'go_cellular_component', 'go_molecular_function'
            ])
            empty_df.to_csv(output_file, sep='\\t', index=False)
            return
        
        # Get unique IDs, filtering out NaN/empty
        unique_ids = [
            str(uid).strip() for uid in mapping_df['uniprot_id'].unique()
            if pd.notna(uid) and str(uid).strip() != ''
        ]
        
        print(f"Fetching annotations for {len(unique_ids)} unique UniProt IDs (batch mode)...")
        
        # Batch fetch all annotations
        results = batch_fetch_uniprot(unique_ids)
        
        print(f"Successfully retrieved {len(results)} annotations from UniProt")
        
        # Build output: map results back to sequence IDs
        annotations = []
        found = 0
        not_found = 0
        
        for _, row in mapping_df.iterrows():
            seq_id = row['sequence_id']
            uniprot_id = str(row['uniprot_id']).strip()
            
            if uniprot_id in results:
                info = results[uniprot_id]
                found += 1
                annotations.append({
                    'sequence_id': seq_id,
                    'uniprot_id': uniprot_id,
                    'protein_name': info['protein_name'],
                    'gene_names': info['gene_names'],
                    'organism': info['organism'],
                    'length': info['length'],
                    'go_biological_process': info['go_biological_process'],
                    'go_cellular_component': info['go_cellular_component'],
                    'go_molecular_function': info['go_molecular_function']
                })
            else:
                not_found += 1
                annotations.append({
                    'sequence_id': seq_id,
                    'uniprot_id': uniprot_id,
                    'protein_name': 'Unknown protein',
                    'gene_names': '',
                    'organism': 'Unknown organism',
                    'length': '0',
                    'go_biological_process': '',
                    'go_cellular_component': '',
                    'go_molecular_function': ''
                })
        
        # Save results
        annotations_df = pd.DataFrame(annotations)
        
        if len(annotations_df) > 0:
            annotations_df.to_csv(output_file, sep='\\t', index=False)
            print(f"\\nAnnotation complete:")
            print(f"  Total sequences: {len(annotations_df)}")
            print(f"  Annotated: {found}")
            print(f"  Not found: {not_found}")
            print(f"  Saved to: {output_file}")
        else:
            print("No annotations could be retrieved")
            empty_df = pd.DataFrame(columns=[
                'sequence_id', 'uniprot_id', 'protein_name', 'gene_names',
                'organism', 'length', 'go_biological_process', 
                'go_cellular_component', 'go_molecular_function'
            ])
            empty_df.to_csv(output_file, sep='\\t', index=False)
            
    except Exception as e:
        print(f"Error in protein annotation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 annotate_proteins.py <mapping.tsv> <output.tsv>")
        sys.exit(1)
    
    mapping_file = sys.argv[1]
    output_file = sys.argv[2]
    
    annotate_proteins(mapping_file, output_file)
PYEOF

    chmod +x annotate_proteins.py
    python3 annotate_proteins.py ${uniprot_mapping} protein_annotations.tsv
    """
}

process AnnotateGOTerms {
    tag "go_annotation"
    publishDir "${params.outdir}/annotation", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path protein_annotations

    output:
    path "final_annotations.tsv", emit: final_annotations

    script:
    """
    # Create inline Python script for GO term processing
    cat << 'EOF' > annotate_go_terms.py
#!/usr/bin/env python3

import sys
import pandas as pd

def process_go_terms(input_file, output_file):
    print(f"Processing GO terms from: {input_file}")
    
    try:
        df = pd.read_csv(input_file, sep='\\t')
        print(f"Loaded {len(df)} protein annotations")
        
        if len(df) == 0:
            print("Warning: No protein annotations found")
            # Create empty output with expected headers
            empty_df = pd.DataFrame(columns=[
                'sequence_id', 'uniprot_id', 'protein_name', 'gene_names',
                'organism', 'length', 'GO_Biological_Process', 
                'GO_Cellular_Component', 'GO_Molecular_Function'
            ])
            empty_df.to_csv(output_file, sep='\\t', index=False)
            return
        
        # Rename GO columns to match expected output format
        column_mapping = {
            'go_biological_process': 'GO_Biological_Process',
            'go_cellular_component': 'GO_Cellular_Component', 
            'go_molecular_function': 'GO_Molecular_Function'
        }
        
        df = df.rename(columns=column_mapping)
        
        # Clean up GO terms - remove empty entries and clean formatting
        for go_col in ['GO_Biological_Process', 'GO_Cellular_Component', 'GO_Molecular_Function']:
            if go_col in df.columns:
                # Replace NaN with empty string
                df[go_col] = df[go_col].fillna('')
                # Clean up any extra whitespace
                df[go_col] = df[go_col].astype(str).str.strip()
                # Replace 'nan' string with empty string
                df[go_col] = df[go_col].replace('nan', '')
        
        # Save final annotations
        df.to_csv(output_file, sep='\\t', index=False)
        
        # Print summary statistics
        total_sequences = len(df)
        with_bp = len(df[df['GO_Biological_Process'].str.len() > 0])
        with_cc = len(df[df['GO_Cellular_Component'].str.len() > 0]) 
        with_mf = len(df[df['GO_Molecular_Function'].str.len() > 0])
        
        print(f"Final annotation summary:")
        print(f"  Total sequences: {total_sequences}")
        print(f"  With Biological Process GO terms: {with_bp} ({with_bp/total_sequences*100:.1f}%)")
        print(f"  With Cellular Component GO terms: {with_cc} ({with_cc/total_sequences*100:.1f}%)")
        print(f"  With Molecular Function GO terms: {with_mf} ({with_mf/total_sequences*100:.1f}%)")
        
        print(f"Final annotations saved to: {output_file}")
        
    except Exception as e:
        print(f"Error processing GO terms: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 annotate_go_terms.py <input.tsv> <output.tsv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_go_terms(input_file, output_file)
EOF

    # Make script executable and run it
    chmod +x annotate_go_terms.py
    python3 annotate_go_terms.py ${protein_annotations} final_annotations.tsv
    """
}

process GenerateAnnotationReport {
    tag "annotation_report"
    publishDir "${params.outdir}/annotation", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path final_annotations

    output:
    path "annotation_summary.html", emit: report
    path "*.png", optional: true, emit: plots

    script:
    """
    # Create inline R script for annotation reporting
    cat << 'EOF' > generate_annotation_report.R
#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

cat("Generating annotation report from:", input_file, "\\n")

# Read annotation data
tryCatch({
    data <- read_tsv(input_file, show_col_types = FALSE)
    
    if (nrow(data) == 0) {
        cat("No annotation data found. Creating minimal report.\\n")
        
        html_content <- paste(
            "<html><head><title>SCEPTR Annotation Report</title></head>",
            "<body><h1>Annotation Summary</h1>",
            "<p>No protein annotations were generated.</p>",
            "<p>This may indicate:</p>",
            "<ul>",
            "<li>No significant BLAST hits against UniProt</li>",
            "<li>Network issues during annotation retrieval</li>",
            "<li>Very stringent BLAST parameters</li>",
            "</ul></body></html>",
            sep = "\\n"
        )
        
        writeLines(html_content, "annotation_summary.html")
        quit(save = "no", status = 0)
    }
    
    cat("Creating annotation summary visualizations...\\n")
    
    # Calculate summary statistics
    total_seqs <- nrow(data)
    with_bp <- sum(nchar(data\$GO_Biological_Process) > 0, na.rm = TRUE)
    with_cc <- sum(nchar(data\$GO_Cellular_Component) > 0, na.rm = TRUE)
    with_mf <- sum(nchar(data\$GO_Molecular_Function) > 0, na.rm = TRUE)
    with_uniprot <- sum(!is.na(data\$uniprot_id) & nchar(data\$uniprot_id) > 0)
    
    # Create annotation coverage plot
    coverage_data <- data.frame(
        Category = c("UniProt ID", "Biological Process", "Cellular Component", "Molecular Function"),
        Count = c(with_uniprot, with_bp, with_cc, with_mf),
        Percentage = c(with_uniprot/total_seqs*100, with_bp/total_seqs*100, 
                      with_cc/total_seqs*100, with_mf/total_seqs*100)
    )
    
    p1 <- ggplot(coverage_data, aes(x = reorder(Category, Count), y = Count)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        geom_text(aes(label = paste0(Count, " (", round(Percentage, 1), "%)")), 
                 hjust = -0.1) +
        coord_flip() +
        labs(title = "Annotation Coverage Summary",
             x = "Annotation Type", 
             y = "Number of Sequences") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
    
    ggsave("annotation_coverage.png", p1, width = 10, height = 6)
    
    # Create HTML report
    html_content <- paste(
        "<html><head><title>SCEPTR Annotation Report</title></head>",
        "<body>",
        "<h1>SCEPTR Functional Annotation Summary</h1>",
        "<h2>Overview</h2>",
        sprintf("<p><strong>Total sequences processed:</strong> %d</p>", total_seqs),
        sprintf("<p><strong>Sequences with UniProt annotations:</strong> %d (%.1f%%)</p>", 
                with_uniprot, with_uniprot/total_seqs*100),
        "<h2>GO Term Coverage</h2>",
        sprintf("<p><strong>Biological Process terms:</strong> %d (%.1f%%)</p>", 
                with_bp, with_bp/total_seqs*100),
        sprintf("<p><strong>Cellular Component terms:</strong> %d (%.1f%%)</p>", 
                with_cc, with_cc/total_seqs*100),
        sprintf("<p><strong>Molecular Function terms:</strong> %d (%.1f%%)</p>", 
                with_mf, with_mf/total_seqs*100),
        "<h2>Visualizations</h2>",
        "<p>Annotation coverage plot: <a href='annotation_coverage.png'>annotation_coverage.png</a></p>",
        "<h2>Files Generated</h2>",
        "<p>Final annotations file: final_annotations.tsv</p>",
        "<p>This file contains sequence IDs, UniProt annotations, and GO terms ready for enrichment analysis.</p>",
        "</body></html>",
        sep = "\\n"
    )
    
    writeLines(html_content, "annotation_summary.html")
    cat("Annotation report generated successfully.\\n")
    
}, error = function(e) {
    cat("Error in annotation report generation:", conditionMessage(e), "\\n")
    
    # Create error report
    html_content <- paste(
        "<html><head><title>SCEPTR Annotation Report</title></head>",
        "<body><h1>Annotation Report</h1>",
        "<p>Error in report generation: ", conditionMessage(e), "</p>",
        "</body></html>",
        sep = ""
    )
    writeLines(html_content, "annotation_summary.html")
})
EOF

    # Make R script executable and run it
    chmod +x generate_annotation_report.R
    Rscript generate_annotation_report.R ${final_annotations} .
    """
}
