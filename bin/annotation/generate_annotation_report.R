#!/usr/bin/env Rscript

# Script to generate annotation summary report
# Usage: Rscript generate_annotation_report.R <final_annotations.tsv> <output_directory>

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript generate_annotation_report.R <final_annotations.tsv> [output_directory]")
}

input_file <- args[1]
output_dir <- if (length(args) > 1) args[2] else "."

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# Read annotations data
message("Reading annotations from ", input_file)
data <- tryCatch({
    read_tsv(input_file, show_col_types = FALSE)
}, error = function(e) {
    message("Error reading input file: ", e$message)
    # Create minimal DataFrame for report generation
    data.frame(
        Sequence_ID = character(0),
        UniProt_ID = character(0),
        Protein_Name = character(0),
        GO_BP = character(0),
        GO_CC = character(0),
        GO_MF = character(0)
    )
})

# Count annotations
total_sequences <- nrow(data)
message("Total sequences: ", total_sequences)

sequences_with_uniprot <- sum(!is.na(data$UniProt_ID))
sequences_with_protein_name <- sum(!is.na(data$Protein_Name))
sequences_with_function <- sum(!is.na(data$Function))
sequences_with_go_bp <- sum(!is.na(data$GO_BP))
sequences_with_go_cc <- sum(!is.na(data$GO_CC))
sequences_with_go_mf <- sum(!is.na(data$GO_MF))
sequences_with_any_go <- sum(!is.na(data$GO_BP) | !is.na(data$GO_CC) | !is.na(data$GO_MF))

# Calculate percentages
pct_with_uniprot <- round(sequences_with_uniprot / total_sequences * 100, 2)
pct_with_protein_name <- round(sequences_with_protein_name / total_sequences * 100, 2)
pct_with_function <- round(sequences_with_function / total_sequences * 100, 2)
pct_with_go_bp <- round(sequences_with_go_bp / total_sequences * 100, 2)
pct_with_go_cc <- round(sequences_with_go_cc / total_sequences * 100, 2)
pct_with_go_mf <- round(sequences_with_go_mf / total_sequences * 100, 2)
pct_with_any_go <- round(sequences_with_any_go / total_sequences * 100, 2)

message("Sequences with UniProt IDs: ", sequences_with_uniprot, " (", pct_with_uniprot, "%)")
message("Sequences with any GO terms: ", sequences_with_any_go, " (", pct_with_any_go, "%)")

# Create a data frame for the summary plot
plot_data <- data.frame(
    Category = c("UniProt ID", "Protein Name", "Function", "GO: Biological Process", 
                "GO: Cellular Component", "GO: Molecular Function", "Any GO Annotation"),
    Count = c(sequences_with_uniprot, sequences_with_protein_name, sequences_with_function, 
            sequences_with_go_bp, sequences_with_go_cc, sequences_with_go_mf, sequences_with_any_go),
    Percentage = c(pct_with_uniprot, pct_with_protein_name, pct_with_function, 
                  pct_with_go_bp, pct_with_go_cc, pct_with_go_mf, pct_with_any_go)
)

# Create annotation coverage plot
p1 <- ggplot(plot_data, aes(x = reorder(Category, -Percentage), y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Annotation Coverage", x = "", y = "Percentage of Sequences (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Blues")

# Save plot
plot_file <- file.path(output_dir, "annotation_coverage.png")
message("Saving plot to: ", plot_file)
ggsave(plot_file, p1, width = 10, height = 6)

# Generate HTML report
html_file <- file.path(output_dir, "annotation_summary.html")
message("Generating HTML report: ", html_file)

html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
    <title>SCEPTR Functional Annotation Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
        h1, h2, h3 { color: #2c3e50; }
        .summary { background-color: #f8f9fa; padding: 20px; border-radius: 5px; margin-bottom: 30px; }
        .plot { margin: 30px 0; background-color: #fff; padding: 20px; border-radius: 5px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .success { color: #28a745; }
        .warning { color: #ffc107; }
        .danger { color: #dc3545; }
    </style>
</head>
<body>
    <h1>SCEPTR Functional Annotation Report</h1>
    
    <div class="summary">
        <h2>Summary</h2>
        <p>Total sequences analysed: ', total_sequences, '</p>
        
        <table>
            <tr>
                <th>Annotation Type</th>
                <th>Count</th>
                <th>Percentage</th>
            </tr>
            <tr>
                <td>Sequences with UniProt matches</td>
                <td>', sequences_with_uniprot, '</td>
                <td>', pct_with_uniprot, '%</td>
            </tr>
            <tr>
                <td>Sequences with Protein Names</td>
                <td>', sequences_with_protein_name, '</td>
                <td>', pct_with_protein_name, '%</td>
            </tr>
            <tr>
                <td>Sequences with Function Annotations</td>
                <td>', sequences_with_function, '</td>
                <td>', pct_with_function, '%</td>
            </tr>
            <tr>
                <td>Sequences with GO Biological Process</td>
                <td>', sequences_with_go_bp, '</td>
                <td>', pct_with_go_bp, '%</td>
            </tr>
            <tr>
                <td>Sequences with GO Cellular Component</td>
                <td>', sequences_with_go_cc, '</td>
                <td>', pct_with_go_cc, '%</td>
            </tr>
            <tr>
                <td>Sequences with GO Molecular Function</td>
                <td>', sequences_with_go_mf, '</td>
                <td>', pct_with_go_mf, '%</td>
            </tr>
            <tr>
                <td>Sequences with any GO Annotation</td>
                <td>', sequences_with_any_go, '</td>
                <td>', pct_with_any_go, '%</td>
            </tr>
        </table>
    </div>
    
    <div class="plot">
        <h2>Annotation Coverage</h2>
        <p>The following plot shows the percentage of sequences with each type of annotation:</p>
        <img src="annotation_coverage.png" alt="Annotation Coverage" width="800">
    </div>
    
    <div class="summary">
        <h2>Interpretation</h2>
        <p>This report summarises the functional annotation process results:</p>
        <ul>
            <li><strong>UniProt matches:</strong> The percentage of sequences that matched entries in the UniProt database.</li>
            <li><strong>Protein Names:</strong> The percentage of sequences with assigned protein names from UniProt.</li>
            <li><strong>Function Annotations:</strong> The percentage of sequences with functional descriptions.</li>
            <li><strong>GO Terms:</strong> The percentage of sequences annotated with Gene Ontology terms in each category.</li>
        </ul>
        
        <h3>Annotation Quality Assessment</h3>
        <p>', 
        if(pct_with_uniprot > 80) {
            '<span class="success">High annotation coverage (>80%)</span>'
        } else if(pct_with_uniprot > 50) {
            '<span class="warning">Moderate annotation coverage (50-80%)</span>'
        } else {
            '<span class="danger">Low annotation coverage (<50%)</span>'
        },
        ' - ', pct_with_uniprot, '% of sequences have UniProt matches.</p>
        
        <p>',
        if(pct_with_any_go > 70) {
            '<span class="success">Good GO term coverage (>70%)</span>'
        } else if(pct_with_any_go > 40) {
            '<span class="warning">Moderate GO term coverage (40-70%)</span>'
        } else {
            '<span class="danger">Low GO term coverage (<40%)</span>'
        },
        ' - ', pct_with_any_go, '% of sequences have at least one GO term.</p>
    </div>
    
    <p><em>Report generated by SCEPTR on ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</em></p>
</body>
</html>
')

# Write HTML report
writeLines(html_content, html_file)

message("Annotation report generated successfully.")
