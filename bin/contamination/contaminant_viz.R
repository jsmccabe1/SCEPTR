#!/usr/bin/env Rscript

# SCEPTR DeCon Contamination Visualisation
# Generates HTML report and diagnostic plots from contamination filtering results.
# Reads filtering parameters from JSON to ensure report values match actual thresholds.

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript contaminant_viz.R <contaminant_details.csv> [output_dir]")
}

input_file <- args[1]
output_dir <- if (length(args) > 1) args[2] else "."
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read parameters from JSON (written by filter_contaminants.py) ----
params <- list(
  bacterial_threshold = 70.0,
  viral_threshold = 50.0,
  fungal_threshold = 70.0,
  eukaryotic_threshold = 90.0,
  coverage_threshold = 30.0,
  evalue_threshold = 1e-3
)

# Try to read JSON params - use base R parsing to avoid jsonlite dependency
read_json_params <- function(filepath) {
  tryCatch({
    lines <- readLines(filepath, warn = FALSE)
    json_str <- paste(lines, collapse = "")
    # Simple key-value extraction for flat JSON
    keys <- regmatches(json_str, gregexpr('"([^"]+)"\\s*:', json_str))[[1]]
    keys <- gsub('[":[:space:]]', '', keys)
    vals <- regmatches(json_str, gregexpr(':\\s*([^,}]+)', json_str))[[1]]
    vals <- gsub('^:\\s*', '', vals)
    vals <- gsub('"', '', vals)
    vals <- trimws(vals)
    result <- as.list(setNames(as.numeric(vals), keys))
    result[!is.na(result)]
  }, error = function(e) NULL)
}

params_file <- file.path(dirname(input_file), "filtering_params.json")
if (!file.exists(params_file)) params_file <- "filtering_params.json"
if (file.exists(params_file)) {
  parsed <- read_json_params(params_file)
  if (!is.null(parsed)) {
    for (nm in names(parsed)) params[[nm]] <- parsed[[nm]]
    message("Read filtering parameters from: ", params_file)
  }
} else {
  message("Warning: filtering_params.json not found, using defaults")
}

# ---- Read contaminant data ----
message("Reading contaminant data from ", input_file)
contam_data <- tryCatch({
  read_csv(input_file, show_col_types = FALSE)
}, error = function(e) {
  message("Error reading input file: ", e$message)
  return(NULL)
})

if (is.null(contam_data) || nrow(contam_data) == 0) {
  message("No contaminant data found. Creating minimal report.")
  html_content <- paste0('<!DOCTYPE html>
<html><head><title>SCEPTR DeCon Report</title>
<style>body{font-family:Arial,sans-serif;margin:40px;}
.success{background:#d4edda;padding:20px;border-radius:5px;border-left:5px solid #28a745;}</style>
</head><body>
<h1>SCEPTR DeCon Contamination Analysis</h1>
<div class="success"><h2>No Contamination Detected</h2>
<p>No sequences exceeded the filtering thresholds.</p></div>
<p><em>Generated ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</em></p>
</body></html>')
  writeLines(html_content, file.path(output_dir, "contaminant_analysis_report.html"))
  quit(status = 0)
}

message("Found ", nrow(contam_data), " contaminant entries")

# ---- Ensure required columns ----
if (!"category" %in% colnames(contam_data)) {
  stop("Input CSV must contain a 'category' column")
}

# ---- Category colours ----
category_colors <- c(
  "Bacteria" = "#E74C3C",
  "Viruses"  = "#9B59B6",
  "Fungi"    = "#F39C12",
  "Mammals"  = "#E91E63",
  "Plants"   = "#27AE60",
  "Other"    = "#3498DB"
)

# Category summary
category_summary <- contam_data %>%
  group_by(category) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

# ---- Plot 1: Category distribution ----
p1 <- ggplot(category_summary, aes(x = reorder(category, -count), y = count, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, size = 4) +
  scale_fill_manual(values = category_colors) +
  labs(title = "Contaminant Sequences by Category",
       x = NULL, y = "Number of sequences") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "contaminant_category_distribution.png"), p1, width = 8, height = 5, dpi = 150)

# ---- Plot 2: Identity distribution by category ----
if (nrow(contam_data) >= 3) {
  p2 <- ggplot(contam_data, aes(x = percent_identity, fill = category)) +
    geom_histogram(binwidth = 2, alpha = 0.8, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = category_colors) +
    labs(title = "Percent Identity Distribution",
         x = "Percent Identity (%)", y = "Count") +
    facet_wrap(~category, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "contaminant_percent_identity_by_category.png"), p2, width = 10, height = 7, dpi = 150)
} else {
  # With very few points, use a dot plot instead
  p2 <- ggplot(contam_data, aes(x = category, y = percent_identity, colour = category)) +
    geom_jitter(width = 0.2, size = 3) +
    scale_colour_manual(values = category_colors) +
    labs(title = "Percent Identity by Category",
         x = NULL, y = "Percent Identity (%)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "contaminant_percent_identity_by_category.png"), p2, width = 8, height = 5, dpi = 150)
}

# ---- Plot 3: E-value distribution ----
p3 <- ggplot(contam_data, aes(x = category, y = -log10(e_value), fill = category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  scale_fill_manual(values = category_colors) +
  labs(title = "E-value Distribution by Category",
       x = NULL, y = expression(-log[10](E-value))) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "contaminant_evalue_distribution.png"), p3, width = 8, height = 5, dpi = 150)

# ---- Plot 4: Identity vs E-value scatter with threshold lines ----
# Build threshold data for vertical lines
threshold_df <- data.frame(
  category = c("Bacteria", "Viruses", "Fungi", "Mammals", "Plants"),
  threshold = c(params$bacterial_threshold, params$viral_threshold,
                params$fungal_threshold, params$eukaryotic_threshold,
                params$eukaryotic_threshold)
) %>% distinct(threshold, .keep_all = TRUE)

p4 <- ggplot(contam_data, aes(x = percent_identity, y = -log10(e_value), colour = category)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = category_colors) +
  geom_hline(yintercept = -log10(params$evalue_threshold),
             linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_vline(data = threshold_df, aes(xintercept = threshold),
             linetype = "dotted", colour = "grey50", linewidth = 0.5) +
  labs(title = "Contaminant Hits: Identity vs E-value",
       x = "Percent Identity (%)",
       y = expression(-log[10](E-value)),
       colour = "Category") +
  annotate("text", x = max(contam_data$percent_identity, na.rm = TRUE),
           y = -log10(params$evalue_threshold) + 0.5,
           label = paste0("E-value = ", params$evalue_threshold),
           hjust = 1, size = 3, colour = "grey40") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "contaminant_identity_vs_evalue.png"), p4, width = 10, height = 7, dpi = 150)

# ---- Plot 5: Coverage vs Identity (new diagnostic plot) ----
if ("query_coverage" %in% colnames(contam_data)) {
  p5 <- ggplot(contam_data, aes(x = percent_identity, y = query_coverage, colour = category)) +
    geom_point(alpha = 0.8, size = 3) +
    scale_colour_manual(values = category_colors) +
    geom_hline(yintercept = params$coverage_threshold,
               linetype = "dashed", colour = "grey50") +
    labs(title = "Query Coverage vs Percent Identity",
         x = "Percent Identity (%)",
         y = "Query Coverage (%)",
         colour = "Category") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "contaminant_coverage_vs_identity.png"), p5, width = 10, height = 7, dpi = 150)
}

# ---- Generate HTML report ----
# Build threshold table rows
threshold_table_rows <- paste0(
  '<tr><td style="color:', category_colors["Bacteria"], '">Bacteria</td>',
  '<td>', params$bacterial_threshold, '%</td></tr>\n',
  '<tr><td style="color:', category_colors["Viruses"], '">Viruses</td>',
  '<td>', params$viral_threshold, '%</td></tr>\n',
  '<tr><td style="color:', category_colors["Fungi"], '">Fungi</td>',
  '<td>', params$fungal_threshold, '%</td></tr>\n',
  '<tr><td style="color:', category_colors["Mammals"], '">Mammals</td>',
  '<td>', params$eukaryotic_threshold, '%</td></tr>\n',
  '<tr><td style="color:', category_colors["Plants"], '">Plants</td>',
  '<td>', params$eukaryotic_threshold, '%</td></tr>\n'
)

# Build category summary rows
summary_rows <- ""
for (i in seq_len(nrow(category_summary))) {
  pct <- round(category_summary$count[i] / sum(category_summary$count) * 100, 1)
  cat_colour <- category_colors[category_summary$category[i]]
  if (is.na(cat_colour)) cat_colour <- "#333"
  summary_rows <- paste0(summary_rows,
    '<tr><td style="color:', cat_colour, ';font-weight:bold">',
    category_summary$category[i], '</td><td>',
    category_summary$count[i], '</td><td>', pct, '%</td></tr>\n')
}

# Coverage plot reference (only if column exists)
coverage_plot_html <- ""
if ("query_coverage" %in% colnames(contam_data)) {
  coverage_plot_html <- '
    <div class="plot-card">
        <h3>Query Coverage vs Percent Identity</h3>
        <img src="contaminant_coverage_vs_identity.png" alt="Coverage vs Identity" width="800">
        <p>Sequences below the dashed line did not meet the minimum coverage threshold and would not be flagged.</p>
    </div>'
}

html_content <- paste0('<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>SCEPTR DeCon Analysis</title>
    <style>
        * { box-sizing: border-box; }
        body { font-family: "Segoe UI", Arial, sans-serif; margin: 0; padding: 0; background: #f5f6fa; color: #2c3e50; }
        .header { background: linear-gradient(135deg, #2c3e50, #3498db); color: white; padding: 40px; }
        .header h1 { margin: 0 0 8px 0; font-size: 28px; }
        .header p { margin: 0; opacity: 0.85; font-size: 14px; }
        .content { max-width: 1100px; margin: 0 auto; padding: 30px; }
        .card { background: white; border-radius: 8px; padding: 24px; margin-bottom: 24px;
                box-shadow: 0 1px 4px rgba(0,0,0,0.08); }
        .plot-card { background: white; border-radius: 8px; padding: 24px; margin-bottom: 24px;
                     box-shadow: 0 1px 4px rgba(0,0,0,0.08); text-align: center; }
        .plot-card img { max-width: 100%; height: auto; border-radius: 4px; }
        .plot-card p { color: #666; font-size: 13px; margin-top: 12px; }
        h2 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 8px; }
        h3 { color: #34495e; }
        table { border-collapse: collapse; width: 100%; margin: 12px 0; }
        th { background: #f8f9fa; padding: 10px 12px; text-align: left; font-weight: 600;
             border-bottom: 2px solid #dee2e6; }
        td { padding: 8px 12px; border-bottom: 1px solid #eee; }
        .method-box { background: #f8f9fa; border-left: 4px solid #3498db; padding: 16px 20px;
                      border-radius: 0 8px 8px 0; margin: 16px 0; }
        .method-box p { margin: 6px 0; }
        .footer { text-align: center; color: #999; font-size: 12px; padding: 20px; }
    </style>
</head>
<body>
    <div class="header">
        <h1>SCEPTR DeCon Contamination Analysis</h1>
        <p>Taxonomy-aware contamination detection and filtering</p>
    </div>
    <div class="content">

    <div class="card">
        <h2>Filtering Parameters</h2>
        <p>Identity thresholds are applied per-category. Coverage and E-value filters apply globally.</p>
        <table>
            <tr><th>Category</th><th>Identity Threshold</th></tr>
            ', threshold_table_rows, '
        </table>
        <p style="margin-top:12px">
            <strong>Minimum query coverage:</strong> ', params$coverage_threshold, '%<br>
            <strong>Maximum E-value:</strong> ', params$evalue_threshold, '
        </p>
    </div>

    <div class="card">
        <h2>Summary</h2>
        <p>Total contaminant sequences identified: <strong>', nrow(contam_data), '</strong></p>
        <table>
            <tr><th>Category</th><th>Count</th><th>Percentage</th></tr>
            ', summary_rows, '
        </table>
    </div>

    <h2>Visualisations</h2>

    <div class="plot-card">
        <h3>Contaminant Distribution by Category</h3>
        <img src="contaminant_category_distribution.png" alt="Category Distribution" width="800">
    </div>

    <div class="plot-card">
        <h3>Percent Identity Distribution</h3>
        <img src="contaminant_percent_identity_by_category.png" alt="Identity Distribution" width="800">
    </div>

    <div class="plot-card">
        <h3>E-value Distribution by Category</h3>
        <img src="contaminant_evalue_distribution.png" alt="E-value Distribution" width="800">
    </div>

    <div class="plot-card">
        <h3>Identity vs E-value</h3>
        <img src="contaminant_identity_vs_evalue.png" alt="Identity vs E-value" width="800">
        <p>Dotted vertical lines show category-specific identity thresholds. Dashed horizontal line shows the E-value cutoff.</p>
    </div>

    ', coverage_plot_html, '

    <div class="card">
        <h2>Methods</h2>
        <div class="method-box">
            <p>Predicted protein sequences were searched against a curated contaminant database
            using DIAMOND BLASTp in sensitive mode. Hits were classified into taxonomic categories
            (Bacteria, Viruses, Fungi, Mammals, Plants) using word-boundary pattern matching on
            subject descriptions. Category-specific identity thresholds were applied to account for
            differing evolutionary distances: prokaryotic contaminants are detectable at lower identity
            due to ancient divergence from eukaryotes, while eukaryotic hits require near-identical
            matches to distinguish genuine contamination from conserved orthologs. All hits were
            additionally filtered by minimum query coverage and maximum E-value.</p>
            <p><strong>Conserved protein family protection:</strong> Hits matching known ultra-conserved
            protein families (ubiquitin, histones, tubulins, actins, ribosomal proteins, elongation
            factors, heat shock proteins, calmodulin, and cyclophilins) are subject to an elevated
            identity threshold (&ge;99%). These families show &gt;90% cross-kingdom identity due to
            ancient functional constraints (Ciccarelli et al., 2006), and moderate-identity hits
            represent genuine orthologs rather than contamination.</p>
        </div>
    </div>

    </div>
    <div class="footer">
        Generated by SCEPTR DeCon &middot; ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '
    </div>
</body>
</html>')

writeLines(html_content, file.path(output_dir, "contaminant_analysis_report.html"))
message("Report complete: ", file.path(output_dir, "contaminant_analysis_report.html"))
