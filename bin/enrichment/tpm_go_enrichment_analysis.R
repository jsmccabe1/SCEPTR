#!/usr/bin/env Rscript

# SCEPTR TPM-Based GO Enrichment Analysis
#
# Tier-based enrichment analysis for single-sample transcriptomics.
# Uses topGO (weight01 algorithm) with expression-tier stratification
# to characterise expression-dependent functional enrichment.
#
# Statistical enhancements:
#   - Fold enrichment with exact binomial confidence intervals
#   - Bootstrap stability analysis with configurable iterations
#   - Global FDR correction across all tier x category combinations
#   - Effect size classification (heuristic, user-configurable)
#
# Note on weight01 + BH correction:
#   topGO's weight01 algorithm produces decorrelated p-values that account
#   for the GO DAG hierarchy. Applying BH FDR on top is slightly conservative
#   (Alexa & Rahnenführer, 2009). We retain BH for consistency with standard
#   practice, noting that reported FDR values may be somewhat conservative.

options(warn = 1)
options(error = function() {
  cat("ERROR: ", geterrmessage(), "\n")
  cat("Stack trace:\n")
  print(sys.calls())
  quit(status = 1)
})

# Library loading
safe_library_load <- function(lib_name) {
  tryCatch({
    suppressPackageStartupMessages(library(lib_name, character.only = TRUE))
    TRUE
  }, error = function(e) {
    message(paste("Error loading package:", lib_name, "-", e$message))
    FALSE
  })
}

required_libs <- c("optparse", "dplyr", "ggplot2", "stringr", "readr", "reshape2", "RColorBrewer")
missing_libs <- c()
for (lib in required_libs) {
  if (!safe_library_load(lib)) missing_libs <- c(missing_libs, lib)
}
if (!safe_library_load("topGO")) missing_libs <- c(missing_libs, "topGO")

if (length(missing_libs) > 0) {
  message("Missing required packages: ", paste(missing_libs, collapse = ", "))
  writeLines("GO.ID\tTerm\tCategory\tTier\tFold_Enrichment\tp.value\tFDR",
             "sceptr_go_enrichment_results.tsv")
  writeLines('<html><body><h1>Analysis failed: Missing R packages</h1></body></html>',
             "sceptr_tpm_enrichment.html")
  quit(status = 1)
}

# ---- Command line options ----
option_list <- list(
  make_option("--go_data", type = "character", help = "Path to long-format GO terms CSV"),
  make_option("--expression", type = "character", help = "Path to expression data TSV"),
  make_option("--output_prefix", type = "character", default = "sceptr"),
  make_option("--tpm_threshold", type = "double", default = 1.0),
  make_option("--tiers", type = "character", default = "50,100,250,500"),
  make_option("--pvalue", type = "double", default = 0.05),
  make_option("--bootstrap_n", type = "integer", default = 50,
              help = "Bootstrap iterations (default: 50)"),
  make_option("--min_fold_enrichment", type = "double", default = 1.5),
  make_option("--confidence_level", type = "double", default = 0.95)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Statistical functions ----

calculate_fold_enrichment_ci <- function(observed, expected, total_annotated,
                                          total_genes, confidence = 0.95) {
  if (expected == 0 || total_annotated == 0 || observed == 0) {
    return(c(fold_enrichment = 0, ci_lower = 0, ci_upper = Inf))
  }
  
  fold_enrichment <- observed / expected
  
  n_trial <- min(total_annotated, total_genes)
  
  if (observed >= n_trial) {
    ci_result <- binom.test(observed, n_trial, conf.level = confidence)
    p_lower <- ci_result$conf.int[1]
    return(c(fold_enrichment = fold_enrichment,
             ci_lower = (p_lower * n_trial) / expected,
             ci_upper = Inf))
  }
  
  ci_result <- binom.test(observed, n_trial, conf.level = confidence)
  p_lower <- ci_result$conf.int[1]
  p_upper <- ci_result$conf.int[2]
  
  ci_lower <- (p_lower * n_trial) / expected
  ci_upper <- (p_upper * n_trial) / expected
  
  return(c(fold_enrichment = fold_enrichment, ci_lower = ci_lower, ci_upper = ci_upper))
}


bootstrap_stability_analysis <- function(expr_data, go_data, tier, category, n_bootstrap = 50) {
  message(paste("  Bootstrap stability:", category, "tier", tier,
                "(", n_bootstrap, "iterations)"))
  
  go_category <- go_data[go_data$Category == category, ]
  if (nrow(go_category) == 0) return(data.frame())
  
  unique_go_ids <- unique(go_category$GO_ID)
  n_go <- length(unique_go_ids)
  original_n <- nrow(expr_data)
  
  if (tier > original_n) return(data.frame())
  
  go_gene_sets <- split(go_category$UniProt_ID, go_category$GO_ID)
  go_term_sizes <- vapply(go_gene_sets, length, integer(1))
  
  fold_matrix <- matrix(NA_real_, nrow = n_go, ncol = n_bootstrap)
  rownames(fold_matrix) <- unique_go_ids
  
  all_ids <- expr_data$UniProt_ID
  all_tpm <- expr_data$TPM
  
  for (i in seq_len(n_bootstrap)) {
    boot_idx <- sample.int(original_n, size = original_n, replace = TRUE)
    boot_tpm <- all_tpm[boot_idx]
    boot_ids <- all_ids[boot_idx]
    
    top_idx <- order(boot_tpm, decreasing = TRUE)[seq_len(min(tier, length(boot_tpm)))]
    tier_genes <- unique(boot_ids[top_idx])
    
    actual_tier_size <- length(tier_genes)
    
    for (j in seq_len(n_go)) {
      go_id <- unique_go_ids[j]
      observed <- sum(go_gene_sets[[go_id]] %in% tier_genes)
      expected <- (go_term_sizes[go_id] * actual_tier_size) / original_n
      
      if (expected > 0) {
        fold_matrix[j, i] <- observed / expected
      }
    }
  }
  
  has_data <- rowSums(!is.na(fold_matrix)) > 0
  if (!any(has_data)) return(data.frame())
  
  fold_sub <- fold_matrix[has_data, , drop = FALSE]
  
  means <- rowMeans(fold_sub, na.rm = TRUE)
  sds <- apply(fold_sub, 1, sd, na.rm = TRUE)
  cvs <- ifelse(means > 0, sds / means, NA_real_)
  
  data.frame(
    GO_ID = rownames(fold_sub),
    Bootstrap_Mean_FC = round(means, 3),
    Bootstrap_SD_FC = round(sds, 3),
    Bootstrap_CV = round(cvs, 3),
    Stability_Score = round(1 / (1 + cvs), 3),
    Bootstrap_N = n_bootstrap,
    stringsAsFactors = FALSE
  )
}


# ---- Read and validate input data ----

make_empty_outputs <- function(prefix, msg) {
  writeLines(paste0("GO.ID\tTerm\tCategory\tTier\tFold_Enrichment\tp.value\tFDR"),
             paste0(prefix, "_go_enrichment_results.tsv"))
  writeLines(paste0('<html><body><h1>GO Enrichment Analysis</h1><p>', msg, '</p></body></html>'),
             paste0(prefix, "_tpm_enrichment.html"))
}

message("Reading input files...")

go_data <- tryCatch(read_csv(opt$go_data, show_col_types = FALSE), error = function(e) NULL)
if (is.null(go_data) || nrow(go_data) == 0) {
  make_empty_outputs(opt$output_prefix, "No GO term data found")
  quit(status = 1)
}

expr_data <- tryCatch({
  if (grepl("\\.csv$", opt$expression)) read_csv(opt$expression, show_col_types = FALSE)
  else read_tsv(opt$expression, show_col_types = FALSE)
}, error = function(e) NULL)

if (is.null(expr_data) || nrow(expr_data) == 0) {
  make_empty_outputs(opt$output_prefix, "No expression data found")
  quit(status = 1)
}

message(paste("GO data:", nrow(go_data), "rows, columns:", paste(colnames(go_data), collapse = ", ")))
message(paste("Expression data:", nrow(expr_data), "rows, columns:", paste(colnames(expr_data), collapse = ", ")))

# ---- Column standardisation ----

standardise_col <- function(df, target, alternatives) {
  if (target %in% colnames(df)) return(df)
  match <- colnames(df)[tolower(colnames(df)) %in% tolower(alternatives)]
  if (length(match) > 0) {
    message(paste("  Renaming", match[1], "->", target))
    df[[target]] <- df[[match[1]]]
  }
  df
}

go_data <- standardise_col(go_data, "GO_ID", c("go_id", "goid"))
go_data <- standardise_col(go_data, "UniProt_ID", c("uniprot_id", "uniprotid", "protein_id"))
go_data <- standardise_col(go_data, "GO_Term", c("go_term", "goterm", "term"))
go_data <- standardise_col(go_data, "Category", c("category", "go_category", "gocategory"))

expr_data <- standardise_col(expr_data, "UniProt_ID", c("uniprot_id", "uniprotid", "protein_id"))
expr_data <- standardise_col(expr_data, "TPM", c("tpm", "fpkm", "expression"))

for (col in c("GO_ID", "UniProt_ID", "GO_Term", "Category")) {
  if (!col %in% colnames(go_data)) stop(paste("Missing required GO data column:", col))
}
for (col in c("UniProt_ID", "TPM")) {
  if (!col %in% colnames(expr_data)) stop(paste("Missing required expression column:", col))
}

# ---- Filter and prepare ----

message(paste("Filtering with TPM threshold:", opt$tpm_threshold))
expr_filtered <- expr_data %>%
  filter(TPM >= opt$tpm_threshold, !is.na(UniProt_ID)) %>%
  arrange(desc(TPM))

n_before <- nrow(expr_filtered)
expr_filtered <- expr_filtered %>% distinct(UniProt_ID, .keep_all = TRUE)
n_after <- nrow(expr_filtered)
if (n_before > n_after) {
  message(paste("Deduplicated UniProt_IDs:", n_before, "->", n_after))
}

total_genes <- nrow(expr_filtered)
message(paste("Genes above TPM threshold:", total_genes))

if (total_genes == 0) {
  make_empty_outputs(opt$output_prefix,
                     paste0("No genes pass TPM threshold of ", opt$tpm_threshold))
  quit(status = 0)
}

go_filtered <- go_data %>% filter(UniProt_ID %in% expr_filtered$UniProt_ID)
message(paste("GO term mappings for expressed genes:", nrow(go_filtered)))

if (nrow(go_filtered) == 0) {
  make_empty_outputs(opt$output_prefix, "No GO terms mapped to expressed genes")
  quit(status = 0)
}

tier_values <- as.numeric(unlist(strsplit(opt$tiers, ",")))
categories <- unique(go_filtered$Category)
message(paste("Categories:", paste(categories, collapse = ", ")))
message(paste("Tiers:", paste(tier_values, collapse = ", ")))

# ---- Enrichment analysis ----

perform_tier_enrichment <- function(go_data_filtered, expr_filtered, tier, category_name) {
  if (tier > nrow(expr_filtered)) {
    message(paste("  Tier", tier, "exceeds gene count, skipping"))
    return(NULL)
  }
  
  message(paste("  Enrichment: top", tier, "genes,", category_name))
  
  go_category <- go_data_filtered %>% filter(Category == category_name)
  if (nrow(go_category) == 0) return(NULL)
  
  gene_to_go <- split(go_category$GO_ID, go_category$UniProt_ID)
  
  gene_list <- rep(0, nrow(expr_filtered))
  names(gene_list) <- expr_filtered$UniProt_ID
  gene_list[seq_len(min(tier, length(gene_list)))] <- 1
  gene_list <- factor(gene_list)
  
  tryCatch({
    go_data_obj <- new("topGOdata",
                       ontology = category_name,
                       allGenes = gene_list,
                       annot = annFUN.gene2GO,
                       gene2GO = gene_to_go)
    
    result <- runTest(go_data_obj, algorithm = "weight01", statistic = "fisher")
    
    if (length(result@score) == 0) return(NULL)
    
    go_results <- GenTable(go_data_obj,
                           weight = result,
                           orderBy = "weight",
                           ranksOf = "weight",
                           topNodes = min(100, length(result@score)),
                           numChar = 100)
    
    # Coerce topGO character columns to numeric
    go_results$Expected <- as.numeric(gsub("< *", "", go_results$Expected))
    go_results$Significant <- as.numeric(go_results$Significant)
    go_results$Annotated <- as.numeric(go_results$Annotated)
    
    go_results$p.value <- result@score[go_results$GO.ID]
    
    # NOTE: FDR will be applied globally AFTER combining all results
    
    ci_results <- t(sapply(seq_len(nrow(go_results)), function(i) {
      calculate_fold_enrichment_ci(
        observed = go_results$Significant[i],
        expected = go_results$Expected[i],
        total_annotated = go_results$Annotated[i],
        total_genes = total_genes,
        confidence = opt$confidence_level
      )
    }))
    
    go_results$Fold_Enrichment <- ci_results[, "fold_enrichment"]
    go_results$FC_CI_Lower <- ci_results[, "ci_lower"]
    go_results$FC_CI_Upper <- ci_results[, "ci_upper"]
    
    if (opt$bootstrap_n > 0) {
      stability <- bootstrap_stability_analysis(
        expr_filtered, go_category, tier, category_name, opt$bootstrap_n
      )
      if (nrow(stability) > 0) {
        go_results <- go_results %>%
          left_join(stability, by = c("GO.ID" = "GO_ID"))
      }
    }
    
    go_results$Effect_Size_Category <- case_when(
      go_results$Fold_Enrichment >= 3.0 ~ "Large",
      go_results$Fold_Enrichment >= 2.0 ~ "Medium",
      go_results$Fold_Enrichment >= opt$min_fold_enrichment ~ "Small",
      TRUE ~ "Negligible"
    )
    
    go_results$Practically_Significant <- (
      go_results$Fold_Enrichment >= opt$min_fold_enrichment &
      go_results$FC_CI_Lower > 1.0
    )
    
    go_results$Category <- category_name
    go_results$Tier <- tier
    
    return(go_results)
  }, error = function(e) {
    message(paste("  topGO error:", e$message))
    return(NULL)
  })
}

message("\nStarting GO enrichment analysis...")
results_list <- list()
idx <- 0

for (category in categories) {
  for (tier in tier_values) {
    res <- perform_tier_enrichment(go_filtered, expr_filtered, tier, category)
    if (!is.null(res) && nrow(res) > 0) {
      idx <- idx + 1
      results_list[[idx]] <- res
    }
  }
}

all_results <- if (length(results_list) > 0) do.call(rbind, results_list) else data.frame()

# ---- GLOBAL FDR correction ----
if (nrow(all_results) > 0) {
  message(paste("\nGlobal BH FDR correction across", nrow(all_results), "tests"))
  all_results$FDR <- p.adjust(all_results$p.value, method = "BH")
  
  priority_cols <- c("GO.ID", "Term", "Category", "Tier",
                     "Fold_Enrichment", "FC_CI_Lower", "FC_CI_Upper",
                     "Effect_Size_Category", "Practically_Significant",
                     "Annotated", "Significant", "Expected",
                     "p.value", "FDR",
                     "Bootstrap_Mean_FC", "Bootstrap_SD_FC", "Stability_Score", "Bootstrap_N")
  
  available_cols <- intersect(priority_cols, colnames(all_results))
  extra_cols <- setdiff(colnames(all_results), priority_cols)
  all_results <- all_results[, c(available_cols, extra_cols)]
  all_results <- all_results %>% arrange(p.value, desc(Fold_Enrichment))
}

# ---- Save results ----
output_file <- paste0(opt$output_prefix, "_go_enrichment_results.tsv")
write.table(all_results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
message(paste("Saved", nrow(all_results), "results to", output_file))

# ---- Visualisations ----
message("Generating visualisations...")
plot_files <- list()

if (nrow(all_results) > 0) {
  significant <- all_results %>% filter(p.value <= opt$pvalue)
  
  if (nrow(significant) > 0) {
    for (cat in categories) {
      cat_results <- significant %>% filter(Category == cat)
      if (nrow(cat_results) == 0) next
      
      for (t in tier_values) {
        tier_results <- cat_results %>%
          filter(Tier == t) %>%
          arrange(p.value) %>%
          head(15)
        
        if (nrow(tier_results) == 0) next
        
        p <- ggplot(tier_results, aes(x = reorder(Term, -log10(p.value)), y = -log10(p.value))) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          labs(title = paste0("Top GO Terms: ", cat, " - Top ", t, " Genes"),
               x = "GO Term", y = "-log10(p-value)") +
          theme_minimal() +
          theme(plot.title = element_text(size = 14, face = "bold"),
                axis.text.y = element_text(size = 9),
                axis.title = element_text(size = 12),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank())
        
        fn_png <- paste0(opt$output_prefix, "_", cat, "_top", t, ".png")
        fn_svg <- paste0(opt$output_prefix, "_", cat, "_top", t, ".svg")
        ggsave(fn_png, p, width = 10, height = 6, dpi = 300)
        ggsave(fn_svg, p, width = 10, height = 6)
        plot_files[[paste0(cat, "_", t)]] <- list(png = fn_png, svg = fn_svg)
      }
    }
  }
}

# ---- HTML Report ----
html_file <- paste0(opt$output_prefix, "_tpm_enrichment.html")

# Helper: encode a PNG file as base64 data URI
img_to_base64 <- function(filepath) {
  if (!file.exists(filepath)) return(NULL)
  tryCatch({
    raw <- readBin(filepath, "raw", file.info(filepath)$size)
    paste0("data:image/png;base64,", base64enc::base64encode(raw))
  }, error = function(e) {
    # Fallback without base64enc package: use system base64
    tryCatch({
      b64 <- system2("base64", args = c("-w", "0", filepath), stdout = TRUE)
      paste0("data:image/png;base64,", paste(b64, collapse = ""))
    }, error = function(e2) NULL)
  })
}

# Format p-value for display
fmt_pval <- function(p) {
  if (is.na(p)) return("N/A")
  if (p < 0.001) sprintf("%.2e", p) else sprintf("%.4f", p)
}

if (nrow(all_results) > 0) {
  significant <- all_results %>% filter(p.value <= opt$pvalue)
  practically_significant <- all_results %>% filter(Practically_Significant == TRUE)
  n_large <- sum(all_results$Effect_Size_Category == "Large", na.rm = TRUE)
  n_medium <- sum(all_results$Effect_Size_Category == "Medium", na.rm = TRUE)

  # ---- CSS ----
  css <- '
:root {
    --primary: #2d5a87;
    --primary-light: #4a8bc2;
    --primary-dark: #1d3d5c;
    --accent: #e8913a;
    --bg: #f0f2f5;
    --card: #ffffff;
    --text: #2c3e50;
    --text-light: #5a6c7d;
    --border: #e0e4e8;
    --enriched: #27ae60;
    --depleted: #c0392b;
    --sig-bg: #fef9e7;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: "Segoe UI", -apple-system, BlinkMacSystemFont, sans-serif;
    line-height: 1.6; color: var(--text);
    max-width: 1300px; margin: 0 auto; padding: 24px;
    background: var(--bg);
}
.header {
    background: linear-gradient(135deg, var(--primary), var(--primary-dark));
    color: white; padding: 36px 40px; border-radius: 12px;
    margin-bottom: 28px; position: relative; overflow: hidden;
}
.header::after {
    content: ""; position: absolute; top: -40%; right: -10%;
    width: 300px; height: 300px; border-radius: 50%;
    background: rgba(255,255,255,0.06);
}
.header h1 { font-size: 1.9em; font-weight: 700; margin-bottom: 6px; position: relative; z-index: 1; }
.header .subtitle { opacity: 0.85; font-size: 1.05em; position: relative; z-index: 1; }
.header .meta { opacity: 0.65; font-size: 0.85em; margin-top: 10px; position: relative; z-index: 1; }

.section {
    background: var(--card); padding: 28px 32px;
    margin-bottom: 24px; border-radius: 10px;
    box-shadow: 0 1px 6px rgba(0,0,0,0.06);
}
.section h2 {
    color: var(--primary); font-size: 1.35em; font-weight: 700;
    border-bottom: 2px solid var(--primary-light);
    padding-bottom: 10px; margin-bottom: 20px;
}
.section h3 { color: var(--text); font-size: 1.1em; margin: 24px 0 12px 0; }
.section p { margin-bottom: 12px; }

.stats-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
    gap: 16px; margin: 16px 0 8px 0;
}
.stat-card {
    background: var(--bg); padding: 18px 16px;
    border-radius: 8px; text-align: center;
    border-left: 4px solid var(--primary);
    transition: transform 0.15s;
}
.stat-card:hover { transform: translateY(-2px); }
.stat-card.accent { border-left-color: var(--accent); }
.stat-card .label {
    font-size: 0.78em; text-transform: uppercase;
    letter-spacing: 0.5px; color: var(--text-light); font-weight: 600;
}
.stat-card .value {
    font-size: 1.8em; font-weight: 700; color: var(--text);
    margin: 6px 0 2px 0;
}
.stat-card .detail { font-size: 0.82em; color: var(--text-light); }

table {
    width: 100%; border-collapse: collapse;
    margin: 16px 0; font-size: 0.88em;
}
th, td { padding: 10px 14px; text-align: left; border-bottom: 1px solid var(--border); }
th {
    background: var(--primary); color: white;
    font-weight: 600; font-size: 0.85em;
    text-transform: uppercase; letter-spacing: 0.3px;
}
tr:hover { background: #f8fafb; }
.large-effect { background-color: #d5f5e3; }
.medium-effect { background-color: var(--sig-bg); }
.small-effect { background-color: #fdebd0; }
.stable { color: var(--enriched); font-weight: 700; }
.unstable { color: var(--text-light); }

.fig-container { text-align: center; margin: 24px 0; }
.fig-container img {
    max-width: 100%; height: auto; border-radius: 8px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08);
}
.fig-caption {
    color: var(--text-light); font-size: 0.88em;
    font-style: italic; margin-top: 8px;
}
.two-col {
    display: grid; grid-template-columns: 1fr 1fr; gap: 24px;
}
@media (max-width: 850px) { .two-col { grid-template-columns: 1fr; } }

.methods-box {
    background: #f8fafb; padding: 20px 24px;
    border-radius: 8px; border-left: 4px solid var(--accent);
    font-size: 0.9em; line-height: 1.7; color: var(--text-light);
}
.methods-box p { margin-bottom: 10px; }
.methods-box strong { color: var(--text); }
.footer {
    text-align: center; padding: 20px; color: var(--text-light);
    font-size: 0.82em;
}
'

  # ---- Build stat cards ----
  stat_cards <- paste0('
    <div class="stats-grid">
        <div class="stat-card">
            <div class="label">Genes Analysed</div>
            <div class="value">', format(total_genes, big.mark = ","), '</div>
        </div>
        <div class="stat-card">
            <div class="label">GO Terms Tested</div>
            <div class="value">', format(nrow(all_results), big.mark = ","), '</div>
        </div>
        <div class="stat-card accent">
            <div class="label">Statistically Significant</div>
            <div class="value">', nrow(significant), '</div>
            <div class="detail">p &le; ', opt$pvalue, '</div>
        </div>
        <div class="stat-card accent">
            <div class="label">Practically Significant</div>
            <div class="value">', nrow(practically_significant), '</div>
            <div class="detail">FC &ge; ', opt$min_fold_enrichment, ', CI &gt; 1.0</div>
        </div>
        <div class="stat-card">
            <div class="label">Large Effect (FC &ge; 3)</div>
            <div class="value">', n_large, '</div>
        </div>
        <div class="stat-card">
            <div class="label">Bootstrap Iterations</div>
            <div class="value">', opt$bootstrap_n, '</div>
        </div>
    </div>')

  # ---- Build practically significant table ----
  prac_sig_table <- ""
  if (nrow(practically_significant) > 0) {
    rows_html <- ""
    for (i in seq_len(min(30, nrow(practically_significant)))) {
      row <- practically_significant[i, ]
      stab_class <- if (!is.na(row$Stability_Score) && row$Stability_Score > 0.7) "stable" else "unstable"
      eff_class <- tolower(paste0(row$Effect_Size_Category, "-effect"))
      stab_str <- if (!is.na(row$Stability_Score)) sprintf("%.3f", row$Stability_Score) else "N/A"

      rows_html <- paste0(rows_html, '
            <tr class="', eff_class, '">
                <td><strong>', row$Term, '</strong><br><small style="color:var(--text-light)">', row$GO.ID, '</small></td>
                <td>', row$Category, '</td>
                <td>Top ', row$Tier, '</td>
                <td>', sprintf("%.2f", row$Fold_Enrichment), 'x</td>
                <td>', sprintf("%.2f", row$FC_CI_Lower), ' &ndash; ', sprintf("%.2f", row$FC_CI_Upper), '</td>
                <td>', row$Effect_Size_Category, '</td>
                <td class="', stab_class, '">', stab_str, '</td>
                <td>', fmt_pval(row$p.value), '</td>
                <td>', fmt_pval(row$FDR), '</td>
            </tr>')
    }

    remaining <- nrow(practically_significant) - 30
    remaining_note <- if (remaining > 0) paste0(
      '<p style="font-size:0.85em;color:var(--text-light);margin-top:8px">',
      'Showing top 30 of ', nrow(practically_significant), ' results. ',
      'See TSV output for complete data.</p>') else ""

    prac_sig_table <- paste0('
    <div class="section">
        <h2>Practically Significant Results</h2>
        <p>GO terms with fold enrichment &ge; ', opt$min_fold_enrichment,
        ' and 95% CI lower bound &gt; 1.0, ranked by p-value:</p>
        <table>
            <thead>
            <tr><th>GO Term</th><th>Category</th><th>Tier</th><th>Fold Change</th>
                <th>95% CI</th><th>Effect Size</th><th>Stability</th><th>p-value</th><th>FDR</th></tr>
            </thead>
            <tbody>', rows_html, '
            </tbody>
        </table>
        ', remaining_note, '
        <p style="font-size:0.82em;color:var(--text-light);">
            Row colours: <span style="background:#d5f5e3;padding:2px 8px;border-radius:3px">Large effect (&ge;3x)</span>
            <span style="background:#fef9e7;padding:2px 8px;border-radius:3px">Medium (&ge;2x)</span>
            <span style="background:#fdebd0;padding:2px 8px;border-radius:3px">Small (&ge;', opt$min_fold_enrichment, 'x)</span>
        </p>
    </div>')
  }

  # ---- Build visualisation section with base64-embedded images ----
  viz_html <- ""
  fig_num <- 1

  if (length(plot_files) > 0) {
    viz_html <- '
    <div class="section">
        <h2>Enrichment Visualisations</h2>
        <p>Top enriched GO terms per category and expression tier, ranked by statistical significance.</p>'

    for (cat in categories) {
      cat_has_plots <- FALSE
      cat_html <- ""

      for (t in tier_values) {
        plot_key <- paste0(cat, "_", t)
        if (!plot_key %in% names(plot_files)) next

        png_path <- plot_files[[plot_key]]$png
        b64 <- img_to_base64(png_path)

        if (!is.null(b64)) {
          cat_has_plots <- TRUE
          cat_label <- switch(cat,
            "BP" = "Biological Process",
            "MF" = "Molecular Function",
            "CC" = "Cellular Component",
            cat)

          cat_html <- paste0(cat_html, '
        <div class="fig-container">
            <img src="', b64, '" alt="', cat_label, ' Top ', t, '">
            <p class="fig-caption">Figure ', fig_num, ': ', cat_label, ' &mdash; Top ', t, ' genes</p>
        </div>')
          fig_num <- fig_num + 1
        }
      }

      if (cat_has_plots) {
        cat_label <- switch(cat,
          "BP" = "Biological Process",
          "MF" = "Molecular Function",
          "CC" = "Cellular Component",
          cat)
        viz_html <- paste0(viz_html, '
        <h3>', cat_label, '</h3>
        <div class="two-col">', cat_html, '</div>')
      }
    }

    viz_html <- paste0(viz_html, '
    </div>')
  }

  # ---- Build methods section ----
  methods_html <- paste0('
    <div class="section">
        <h2>Methods</h2>
        <div class="methods-box">
            <p><strong>Enrichment algorithm:</strong> GO enrichment was tested using topGO with the
            weight01 algorithm, which accounts for the hierarchical structure of the Gene Ontology
            DAG by decorrelating p-values across parent/child terms
            (Alexa &amp; Rahnenf&uuml;hrer, 2009).</p>

            <p><strong>Expression tiers:</strong> Genes were ranked by TPM (Transcripts Per Million)
            from Salmon quantification and stratified into tiers (top N genes). Each tier was tested
            independently against the full expressed gene set as background.</p>

            <p><strong>Fold enrichment:</strong> Ratio of observed to expected genes annotated with
            each GO term in the tier. Expected counts are derived from the proportion of annotated
            genes in the full dataset.</p>

            <p><strong>Confidence intervals:</strong> 95% exact binomial (Clopper-Pearson) intervals
            on the observed proportion, transformed to the fold enrichment scale.</p>

            <p><strong>FDR correction:</strong> Benjamini-Hochberg correction applied
            <strong>globally</strong> across all ', nrow(all_results), ' tests
            (tier &times; category combinations). Because weight01 p-values are already
            decorrelated for the GO hierarchy, BH correction may be somewhat conservative.</p>

            <p><strong>Bootstrap stability:</strong> ', opt$bootstrap_n, ' iterations of resampling
            genes with replacement and recomputing tier membership. The stability score (1/(1+CV))
            maps coefficient of variation to [0,1]; values &gt; 0.7 indicate stable fold enrichment
            estimates.</p>

            <p><strong>Practical significance:</strong> A GO term is practically significant when
            fold enrichment &ge; ', opt$min_fold_enrichment, ' AND the lower 95% CI bound exceeds
            1.0, ensuring the enrichment is both large enough to be biologically meaningful and
            statistically robust.</p>
        </div>
    </div>

    <div class="section">
        <h2>Interpretation Guide</h2>
        <div class="methods-box">
            <p><strong>Fold enrichment &gt; 1:</strong> GO term is over-represented in the expression
            tier. A value of 2.0 means twice as many annotated genes as expected by chance.</p>

            <p><strong>Effect size categories:</strong> Large (&ge;3x), Medium (&ge;2x),
            Small (&ge;', opt$min_fold_enrichment, 'x). These are heuristic categories for
            convenience, not established statistical thresholds.</p>

            <p><strong>Stability score:</strong> Values &gt; 0.7 (shown in green) indicate that the
            fold enrichment estimate is robust to resampling. Low stability may indicate that a few
            genes drive the enrichment signal.</p>

            <p><strong>Limitations:</strong> This is a single-sample analysis. Results characterise
            expression-dependent functional enrichment patterns but require independent validation.
            The tier-based approach assumes that highly expressed genes reflect the dominant biological
            processes in the sample. Annotation bias may affect results &mdash; well-studied pathways
            may appear more enriched due to more complete GO annotation.</p>
        </div>
    </div>')

  # ---- Assemble full report ----
  html_content <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SCEPTR GO Enrichment Report</title>
    <style>', css, '</style>
</head>
<body>

<div class="header">
    <h1>SCEPTR GO Enrichment Report</h1>
    <div class="subtitle">Expression-weighted Gene Ontology enrichment analysis</div>
    <div class="meta">Generated ', format(Sys.time(), "%Y-%m-%d %H:%M"), ' &middot; SCEPTR TPMGOEnrichment</div>
</div>

<div class="section">
    <h2>Analysis Summary</h2>
    ', stat_cards, '
</div>

', prac_sig_table, '

', viz_html, '

', methods_html, '

<div class="footer">
    SCEPTR v1.0.0 &middot; TPMGOEnrichment &middot; topGO weight01 + global BH FDR
    &middot; ', format(Sys.time(), "%Y-%m-%d"), '
</div>

</body>
</html>')

  writeLines(html_content, html_file)
  message(paste("Generated HTML report:", html_file))
} else {
  writeLines(paste0('<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>SCEPTR GO Enrichment</title>
<style>body{font-family:"Segoe UI",sans-serif;max-width:800px;margin:40px auto;padding:24px;color:#2c3e50;}
.card{background:#fff;border-radius:10px;padding:28px;box-shadow:0 1px 6px rgba(0,0,0,0.06);}
h1{color:#2d5a87;}</style></head>
<body><div class="card"><h1>SCEPTR GO Enrichment Report</h1>
<p>No GO enrichment results were detected. This may indicate:</p>
<ul><li>Insufficient GO annotations for expressed genes</li>
<li>No significant enrichment at the specified thresholds</li>
<li>Expression tiers too small for meaningful analysis</li></ul>
<p style="color:#5a6c7d;margin-top:16px">Generated ', format(Sys.time(), "%Y-%m-%d %H:%M"),
' &middot; SCEPTR v1.0.0</p></div></body></html>'), html_file)
}

message("\nTPM-based GO enrichment analysis completed.")
message(paste("Results:", output_file))
message(paste("Report:", html_file))
message(paste("Plots:", length(plot_files), "sets"))
