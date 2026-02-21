#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Use projectDir directly - it's always correct
sceptr_base = projectDir

log.info "Base directory: ${sceptr_base}"
log.info "Script file: ${workflow.scriptFile}"
log.info "Launch directory: ${workflow.launchDir}"

/*
 * SCEPTR: Single-sample Characterisation of Expression Profiles in Transcriptomics
 * Expression-Weighted Functional Profiling Pipeline
 * Version 1.0.0
 * 
 * REQUIRES: RNA-seq reads + transcriptome assembly for expression quantification
 */

// All parameter defaults are defined in nextflow.config (single source of truth).
// No parameter reassignment here - config values are used directly.

// Function to resolve category files based on user selection
def resolveCategoryFiles() {
    // Host path for file existence checks (runs on host, not in container)
    def categories_host = "${projectDir}/modules/explot/categories"
    // Container path for process execution
    def categories_container = "${params.modules_dir}/explot/categories"
    
    log.info "Starting category file resolution..."
    log.info "Categories host directory: ${categories_host}"
    log.info "Categories container directory: ${categories_container}"
    log.info "Selected category set: ${params.category_set}"
    
    // If custom files provided, use them
    if (params.category_set == "custom") {
        if (!params.custom_functional_categories || !params.custom_cellular_categories) {
            error """
            When using --category_set custom, you must provide both:
              --custom_functional_categories /path/to/functional.json
              --custom_cellular_categories /path/to/cellular.json
            """
        }
        
        if (!file(params.custom_functional_categories).exists()) {
            error "Custom functional categories file not found: ${params.custom_functional_categories}"
        }
        if (!file(params.custom_cellular_categories).exists()) {
            error "Custom cellular categories file not found: ${params.custom_cellular_categories}"
        }
        
        params.explot_functional_categories = params.custom_functional_categories
        params.explot_cellular_categories = params.custom_cellular_categories
        log.info "Using custom category files:"
        log.info "  Functional: ${params.explot_functional_categories}"
        log.info "  Cellular: ${params.explot_cellular_categories}"
        return
    }
    
    // Validate category set
    def available_sets = params.available_category_sets.split(',')
    if (!available_sets.contains(params.category_set)) {
        error """
        Invalid category set: ${params.category_set}
        Available sets: ${available_sets.join(', ')}
        
        Category Set Descriptions:
          general               - Universal categories, works for any organism
          parasite_protozoan    - Single-celled parasites (Toxoplasma, Plasmodium, etc.)
          parasite_metazoan     - Multicellular parasites (helminths, arthropods)
          protist_dinoflagellate - Dinoflagellates and related protists (blooms, toxins)
          model_organism        - Well-studied species (mouse, human, fly, worm, yeast)
          plant                 - Plant-specific processes (photosynthesis, etc.)
          bacteria              - Prokaryotic systems (virulence, cell wall, motility, AMR)
          bacteria_gram_negative - Gram-negative bacteria (LPS, outer membrane, T3SS/T6SS, porins)
          bacteria_gram_positive - Gram-positive bacteria (teichoic acids, sortase, sporulation, competence)
          virus                 - Viral life cycle (replication, entry, immune evasion)
          vertebrate_host       - Host response to infection (immunity, inflammation)
          vertebrate_host_hallmark - High-resolution host response (28 Hallmark-inspired pathways)
          custom                - User-provided JSON files

        Usage examples:
          --category_set general
          --category_set parasite_protozoan
          --category_set bacteria
          --category_set bacteria_gram_negative
          --category_set bacteria_gram_positive
          --category_set virus
          --category_set vertebrate_host
          --category_set vertebrate_host_hallmark
          --category_set custom --custom_functional_categories my_func.json --custom_cellular_categories my_cell.json
        """
    }
    
    // Resolve category files based on set
    // Check existence on host filesystem
    def functional_host = "${categories_host}/functional/${params.category_set}_functional_categories.json"
    def cellular_host = "${categories_host}/cellular/${params.category_set}_cellular_categories.json"
    
    // Container paths for process execution
    def functional_container = "${categories_container}/functional/${params.category_set}_functional_categories.json"
    def cellular_container = "${categories_container}/cellular/${params.category_set}_cellular_categories.json"
    
    // Check if files exist on host
    if (!file(functional_host).exists()) {
        error "Functional categories file not found: ${functional_host}"
    }
    if (!file(cellular_host).exists()) {
        error "Cellular categories file not found: ${cellular_host}"
    }
    
    // Set container paths for processes to use
    params.explot_functional_categories = functional_container
    params.explot_cellular_categories = cellular_container
    
    log.info "Successfully resolved category files:"
    log.info "  Functional: ${functional_container}"
    log.info "  Cellular: ${cellular_container}"
}

// Function to validate required inputs
def validateInputs() {
    def errors = []
    
    log.info "Validating inputs..."
    log.info "Reads pattern: ${params.reads}"
    log.info "Transcripts path: ${params.transcripts}"
    
    // Check if reads pattern is valid
    if (!params.reads) {
        errors.add("ERROR: --reads parameter is required (e.g. --reads 'data/*_{1,2}.fastq.gz')")
    }
    
    // Check transcriptome file
    if (!params.transcripts) {
        errors.add("ERROR: --transcripts parameter is required (e.g. --transcripts assembly.fasta)")
    } else {
        def transcriptsFile = file(params.transcripts)
        if (!transcriptsFile.exists()) {
            errors.add("ERROR: Transcriptome file not found: ${params.transcripts}")
        } else {
            log.info "Transcriptome file found: ${transcriptsFile.name}"
        }
    }
    
    // Check database files
    def contaminantsDb = file(params.contaminants_db)
    if (!contaminantsDb.exists()) {
        errors.add("ERROR: Contaminants database not found: ${params.contaminants_db}")
    }
    
    def uniprotDb = file(params.uniprot_db)
    if (!uniprotDb.exists()) {
        errors.add("ERROR: UniProt database not found: ${params.uniprot_db}")
    }
    
    // Check host transcriptome if provided
    if (params.host_transcriptome && !params.skip_host_filter) {
        def hostFile = file(params.host_transcriptome)
        if (!hostFile.exists()) {
            errors.add("ERROR: Host transcriptome file not found: ${params.host_transcriptome}")
        } else {
            log.info "Host transcriptome file found: ${hostFile.name}"
        }
    }
    
    // Check host proteome if provided
    if (params.host_proteome && !params.skip_host_filter) {
        def hostProteomeFile = file(params.host_proteome)
        if (!hostProteomeFile.exists()) {
            errors.add("ERROR: Host proteome file not found: ${params.host_proteome}")
        } else {
            log.info "Host proteome file found: ${hostProteomeFile.name}"
        }
    }
    
    // Validate that only one of host_transcriptome or host_proteome is provided
    if (params.host_transcriptome && params.host_proteome && !params.skip_host_filter) {
        errors.add("ERROR: Cannot provide both --host_transcriptome and --host_proteome. Please choose one.")
    }
    
    if (errors.size() > 0) {
        errors.each { error -> log.error(error) }
        error """
        SCEPTR requires both RNA-seq reads and a transcriptome assembly for Expression-Weighted Functional Profiling.
        
        Required inputs:
        - RNA-seq reads (--reads): paired-end FASTQ files with pattern like '*_{1,2}.fastq.gz'
        - Transcriptome assembly (--transcripts): FASTA file
        - Database files must exist and be accessible
        
        Optional (for parasite/pathogen studies):
        - Host transcriptome (--host_transcriptome): FASTA nucleotide file for host sequence removal
          (Will be automatically translated to proteins using TransDecoder)
        - Host proteome (--host_proteome): FASTA protein file for host sequence removal
          (Use if you already have translated host proteins - skips translation step)
        
        Note: Provide either --host_transcriptome OR --host_proteome, not both.
        
        Example usage:
          # Standard analysis
          nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --transcripts assembly.fasta
          
          # With host filtering (transcriptome - will be translated automatically):
          nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --transcripts parasite.fasta --host_transcriptome host.fasta
          
          # With host filtering (pre-translated proteome - faster):
          nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --transcripts parasite.fasta --host_proteome host_proteins.fasta
        """
    }
    
    return true
}

// Create essential output directories
def prepareOutputDirs(outputPath) {
    try {
        def outDir = new File(outputPath)
        if (!outDir.exists()) {
            log.info "Creating output directory: ${outputPath}"
            if (!outDir.mkdirs()) {
                log.warn "Warning: Could not create output directory: ${outputPath}"
            }
        }
        
        def pipelineInfoDir = new File("${outputPath}/pipeline_info")
        if (!pipelineInfoDir.exists()) {
            pipelineInfoDir.mkdirs()
        }
    } catch (Exception e) {
        log.warn "Could not create output directories: ${e.message}"
        log.info "Directories will be created by Nextflow processes"
    }
}

// Verify required parameters and file existence
def checkParameters() {
    log.info "Verifying input parameters and files..."
    log.info "Using profile: ${workflow.profile}"
    
    log.info """
    ==============================================
    SCEPTR Pipeline v1.0.0: Expression-Weighted Functional Profiling
    ==============================================
    
    Execution Mode: ${workflow.profile ?: 'local'}
    
    Input Data:
      Reads pattern: ${params.reads}
      Transcripts: ${params.transcripts}
      Output directory: ${params.outdir}
      
    Databases:
      Contaminants DB: ${params.contaminants_db}
      UniProt DB: ${params.uniprot_db}
    """
    
    // Log host filtering status
    if ((params.host_transcriptome || params.host_proteome) && !params.skip_host_filter) {
        def hostSource = params.host_proteome ? "proteome (pre-translated proteins)" : "transcriptome (will be translated)"
        def hostFile = params.host_proteome ?: params.host_transcriptome
        log.info """
      
    Host Filtering: ENABLED (Parasite/Pathogen Mode)
      Host ${hostSource}: ${hostFile}
      Identity threshold: ${params.host_identity}%
      Coverage threshold: ${params.host_coverage}%
      E-value threshold: ${params.host_evalue}
      
      ${params.host_proteome ? '✓ Using pre-translated host proteome (faster)' : '⚠️  Host transcriptome will be translated to proteins using TransDecoder'}
      ⚠️  Host sequences will be removed before general contamination filtering
        """
    } else if ((params.host_transcriptome || params.host_proteome) && params.skip_host_filter) {
        log.info """
      
    Host Filtering: DISABLED (--skip_host_filter flag)
      Note: Host file provided but filtering skipped
        """
    } else {
        log.info """
      
    Host Filtering: DISABLED
      Use --host_transcriptome or --host_proteome for parasite/pathogen studies
        """
    }
    
    log.info """
      
    Category Selection:
      Category set: ${params.category_set}
      Functional categories: ${params.explot_functional_categories}
      Cellular categories: ${params.explot_cellular_categories}
      
    Expression Analysis Parameters:
      TPM threshold: ${params.tpm_threshold}
      Expression tiers: ${params.expression_tiers}
      GO P-value cutoff: ${params.go_pvalue_cutoff}
      Bootstrap iterations: ${params.bootstrap_n}
      Skip ExPlot: ${params.skip_explot}
      
    This pipeline performs Expression-Weighted Functional Profiling,
    analyzing how functional categories are distributed across
    expression tiers to identify biological priorities.
    ==============================================
    """
}

// Include process modules
include { FastQC; MultiQC } from "${sceptr_base}/modules/qc.nf"
include { SalmonIndex; SalmonQuant } from "${sceptr_base}/modules/salmon.nf"
include { TransDecoder; DirectTranslate } from "${sceptr_base}/modules/transdecoder.nf"
include { TranslateHostTranscriptome; BuildHostDatabase; FilterHostSequences; DiamondBlastContaminants; VisualiseContaminants } from "${sceptr_base}/modules/contamination.nf"
include { DiamondBlastUniProt; ExtractUniProtIDs; AnnotateProteins; AnnotateGOTerms; GenerateAnnotationReport } from "${sceptr_base}/modules/annotation.nf"
include { PrepareGOEnrichmentData; MergeExpressionAnnotations; TPMGOEnrichment } from "${sceptr_base}/modules/enrichment/enrichment.nf"
include { ExPlot } from "${sceptr_base}/modules/explot/expression_profiling.nf"
include { Landscape } from "${sceptr_base}/modules/landscape/landscape.nf"
include { CompareConditions } from "${sceptr_base}/modules/comparison.nf"

// MAIN WORKFLOW - Expression-Weighted Functional Profiling
workflow {
    // STEP 1: Resolve category files FIRST (before any validation)
    resolveCategoryFiles()
    
    // STEP 2: Validate inputs
    validateInputs()
    
    // STEP 3: Make sure the output directory exists
    prepareOutputDirs(params.outdir)
    
    // STEP 4: Verify all key files exist and log pipeline info
    checkParameters()
    
    // Create channels for transcripts and databases
    transcripts_ch = Channel.fromPath(params.transcripts)
    contaminants_db_ch = Channel.fromPath(params.contaminants_db)
    uniprot_db_ch = Channel.fromPath(params.uniprot_db)
    
    log.info "Running SCEPTR Expression-Weighted Functional Profiling Pipeline"
    
    // DYNAMIC READ DISCOVERY: Create channels from user-provided reads pattern
    log.info "Discovering read files with pattern: ${params.reads}"
    
    if (params.single_end) {
        // Single-end mode: each file is one sample
        log.info "Mode: Single-end reads"
        read_pairs_ch = Channel
            .fromPath(params.reads, checkIfExists: true)
            .ifEmpty { error "No read files found matching pattern: ${params.reads}" }
            .map { file ->
                def name = file.name.toString()
                def sample = name.replaceAll(/(_R1_001|_R1|_1)?\.f(ast)?q(\.gz)?$/, '')
                log.info "Found single-end reads for sample: ${sample}"
                log.info "  File: ${file.name}"
                return tuple(sample, file)
            }

        // FastQC on individual files
        fastqc_reads_ch = read_pairs_ch
            .map { sample_id, reads -> tuple(sample_id, reads) }
    } else {
        // Paired-end mode (default)
        log.info "Mode: Paired-end reads"
        read_pairs_ch = Channel
            .fromFilePairs(params.reads, size: 2, checkIfExists: true)
            .ifEmpty { error "No read pairs found matching pattern: ${params.reads}" }
            .map { sample_id, files ->
                log.info "Found read pair for sample: ${sample_id}"
                log.info "  Read 1: ${files[0].name}"
                log.info "  Read 2: ${files[1].name}"
                return tuple(sample_id, files)
            }

        // Create FastQC channel by splitting read pairs into individual files
        fastqc_reads_ch = read_pairs_ch
            .flatMap { sample_id, files ->
                [
                    tuple("${sample_id}_R1", files[0]),
                    tuple("${sample_id}_R2", files[1])
                ]
            }
    }
    
    // STEP 1: Quality Control
    log.info "Step 1: Quality Control Analysis"
    fastqc_results = FastQC(fastqc_reads_ch)
    collected_fastqc = fastqc_results.html.mix(fastqc_results.zip).collect()
    multiqc_results = MultiQC(collected_fastqc)
    
    // STEP 2: Expression Quantification
    log.info "Step 2: Expression Quantification with Salmon"
    salmon_index = SalmonIndex(transcripts_ch)
    salmon_quant = SalmonQuant(read_pairs_ch, salmon_index)
    
    // STEP 3: ORF Prediction / CDS Translation
    if (params.proteome) {
        // User provided pre-translated protein sequences - skip translation entirely
        log.info "Step 3: Using pre-translated proteome (--proteome provided, skipping ORF prediction)"
        proteome_ch = Channel.fromPath(params.proteome, checkIfExists: true)
        transdecoder_results = [ predicted_proteome: proteome_ch ]
    } else {
        // Auto-detect: bacteria and virus category sets typically use CDS inputs
        def use_direct_translate = params.skip_transdecoder || params.category_set in ['bacteria', 'bacteria_gram_negative', 'bacteria_gram_positive', 'virus', 'vertebrate_host', 'vertebrate_host_hallmark', 'model_organism']
        if (use_direct_translate) {
            if (!params.skip_transdecoder && params.category_set in ['bacteria', 'bacteria_gram_negative', 'bacteria_gram_positive', 'virus']) {
                log.info "Step 3: Direct CDS Translation (auto-detected from '${params.category_set}' category set)"
                log.info "        Override with --skip_transdecoder false if using a de novo transcriptome assembly"
            } else {
                log.info "Step 3: Direct CDS Translation (TransDecoder skipped)"
            }
            transdecoder_results = DirectTranslate(transcripts_ch)
        } else {
            log.info "Step 3: ORF Prediction with TransDecoder"
            transdecoder_results = TransDecoder(transcripts_ch)
        }
    }
    
    // STEP 4: Contamination Detection (with optional host filtering)
    log.info "Step 4: Contamination Detection and Filtering"
    
    // CONDITIONAL: Host filtering for parasite/pathogen studies
    if ((params.host_transcriptome || params.host_proteome) && !params.skip_host_filter) {
        log.info "  4a: Host Sequence Filtering (Parasite/Pathogen Mode)"
        
        // Determine which host input to use
        if (params.host_proteome) {
            // User provided pre-translated proteome - use directly
            log.info "  4a.1: Using pre-translated host proteome (skipping translation)"
            host_proteome_ch = Channel.fromPath(params.host_proteome, checkIfExists: true)
            
            // Build host database directly from provided proteome
            log.info "  4a.2: Building DIAMOND database from host proteins"
            host_db = BuildHostDatabase(host_proteome_ch)
            
        } else {
            // User provided transcriptome - translate first
            log.info "  4a.1: Translating host transcriptome to proteins with TransDecoder"
            host_transcriptome_ch = Channel.fromPath(params.host_transcriptome, checkIfExists: true)
            host_translation = TranslateHostTranscriptome(host_transcriptome_ch)
            
            // Build host database from translated proteins
            log.info "  4a.2: Building DIAMOND database from host proteins"
            host_db = BuildHostDatabase(host_translation.host_proteome)
        }
        
        // Filter host sequences
        log.info "  4a.3: Filtering host sequences from organism proteome"
        host_filter_results = FilterHostSequences(
            transdecoder_results.predicted_proteome,
            host_db.host_db
        )
        
        // Display host filtering summary
        host_filter_results.report.view { report_file ->
            def report_content = file(report_file).text
            return "\n${report_content}\n"
        }
        
        // Use host-filtered proteome for next stage
        proteome_for_decon = host_filter_results.filtered_proteome
        
        log.info "  4b: General Contamination Detection (after host removal)"
    } else {
        log.info "  Skipping host filtering (not a parasite/pathogen analysis)"
        proteome_for_decon = transdecoder_results.predicted_proteome
    }
    
    // General contamination filtering
    // Auto-skip for bacteria/virus - the contaminant database contains bacterial proteins
    // which would incorrectly flag the target organism's own proteins
    def skip_contam = params.skip_contamination ?: (params.category_set in ['bacteria', 'bacteria_gram_negative', 'bacteria_gram_positive', 'virus', 'vertebrate_host', 'vertebrate_host_hallmark', 'model_organism'])
    if (skip_contam) {
        log.info "  Skipping contaminant filtering (${params.category_set} category - contaminant DB contains bacterial/viral proteins)"
        proteome_for_annotation = proteome_for_decon
    } else {
        decon_results = DiamondBlastContaminants(proteome_for_decon, contaminants_db_ch)
        vis_results = VisualiseContaminants(decon_results.details, decon_results.params_json)
        proteome_for_annotation = decon_results.filtered_proteome
    }
    
    // STEP 5: Functional Annotation
    log.info "Step 5: Functional Annotation with UniProt and GO"
    blast_results = DiamondBlastUniProt(proteome_for_annotation, uniprot_db_ch)
    uniprot_mapping = ExtractUniProtIDs(blast_results.blast_results)
    protein_annotations = AnnotateProteins(uniprot_mapping.uniprot_mapping)
    final_annotations = AnnotateGOTerms(protein_annotations.protein_annotations)
    annotation_report = GenerateAnnotationReport(final_annotations.final_annotations)
    
    // STEP 6: Expression-Weighted Analysis
    log.info "Step 6: Expression-Weighted Functional Analysis"
    go_data_ch = PrepareGOEnrichmentData(final_annotations.final_annotations)
    
    salmon_quant_dirs = salmon_quant.map { sample, quant_file -> 
        def quant_dir = quant_file.getParent()
        log.info "Salmon Quant Directory for ${sample}: ${quant_dir}"
        return quant_dir
    }

    // Merge expression data with annotations
    expression_data_ch = MergeExpressionAnnotations(
        final_annotations.final_annotations, 
        salmon_quant_dirs.collect()
    )
    
    // GO Enrichment Analysis with expression tiers
    go_data_file = go_data_ch.go_terms.first()
    expression_data_file = expression_data_ch.merged_annotations.first()
    
    enrichment_results = TPMGOEnrichment(go_data_file, expression_data_file)
    
    // STEP 7: Expression Profiling (ExPlot)
    if (!params.skip_explot) {
        log.info "Step 7: Expression-Weighted Functional Profiling (ExPlot)"
        try {
            ExPlot(expression_data_ch.merged_annotations.first())
            log.info "Expression-Weighted Functional Profiling completed successfully"
        } catch (Exception e) {
            log.warn "ExPlot step failed: ${e.message}. Continuing with pipeline..."
        }
    } else {
        log.info "Step 7: Skipping Expression Profiling (ExPlot) as requested"
    }

    // STEP 8: Transcriptome Landscape Analysis
    if (!params.skip_landscape) {
        log.info "Step 8: Transcriptome Landscape Analysis"
        try {
            Landscape(expression_data_ch.merged_annotations.first())
            log.info "Transcriptome Landscape analysis completed successfully"
        } catch (Exception e) {
            log.warn "Landscape step failed: ${e.message}. Continuing with pipeline..."
        }
    } else {
        log.info "Step 8: Skipping Transcriptome Landscape as requested"
    }
    
    // Output summary
    if (!skip_contam) {
        decon_results.report.view { report_file ->
            def report_content = file(report_file).text
            return "\nSCEPTR Contamination Detection Results:\n${report_content}\n"
        }
    }
    
    log.info "SCEPTR Expression-Weighted Functional Profiling is complete."
}

// Entry point for running just ExPlot
workflow ExPlotEntry {
    // Resolve category files first
    resolveCategoryFiles()
    
    prepareOutputDirs(params.outdir ?: "results")
    
    Channel
        .fromPath(params.integrated_results ?: "${params.outdir ?: 'results'}/integrated_data/integrated_annotations_expression.tsv")
        .ifEmpty { error "No integrated results file found. Run the full pipeline first to generate expression data." }
        .set { integrated_results_ch }
    
    integrated_results_ch.view { file ->
        log.info "Found integrated results file: ${file}"
        return file
    }
        
    ExPlot(integrated_results_ch)
    
    log.info "SCEPTR ExPlot module complete."
}

// SUBWORKFLOW: Run just decontamination
workflow decon {
    prepareOutputDirs(params.outdir ?: "results")

    contaminants_db_ch = Channel.fromPath(params.contaminants_db)
    transcripts_ch = Channel.fromPath(params.transcripts)

    if (params.proteome) {
        proteome_ch = Channel.fromPath(params.proteome, checkIfExists: true)
        transdecoder_results = [ predicted_proteome: proteome_ch ]
    } else {
        def use_direct_translate = params.skip_transdecoder || params.category_set in ['bacteria', 'bacteria_gram_negative', 'bacteria_gram_positive', 'virus', 'vertebrate_host', 'vertebrate_host_hallmark', 'model_organism']
        if (use_direct_translate) {
            transdecoder_results = DirectTranslate(transcripts_ch)
        } else {
            transdecoder_results = TransDecoder(transcripts_ch)
        }
    }
    
    // Conditional host filtering in decon workflow
    if ((params.host_transcriptome || params.host_proteome) && !params.skip_host_filter) {
        log.info "Running host filtering in decontamination workflow"
        
        if (params.host_proteome) {
            // Use pre-translated proteome
            host_proteome_ch = Channel.fromPath(params.host_proteome, checkIfExists: true)
            host_db = BuildHostDatabase(host_proteome_ch)
        } else {
            // Translate transcriptome
            host_transcriptome_ch = Channel.fromPath(params.host_transcriptome, checkIfExists: true)
            host_translation = TranslateHostTranscriptome(host_transcriptome_ch)
            host_db = BuildHostDatabase(host_translation.host_proteome)
        }
        
        host_filter_results = FilterHostSequences(transdecoder_results.predicted_proteome, host_db.host_db)
        proteome_for_decon = host_filter_results.filtered_proteome
    } else {
        proteome_for_decon = transdecoder_results.predicted_proteome
    }
    
    decon_results = DiamondBlastContaminants(proteome_for_decon, contaminants_db_ch)
    
    try {
        vis_results = VisualiseContaminants(decon_results.details, decon_results.params_json)
    } catch (Exception e) {
        log.warn "Visualisation step failed: ${e.message}. Continuing with pipeline..."
    }
    
    decon_results.report.view { report_file ->
        def report_content = file(report_file).text
        return "\nSCEPTR Contamination Filtering Results:\n${report_content}\n"
    }

    log.info "SCEPTR contamination filtering complete!"
}

// SUBWORKFLOW: Run just annotation
workflow annotate {
    prepareOutputDirs(params.outdir ?: "results")

    uniprot_db_ch = Channel.fromPath(params.uniprot_db)
    transcripts_ch = Channel.fromPath(params.transcripts)

    // Translate transcriptome to proteome first (DiamondBlastUniProt expects protein sequences)
    if (params.proteome) {
        proteome_ch = Channel.fromPath(params.proteome, checkIfExists: true)
        transdecoder_results = [ predicted_proteome: proteome_ch ]
    } else {
        def use_direct_translate = params.skip_transdecoder || params.category_set in ['bacteria', 'bacteria_gram_negative', 'bacteria_gram_positive', 'virus', 'vertebrate_host', 'vertebrate_host_hallmark', 'model_organism']
        if (use_direct_translate) {
            transdecoder_results = DirectTranslate(transcripts_ch)
        } else {
            transdecoder_results = TransDecoder(transcripts_ch)
        }
    }
    
    blast_results = DiamondBlastUniProt(transdecoder_results.predicted_proteome, uniprot_db_ch)
    uniprot_mapping = ExtractUniProtIDs(blast_results.blast_results)
    protein_annotations = AnnotateProteins(uniprot_mapping.uniprot_mapping)
    final_annotations = AnnotateGOTerms(protein_annotations.protein_annotations)
    annotation_report = GenerateAnnotationReport(final_annotations.final_annotations)
    
    log.info "SCEPTR functional annotation complete."
}

// SUBWORKFLOW: Cross-sample comparison
workflow compare {
    resolveCategoryFiles()

    // Validate inputs
    if (!params.condition_a || !params.condition_b) {
        error """
        Cross-sample comparison requires two integrated results files.

        Usage:
          nextflow run main.nf -entry compare \\
            --condition_a results_mock/integrated_data/integrated_annotations_expression.tsv \\
            --condition_b results_infected/integrated_data/integrated_annotations_expression.tsv \\
            --label_a "Mock" --label_b "Infected" \\
            --category_set vertebrate_host \\
            -profile docker
        """
    }

    condition_a_ch = Channel.fromPath(params.condition_a, checkIfExists: true)
    condition_b_ch = Channel.fromPath(params.condition_b, checkIfExists: true)

    // Resolve functional categories file
    def func_cats_path = params.explot_functional_categories ?: "NO_CUSTOM_CATEGORIES"
    def func_cats_ch = func_cats_path != "NO_CUSTOM_CATEGORIES"
        ? Channel.fromPath(func_cats_path, checkIfExists: true)
        : Channel.fromPath("${projectDir}/modules/explot/categories/functional/${params.category_set}_functional_categories.json", checkIfExists: true)

    CompareConditions(
        condition_a_ch,
        condition_b_ch,
        params.label_a ?: "Condition_A",
        params.label_b ?: "Condition_B",
        params.category_set,
        func_cats_ch,
        params.expression_tiers,
        params.n_permutations,
        params.output_prefix,
        params.explot_expand_go,
        params.explot_go_expansion_depth,
        params.comparison_seed
    )

    log.info "SCEPTR cross-sample comparison complete."
}

// Handle workflow completion
workflow.onComplete {
    log.info "Pipeline execution completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
    
    if (!workflow.success) {
        log.error "Pipeline execution failed"
    }
}
