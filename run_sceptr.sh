#!/usr/bin/env bash
#
# SCEPTR: Statistical Characterisation of Expression Profiles in Transcriptomes
# Interactive launcher
#

set -euo pipefail

VERSION="1.0.0"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Colours ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
DIM='\033[2m'
NC='\033[0m'

# ── Helper functions ─────────────────────────────────────────────────────────

banner() {
    echo ""
    echo -e "${CYAN}${BOLD}┌─────────────────────────────────────────────────────────────┐${NC}"
    echo -e "${CYAN}${BOLD}│                                                             │${NC}"
    echo -e "${CYAN}${BOLD}│   SCEPTR v${VERSION}                                            │${NC}"
    echo -e "${CYAN}${BOLD}│   Statistical Characterisation of Expression Profiles        │${NC}"
    echo -e "${CYAN}${BOLD}│   in Transcriptomes                                         │${NC}"
    echo -e "${CYAN}${BOLD}│                                                             │${NC}"
    echo -e "${CYAN}${BOLD}└─────────────────────────────────────────────────────────────┘${NC}"
    echo ""
}

info()    { echo -e "${BLUE}ℹ${NC}  $1"; }
success() { echo -e "${GREEN}✓${NC}  $1"; }
warn()    { echo -e "${YELLOW}⚠${NC}  $1"; }
error()   { echo -e "${RED}✗${NC}  $1"; }

prompt() {
    echo -e -n "${BOLD}$1${NC} " >&2
    read -r REPLY
    echo "$REPLY"
}

# Prompt with default value
prompt_default() {
    local message="$1"
    local default="$2"
    echo -e -n "${BOLD}${message}${NC} ${DIM}[${default}]${NC}: " >&2
    read -r REPLY
    if [[ -z "$REPLY" ]]; then
        echo "$default"
    else
        echo "$REPLY"
    fi
}

# Numbered selection menu
select_option() {
    local prompt_msg="$1"
    shift
    local options=("$@")

    echo -e "\n${BOLD}${prompt_msg}${NC}\n" >&2
    for i in "${!options[@]}"; do
        echo -e "  ${CYAN}$((i+1))${NC}  ${options[$i]}" >&2
    done
    echo "" >&2

    while true; do
        echo -e -n "${BOLD}Select [1-${#options[@]}]:${NC} " >&2
        read -r choice
        if [[ "$choice" =~ ^[0-9]+$ ]] && (( choice >= 1 && choice <= ${#options[@]} )); then
            echo "$((choice-1))"
            return
        fi
        echo -e "${RED}  Invalid choice. Please enter a number between 1 and ${#options[@]}.${NC}" >&2
    done
}

# Yes/no prompt
confirm() {
    local message="$1"
    local default="${2:-n}"
    local hint="y/N"
    [[ "$default" == "y" ]] && hint="Y/n"

    echo -e -n "${BOLD}${message}${NC} ${DIM}[${hint}]${NC}: " >&2
    read -r REPLY
    REPLY="${REPLY:-$default}"
    [[ "${REPLY,,}" == "y" || "${REPLY,,}" == "yes" ]]
}

# Auto-detect paired reads in a directory
detect_reads() {
    local dir="$1"
    local patterns=(
        "${dir}/*_1.fastq.gz"
        "${dir}/*_R1.fastq.gz"
        "${dir}/*_1.fq.gz"
        "${dir}/*_R1.fq.gz"
        "${dir}/*_R1_001.fastq.gz"
    )

    for pattern in "${patterns[@]}"; do
        local files=( $pattern )
        if [[ -f "${files[0]:-}" ]]; then
            # Derive the glob pattern
            local sample="${files[0]}"
            if [[ "$sample" == *"_1.fastq.gz" ]]; then
                echo "${dir}/*_{1,2}.fastq.gz"
            elif [[ "$sample" == *"_R1.fastq.gz" ]]; then
                echo "${dir}/*_{R1,R2}.fastq.gz"
            elif [[ "$sample" == *"_1.fq.gz" ]]; then
                echo "${dir}/*_{1,2}.fq.gz"
            elif [[ "$sample" == *"_R1.fq.gz" ]]; then
                echo "${dir}/*_{R1,R2}.fq.gz"
            elif [[ "$sample" == *"_R1_001.fastq.gz" ]]; then
                echo "${dir}/*_{R1,R2}_001.fastq.gz"
            fi
            return 0
        fi
    done
    return 1
}

# Count read files matching a glob
count_read_files() {
    local pattern="$1"
    # Expand the glob
    local expanded
    expanded=$(eval ls $pattern 2>/dev/null | wc -l)
    echo "$expanded"
}

# ── Category set definitions ─────────────────────────────────────────────────

declare -A CATEGORY_NAMES
declare -A CATEGORY_DESCRIPTIONS

CATEGORY_NAMES=(
    [1]="general"
    [2]="human_host"
    [3]="vertebrate_host"
    [4]="cancer"
    [5]="bacteria"
    [6]="bacteria_gram_negative"
    [7]="bacteria_gram_positive"
    [8]="parasite_protozoan"
    [9]="helminth_nematode"
    [10]="helminth_platyhelminth"
    [11]="fungi"
    [12]="plant"
    [13]="protist_dinoflagellate"
    [14]="insect"
    [15]="custom"
)

CATEGORY_DESCRIPTIONS=(
    [1]="General              │ Universal categories, works for any organism"
    [2]="Human host           │ 33 detailed pathway categories (IFN, NFkB, coagulation, etc.)"
    [3]="Vertebrate host      │ 17 broad categories for mouse, fish, birds, etc."
    [4]="Cancer               │ Hallmarks of cancer, EMT, immune evasion, epigenetics"
    [5]="Bacteria (general)   │ Virulence, cell wall, motility, AMR, iron acquisition"
    [6]="Bacteria (Gram-)     │ LPS, outer membrane, T3SS/T6SS, porins, siderophores"
    [7]="Bacteria (Gram+)     │ Teichoic acids, sortase, sporulation, competence"
    [8]="Parasite (protozoan) │ Plasmodium, Toxoplasma, Leishmania, Trypanosoma, etc."
    [9]="Helminth (nematode)  │ Cuticle, dauer, neuromuscular, ES products"
    [10]="Helminth (platyhelminth) │ Tegument, neoblasts, lifecycle stages, egg biology"
    [11]="Fungi                │ Cell wall, secondary metabolism, sporulation, virulence"
    [12]="Plant                │ Photosynthesis, cell wall, hormones, defence"
    [13]="Dinoflagellate       │ Symbiodiniaceae, HABs, coral symbionts, free-living dinos"
    [14]="Insect               │ Cuticle, metamorphosis, chemosensation, detoxification"
    [15]="Custom               │ Provide your own category JSON files"
)

# ── CLI argument parsing ─────────────────────────────────────────────────────

READS=""
TRANSCRIPTS=""
CATEGORY_SET=""
HOST_REF=""
HOST_TYPE=""
OUTDIR=""
PREFIX=""
PROFILE="docker"
TIERS=""
SINGLE_END=false
SKIP_TRANSDECODER=false
RESUME=false
INTERACTIVE=true
EXTRA_ARGS=""
# Comparison mode
COMPARE_MODE=false
CONDITION_A=""
CONDITION_B=""
LABEL_A=""
LABEL_B=""
N_PERMUTATIONS=""

usage() {
    banner
    echo -e "${BOLD}Usage:${NC}"
    echo "  ./run_sceptr.sh                          # Interactive mode"
    echo "  ./run_sceptr.sh [options]                 # Direct mode (full framework)"
    echo "  ./run_sceptr.sh --compare [options]       # Compare two conditions"
    echo ""
    echo -e "${BOLD}Full Framework Options:${NC}"
    echo "  -r, --reads DIR|GLOB       Path to reads directory or glob pattern"
    echo "  -t, --transcripts FILE     Reference transcriptome FASTA"
    echo "  -c, --category SET         Category set (general, human_host, bacteria, etc.)"
    echo "  -s, --single-end           Single-end reads (default: paired-end)"
    echo "      --skip-transdecoder    Skip TransDecoder, directly translate CDS (for bacteria)"
    echo "  -H, --host FILE            Host reference for decontamination"
    echo "      --host-proteome FILE   Pre-translated host proteome (faster)"
    echo "  -o, --outdir DIR           Output directory (default: results)"
    echo "  -p, --prefix STRING        Output file prefix (default: sceptr)"
    echo "      --tiers LIST           Expression tiers (default: 50,100,250,500)"
    echo "      --profile PROFILE      Nextflow profile (default: docker)"
    echo "      --resume               Resume previous run"
    echo "  -h, --help                 Show this help"
    echo ""
    echo -e "${BOLD}Comparison Options:${NC}"
    echo "      --compare              Run cross-sample comparison mode"
    echo "      --condition-a FILE     Condition A integrated_annotations_expression.tsv"
    echo "      --condition-b FILE     Condition B integrated_annotations_expression.tsv"
    echo "      --label-a STRING       Label for condition A (e.g. Mock)"
    echo "      --label-b STRING       Label for condition B (e.g. Infected)"
    echo "      --permutations INT     Number of permutations (default: 10000)"
    echo ""
    echo -e "${BOLD}Examples:${NC}"
    echo "  ./run_sceptr.sh -r data/reads -t assembly.fasta -c bacteria"
    echo "  ./run_sceptr.sh -r 'data/*_{1,2}.fastq.gz' -t assembly.fasta -c parasite_protozoan -H host.fasta"
    echo "  ./run_sceptr.sh -r 'data/*.fastq.gz' -t cds.fasta -c bacteria --single-end"
    echo ""
    echo -e "${BOLD}Comparison Example:${NC}"
    echo "  ./run_sceptr.sh --compare --condition-a results_mock/integrated_data/integrated_annotations_expression.tsv \\"
    echo "    --condition-b results_infected/integrated_data/integrated_annotations_expression.tsv \\"
    echo "    --label-a Mock --label-b Infected -c vertebrate_host"
    echo ""
    echo -e "${BOLD}Category sets:${NC}"
    for i in $(seq 1 15); do
        echo -e "  ${CYAN}${CATEGORY_NAMES[$i]}${NC} - ${CATEGORY_DESCRIPTIONS[$i]#*│ }"
    done
    echo ""
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--reads)        READS="$2"; INTERACTIVE=false; shift 2 ;;
        -t|--transcripts)  TRANSCRIPTS="$2"; INTERACTIVE=false; shift 2 ;;
        -c|--category)     CATEGORY_SET="$2"; INTERACTIVE=false; shift 2 ;;
        -s|--single-end)   SINGLE_END=true; shift ;;
        --skip-transdecoder) SKIP_TRANSDECODER=true; shift ;;
        -H|--host)         HOST_REF="$2"; HOST_TYPE="transcriptome"; shift 2 ;;
        --host-proteome)   HOST_REF="$2"; HOST_TYPE="proteome"; shift 2 ;;
        -o|--outdir)       OUTDIR="$2"; shift 2 ;;
        -p|--prefix)       PREFIX="$2"; shift 2 ;;
        --tiers)           TIERS="$2"; shift 2 ;;
        --profile)         PROFILE="$2"; shift 2 ;;
        --resume)          RESUME=true; shift ;;
        --compare)         COMPARE_MODE=true; INTERACTIVE=false; shift ;;
        --condition-a)     CONDITION_A="$2"; shift 2 ;;
        --condition-b)     CONDITION_B="$2"; shift 2 ;;
        --label-a)         LABEL_A="$2"; shift 2 ;;
        --label-b)         LABEL_B="$2"; shift 2 ;;
        --permutations)    N_PERMUTATIONS="$2"; shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        *)                 EXTRA_ARGS="$EXTRA_ARGS $1"; shift ;;
    esac
done

# ── Interactive mode ─────────────────────────────────────────────────────────

if $INTERACTIVE; then
    banner

    # Mode selection
    mode_idx=$(select_option "What would you like to do?" \
        "Full framework    - process raw reads to enrichment report" \
        "Method only       - run enrichment profiling on an annotated expression table" \
        "Compare conditions - compare two SCEPTR outputs" \
        "Re-run enrichment  - re-analyse existing SCEPTR results with different settings")

    if (( mode_idx == 2 )); then
        # ── Comparison interactive mode ──
        COMPARE_MODE=true
        echo ""
        echo -e "${BOLD}${CYAN}Cross-Sample Comparison${NC}\n"
        echo "  Compare functional enrichment profiles between two conditions"
        echo -e "  ${DIM}Requires two completed SCEPTR runs.${NC}\n"

        # Condition A
        while true; do
            CONDITION_A=$(prompt "  Path to condition A integrated_annotations_expression.tsv:")
            if [[ -f "$CONDITION_A" ]]; then
                n_genes=$(tail -n +2 "$CONDITION_A" | wc -l)
                success "Found condition A: ${n_genes} genes"
                break
            else
                error "File not found: ${CONDITION_A}"
            fi
        done

        # Condition B
        while true; do
            CONDITION_B=$(prompt "  Path to condition B integrated_annotations_expression.tsv:")
            if [[ -f "$CONDITION_B" ]]; then
                n_genes=$(tail -n +2 "$CONDITION_B" | wc -l)
                success "Found condition B: ${n_genes} genes"
                break
            else
                error "File not found: ${CONDITION_B}"
            fi
        done

        # Labels
        echo ""
        LABEL_A=$(prompt_default "  Label for condition A" "Condition_A")
        LABEL_B=$(prompt_default "  Label for condition B" "Condition_B")

        # Category set
        echo ""
        echo -e "${BOLD}${CYAN}Category Set${NC}"
        cat_options=()
        for i in $(seq 1 15); do
            cat_options+=("${CATEGORY_DESCRIPTIONS[$i]}")
        done
        cat_idx=$(select_option "What type of organism are you studying?" "${cat_options[@]}")
        cat_num=$((cat_idx + 1))
        CATEGORY_SET="${CATEGORY_NAMES[$cat_num]}"
        success "Selected: ${CATEGORY_SET}"

        # Handle custom categories
        if [[ "$CATEGORY_SET" == "custom" ]]; then
            echo ""
            while true; do
                CUSTOM_FUNC=$(prompt "  Path to functional categories JSON:")
                if [[ -f "$CUSTOM_FUNC" ]]; then success "Found: ${CUSTOM_FUNC}"; break
                else error "File not found: ${CUSTOM_FUNC}"; fi
            done
        fi

        # Permutations
        echo ""
        N_PERMUTATIONS=$(prompt_default "  Number of permutations" "10000")

        # Output settings
        echo ""
        echo -e "${BOLD}${CYAN}Output Settings${NC}\n"
        OUTDIR=$(prompt_default "  Output directory" "results_comparison")
        PREFIX=$(prompt_default "  Output file prefix" "sceptr_comparison")

    elif (( mode_idx == 1 )); then
        # ── Method only: enrichment profiling on user's own data ──
        echo ""
        echo -e "${BOLD}${CYAN}Enrichment Profiling Only${NC}\n"
        echo "  Run the SCEPTR statistical method on your own annotated expression data."
        echo -e "  ${DIM}Skips all preprocessing - goes straight to enrichment profiling.${NC}"
        echo ""
        echo "  Your input file should be a tab-separated table with at minimum:"
        echo "    - sequence_id: gene/transcript identifiers"
        echo "    - TPM: expression values (used for ranking)"
        echo "    - protein_name: UniProt-style protein descriptions"
        echo -e "  ${DIM}Optional columns: GO_Biological_Process, GO_Molecular_Function, GO_Cellular_Component${NC}"
        echo -e "  ${DIM}See results/*/integrated_data/integrated_annotations_expression.tsv for format.${NC}"
        echo ""

        while true; do
            INTEGRATED=$(prompt "  Path to annotated expression table (TSV):")
            if [[ -f "$INTEGRATED" ]]; then
                n_genes=$(tail -n +2 "$INTEGRATED" | wc -l)
                # Check it has the minimum required columns
                header=$(head -1 "$INTEGRATED")
                has_seqid=false; has_tpm=false
                [[ "$header" == *"sequence_id"* ]] && has_seqid=true
                [[ "$header" == *"TPM"* ]] && has_tpm=true
                if $has_seqid && $has_tpm; then
                    success "Found ${n_genes} genes with required columns"
                    break
                else
                    error "File must contain 'sequence_id' and 'TPM' columns."
                    echo -e "  ${DIM}Found header: ${header:0:100}...${NC}"
                fi
            else
                error "File not found: ${INTEGRATED}"
            fi
        done

        # Category set
        echo ""
        echo -e "${BOLD}${CYAN}Category Set${NC}"
        cat_options=()
        for i in $(seq 1 15); do
            cat_options+=("${CATEGORY_DESCRIPTIONS[$i]}")
        done
        cat_idx=$(select_option "What type of organism are you studying?" "${cat_options[@]}")
        cat_num=$((cat_idx + 1))
        CATEGORY_SET="${CATEGORY_NAMES[$cat_num]}"
        success "Selected: ${CATEGORY_SET}"

        # Handle custom categories
        if [[ "$CATEGORY_SET" == "custom" ]]; then
            echo ""
            while true; do
                CUSTOM_FUNC=$(prompt "  Path to functional categories JSON:")
                if [[ -f "$CUSTOM_FUNC" ]]; then success "Found: ${CUSTOM_FUNC}"; break
                else error "File not found: ${CUSTOM_FUNC}"; fi
            done
        fi

        # Output settings
        echo ""
        echo -e "${BOLD}${CYAN}Output Settings${NC}\n"
        OUTDIR=$(prompt_default "  Output directory" "results")
        PREFIX=$(prompt_default "  Output file prefix" "sceptr")

        NF_CMD="nextflow run ${SCRIPT_DIR}/main.nf -entry ExPlotEntry"
        NF_CMD+=" --integrated_results \"${INTEGRATED}\""
        NF_CMD+=" --category_set ${CATEGORY_SET}"
        NF_CMD+=" --outdir \"${OUTDIR}\""
        NF_CMD+=" --output_prefix \"${PREFIX}\""
        NF_CMD+=" -profile ${PROFILE}"
        if [[ "$CATEGORY_SET" == "custom" && -n "${CUSTOM_FUNC:-}" ]]; then
            NF_CMD+=" --custom_functional_categories \"${CUSTOM_FUNC}\""
        fi

        echo ""
        echo -e "${CYAN}${BOLD}┌─────────────────────────────────────────────────────────────┐${NC}"
        echo -e "${CYAN}${BOLD}│  Method Only - Summary                                      │${NC}"
        echo -e "${CYAN}${BOLD}└─────────────────────────────────────────────────────────────┘${NC}"
        echo ""
        echo -e "  ${BOLD}Input:${NC}           ${INTEGRATED}"
        echo -e "  ${BOLD}Genes:${NC}           ${n_genes}"
        echo -e "  ${BOLD}Category set:${NC}    ${CATEGORY_SET}"
        echo -e "  ${BOLD}Output:${NC}          ${OUTDIR}/"
        echo ""
        echo -e "  ${DIM}Command: ${NF_CMD}${NC}"
        echo ""
        if confirm "  Launch enrichment profiling?"; then
            echo ""
            info "Launching SCEPTR enrichment profiling..."
            echo ""
            eval $NF_CMD
            exit $?
        else
            info "Cancelled."
            exit 0
        fi

    elif (( mode_idx == 3 )); then
        # ── Re-run enrichment profiling on existing SCEPTR results ──
        echo ""
        echo -e "${BOLD}${CYAN}Re-run Enrichment Profiling${NC}\n"
        echo -e "  ${DIM}Re-analyse existing SCEPTR results with different categories or settings.${NC}\n"
        while true; do
            INTEGRATED=$(prompt "  Path to integrated_annotations_expression.tsv:")
            if [[ -f "$INTEGRATED" ]]; then
                success "Found: ${INTEGRATED}"
                break
            else
                error "File not found: ${INTEGRATED}"
            fi
        done

        # Category set
        cat_options=()
        for i in $(seq 1 15); do
            cat_options+=("${CATEGORY_DESCRIPTIONS[$i]}")
        done
        cat_idx=$(select_option "What type of organism?" "${cat_options[@]}")
        cat_num=$((cat_idx + 1))
        CATEGORY_SET="${CATEGORY_NAMES[$cat_num]}"
        OUTDIR=$(prompt_default "  Output directory" "results")
        PREFIX=$(prompt_default "  Output file prefix" "sceptr")

        NF_CMD="nextflow run ${SCRIPT_DIR}/main.nf -entry ExPlotEntry"
        NF_CMD+=" --integrated_results \"${INTEGRATED}\""
        NF_CMD+=" --category_set ${CATEGORY_SET}"
        NF_CMD+=" --outdir \"${OUTDIR}\""
        NF_CMD+=" --output_prefix \"${PREFIX}\""
        NF_CMD+=" -profile ${PROFILE}"

        echo ""
        echo -e "  ${DIM}Command: ${NF_CMD}${NC}"
        echo ""
        if confirm "  Launch enrichment profiling?"; then
            eval $NF_CMD
            exit $?
        else
            info "Cancelled."
            exit 0
        fi
    fi

    # If we got here in compare mode, skip the full framework interactive flow
    if $COMPARE_MODE; then
        # Skip to validation/run (handled below)
        :
    else

    # Step 1: Reads
    echo -e "${BOLD}${CYAN}Step 1/4: RNA-seq Reads${NC}\n"

    # Ask about read layout
    layout_idx=$(select_option "  What type of reads do you have?" \
        "Paired-end (two files per sample, e.g. *_R1.fastq.gz + *_R2.fastq.gz)" \
        "Single-end (one file per sample)")

    if (( layout_idx == 1 )); then
        SINGLE_END=true
        echo ""
        echo "  Enter the path to your reads directory or a glob pattern."
        echo -e "  ${DIM}Examples: data/reads   or   data/*.fastq.gz${NC}\n"
    else
        echo ""
        echo "  Enter the path to your reads directory or a glob pattern."
        echo -e "  ${DIM}Examples: data/reads   or   data/*_{1,2}.fastq.gz${NC}\n"
    fi

    while true; do
        READS_INPUT=$(prompt "  Reads path:")

        # If it's a directory, try to auto-detect the pattern
        if [[ -d "$READS_INPUT" ]]; then
            if $SINGLE_END; then
                # Look for any fastq files
                n_files=$(ls "${READS_INPUT}"/*.f{ast,}q{,.gz} 2>/dev/null | wc -l)
                if (( n_files > 0 )); then
                    READS="${READS_INPUT}/*.fastq.gz"
                    # Try .fq.gz if no .fastq.gz
                    if (( $(ls "${READS_INPUT}"/*.fastq.gz 2>/dev/null | wc -l) == 0 )); then
                        READS="${READS_INPUT}/*.fq.gz"
                    fi
                    success "Found ${n_files} read files in ${READS_INPUT}"
                    break
                else
                    error "No FASTQ files found in ${READS_INPUT}"
                fi
            else
                detected=$(detect_reads "$READS_INPUT" || true)
                if [[ -n "$detected" ]]; then
                    n_files=$(count_read_files "$detected")
                    success "Auto-detected ${n_files} read files: ${detected}"
                    READS="$detected"
                    break
                else
                    error "No paired-end FASTQ files found in ${READS_INPUT}"
                    echo -e "  ${DIM}Expected patterns: *_1.fastq.gz/*_2.fastq.gz or *_R1.fastq.gz/*_R2.fastq.gz${NC}"
                fi
            fi
        elif [[ "$READS_INPUT" == *"*"* ]] || [[ "$READS_INPUT" == *"{"* ]]; then
            # It's a glob pattern
            n_files=$(count_read_files "$READS_INPUT")
            if (( n_files > 0 )); then
                success "Found ${n_files} read files matching pattern"
                READS="$READS_INPUT"
                break
            else
                error "No files matched pattern: ${READS_INPUT}"
            fi
        elif [[ -f "$READS_INPUT" ]]; then
            if $SINGLE_END; then
                READS="$READS_INPUT"
                success "Using single file: $(basename "$READS_INPUT")"
                break
            else
                # Single file - check if it's one of a pair
                dir=$(dirname "$READS_INPUT")
                detected=$(detect_reads "$dir" || true)
                if [[ -n "$detected" ]]; then
                    success "Auto-detected read pattern: ${detected}"
                    READS="$detected"
                    break
                else
                    READS="$READS_INPUT"
                    break
                fi
            fi
        else
            error "Path not found: ${READS_INPUT}"
        fi
    done

    # Step 2: Transcriptome
    echo ""
    echo -e "${BOLD}${CYAN}Step 2/4: Reference Transcriptome${NC}\n"
    echo "  Path to your transcriptome assembly (FASTA)."
    echo -e "  ${DIM}This is typically a de novo assembly or reference CDS file.${NC}\n"

    while true; do
        TRANSCRIPTS=$(prompt "  Transcriptome FASTA:")
        if [[ -f "$TRANSCRIPTS" ]]; then
            n_seqs=$(grep -c "^>" "$TRANSCRIPTS" 2>/dev/null || echo "?")
            success "Found transcriptome: ${n_seqs} sequences"
            break
        else
            error "File not found: ${TRANSCRIPTS}"
        fi
    done

    # Step 3: Organism type
    echo ""
    echo -e "${BOLD}${CYAN}Step 3/4: Organism Type${NC}"

    cat_options=()
    for i in $(seq 1 15); do
        cat_options+=("${CATEGORY_DESCRIPTIONS[$i]}")
    done

    cat_idx=$(select_option "What type of organism are you studying?" "${cat_options[@]}")
    cat_num=$((cat_idx + 1))
    CATEGORY_SET="${CATEGORY_NAMES[$cat_num]}"
    success "Selected: ${CATEGORY_SET}"

    # Handle custom categories
    if [[ "$CATEGORY_SET" == "custom" ]]; then
        echo ""
        while true; do
            CUSTOM_FUNC=$(prompt "  Path to functional categories JSON:")
            if [[ -f "$CUSTOM_FUNC" ]]; then
                success "Found: ${CUSTOM_FUNC}"
                break
            else
                error "File not found: ${CUSTOM_FUNC}"
            fi
        done
        while true; do
            CUSTOM_CELL=$(prompt "  Path to cellular categories JSON:")
            if [[ -f "$CUSTOM_CELL" ]]; then
                success "Found: ${CUSTOM_CELL}"
                break
            else
                error "File not found: ${CUSTOM_CELL}"
            fi
        done
    fi

    # Auto-enable direct translation for bacteria (CDS inputs)
    if [[ "$CATEGORY_SET" == "bacteria" || "$CATEGORY_SET" == "bacteria_gram_negative" || "$CATEGORY_SET" == "bacteria_gram_positive" ]]; then
        SKIP_TRANSDECODER=true
        echo ""
        echo -e "  ${YELLOW}ℹ  ${CATEGORY_SET} category detected - will use direct CDS translation instead of TransDecoder.${NC}"
        echo -e "  ${DIM}   (Override with --skip_transdecoder false if using a de novo transcriptome assembly)${NC}"
    fi

    # Step 4: Host filtering
    echo ""
    echo -e "${BOLD}${CYAN}Step 4/4: Host Contamination Filtering${NC}\n"

    # Suggest host filtering for relevant category sets
    suggest_host=false
    case "$CATEGORY_SET" in
        parasite_protozoan|helminth_nematode|helminth_platyhelminth)
            suggest_host=true
            echo -e "  ${YELLOW}You selected ${CATEGORY_SET} - host contamination filtering is recommended.${NC}"
            ;;
        *)
            echo "  Optional: remove host sequences from parasite/pathogen data."
            ;;
    esac

    if confirm "  Do you have a host reference to filter against?" "${suggest_host:+y}"; then
        echo ""
        host_type_idx=$(select_option "  What type of host reference do you have?" \
            "Transcriptome (nucleotide FASTA - will be translated automatically)" \
            "Proteome (protein FASTA - faster, skips translation)")

        if (( host_type_idx == 0 )); then
            HOST_TYPE="transcriptome"
        else
            HOST_TYPE="proteome"
        fi

        while true; do
            HOST_REF=$(prompt "  Path to host reference:")
            if [[ -f "$HOST_REF" ]]; then
                n_seqs=$(grep -c "^>" "$HOST_REF" 2>/dev/null || echo "?")
                success "Found host reference: ${n_seqs} sequences"
                break
            else
                error "File not found: ${HOST_REF}"
            fi
        done
    fi

    # Output settings
    echo ""
    echo -e "${BOLD}${CYAN}Output Settings${NC}\n"
    OUTDIR=$(prompt_default "  Output directory" "results")
    PREFIX=$(prompt_default "  Output file prefix" "sceptr")

    fi  # end of non-compare interactive flow
fi

# ── Handle comparison mode ───────────────────────────────────────────────

if $COMPARE_MODE; then
    OUTDIR="${OUTDIR:-results_comparison}"
    PREFIX="${PREFIX:-sceptr_comparison}"
    CATEGORY_SET="${CATEGORY_SET:-general}"
    N_PERMUTATIONS="${N_PERMUTATIONS:-10000}"

    # Validate comparison inputs
    echo ""
    errors=0

    if [[ -z "$CONDITION_A" ]]; then
        error "No condition A specified (use --condition-a)"
        errors=$((errors + 1))
    elif [[ ! -f "$CONDITION_A" ]]; then
        error "Condition A file not found: ${CONDITION_A}"
        errors=$((errors + 1))
    fi

    if [[ -z "$CONDITION_B" ]]; then
        error "No condition B specified (use --condition-b)"
        errors=$((errors + 1))
    elif [[ ! -f "$CONDITION_B" ]]; then
        error "Condition B file not found: ${CONDITION_B}"
        errors=$((errors + 1))
    fi

    if (( errors > 0 )); then
        error "Cannot proceed - please fix the errors above."
        exit 1
    fi

    # Build comparison command
    NF_CMD="nextflow run ${SCRIPT_DIR}/main.nf -entry compare"
    NF_CMD+=" --condition_a \"${CONDITION_A}\""
    NF_CMD+=" --condition_b \"${CONDITION_B}\""
    NF_CMD+=" --category_set ${CATEGORY_SET}"
    NF_CMD+=" --outdir \"${OUTDIR}\""
    NF_CMD+=" --output_prefix \"${PREFIX}\""
    NF_CMD+=" --n_permutations ${N_PERMUTATIONS}"
    NF_CMD+=" -profile ${PROFILE}"

    if [[ -n "$LABEL_A" ]]; then
        NF_CMD+=" --label_a \"${LABEL_A}\""
    fi
    if [[ -n "$LABEL_B" ]]; then
        NF_CMD+=" --label_b \"${LABEL_B}\""
    fi
    if [[ -n "$TIERS" ]]; then
        NF_CMD+=" --expression_tiers \"${TIERS}\""
    fi
    if [[ "$CATEGORY_SET" == "custom" && -n "${CUSTOM_FUNC:-}" ]]; then
        NF_CMD+=" --custom_functional_categories \"${CUSTOM_FUNC}\""
    fi
    if $RESUME; then
        NF_CMD+=" -resume"
    fi
    if [[ -n "$EXTRA_ARGS" ]]; then
        NF_CMD+=" ${EXTRA_ARGS}"
    fi

    # Summary
    echo ""
    echo -e "${CYAN}${BOLD}┌─────────────────────────────────────────────────────────────┐${NC}"
    echo -e "${CYAN}${BOLD}│  Comparison Summary                                         │${NC}"
    echo -e "${CYAN}${BOLD}└─────────────────────────────────────────────────────────────┘${NC}"
    echo ""
    echo -e "  ${BOLD}Condition A:${NC}     ${CONDITION_A}"
    echo -e "  ${BOLD}Label A:${NC}         ${LABEL_A:-Condition_A}"
    echo -e "  ${BOLD}Condition B:${NC}     ${CONDITION_B}"
    echo -e "  ${BOLD}Label B:${NC}         ${LABEL_B:-Condition_B}"
    echo -e "  ${BOLD}Category set:${NC}    ${CATEGORY_SET}"
    echo -e "  ${BOLD}Permutations:${NC}    ${N_PERMUTATIONS}"
    echo -e "  ${BOLD}Output:${NC}          ${OUTDIR}/"
    echo -e "  ${BOLD}Profile:${NC}         ${PROFILE}"
    echo ""
    echo -e "  ${DIM}Command:${NC}"
    echo -e "  ${DIM}${NF_CMD}${NC}"
    echo ""

    if $INTERACTIVE; then
        if ! confirm "  Launch comparison?"; then
            info "Cancelled. You can re-run with the command above."
            exit 0
        fi
    fi

    echo ""
    info "Launching SCEPTR cross-sample comparison..."
    echo ""

    eval $NF_CMD
    exit_code=$?

    echo ""
    if (( exit_code == 0 )); then
        echo -e "${GREEN}${BOLD}┌─────────────────────────────────────────────────────────────┐${NC}"
        echo -e "${GREEN}${BOLD}│  Comparison completed successfully!                         │${NC}"
        echo -e "${GREEN}${BOLD}└─────────────────────────────────────────────────────────────┘${NC}"
        echo ""
        echo -e "  ${BOLD}Results:${NC} ${OUTDIR}/"
        echo ""
        for report in "${OUTDIR}"/*comparison*report*.html "${OUTDIR}"/comparison/*report*.html; do
            if [[ -f "$report" ]]; then
                echo -e "  ${GREEN}→${NC} ${report}"
            fi
        done
        echo ""
        echo -e "  ${DIM}Sláinte a chara!${NC}"
        echo ""
    else
        echo -e "${RED}${BOLD}SCEPTR failed (exit code: ${exit_code})${NC}"
        echo "  Check .nextflow.log for details."
        exit $exit_code
    fi
    exit 0
fi

# ── Set defaults for anything not provided ───────────────────────────────────

OUTDIR="${OUTDIR:-results}"
PREFIX="${PREFIX:-sceptr}"

# ── Validate inputs ──────────────────────────────────────────────────────────

echo ""
errors=0

if [[ -z "$READS" ]]; then
    error "No reads specified (use -r or --reads)"
    errors=$((errors + 1))
fi

if [[ -z "$TRANSCRIPTS" ]]; then
    error "No transcriptome specified (use -t or --transcripts)"
    errors=$((errors + 1))
elif [[ ! -f "$TRANSCRIPTS" ]]; then
    error "Transcriptome file not found: ${TRANSCRIPTS}"
    errors=$((errors + 1))
fi

if [[ -z "$CATEGORY_SET" ]]; then
    CATEGORY_SET="general"
    warn "No category set specified, using 'general'"
fi

if [[ -n "$HOST_REF" && ! -f "$HOST_REF" ]]; then
    error "Host reference not found: ${HOST_REF}"
    errors=$((errors + 1))
fi

# Check Docker is available
if [[ "$PROFILE" == "docker" ]]; then
    if ! command -v docker &> /dev/null; then
        error "Docker not found. Install Docker or use --profile singularity"
        errors=$((errors + 1))
    elif ! docker image inspect sceptr:${VERSION} &> /dev/null 2>&1; then
        warn "Docker image sceptr:${VERSION} not found. Building..."
        if [[ -f "${SCRIPT_DIR}/Dockerfile" ]]; then
            docker build -t sceptr:${VERSION} "${SCRIPT_DIR}" || {
                error "Docker build failed"
                errors=$((errors + 1))
            }
        else
            error "Dockerfile not found. Build the image first: docker build -t sceptr:${VERSION} ."
            errors=$((errors + 1))
        fi
    fi
fi

# Check Nextflow is available
if ! command -v nextflow &> /dev/null; then
    error "Nextflow not found. Install: curl -s https://get.nextflow.io | bash"
    errors=$((errors + 1))
fi

if (( errors > 0 )); then
    echo ""
    error "Cannot proceed - please fix the errors above."
    exit 1
fi

# Auto-enable direct translation for bacteria in CLI mode
if [[ "$CATEGORY_SET" == "bacteria" || "$CATEGORY_SET" == "bacteria_gram_negative" || "$CATEGORY_SET" == "bacteria_gram_positive" ]] && ! $SKIP_TRANSDECODER; then
    SKIP_TRANSDECODER=true
fi

# ── Build Nextflow command ───────────────────────────────────────────────────

NF_CMD="nextflow run ${SCRIPT_DIR}/main.nf"
NF_CMD+=" --reads \"${READS}\""
NF_CMD+=" --transcripts \"${TRANSCRIPTS}\""
NF_CMD+=" --category_set ${CATEGORY_SET}"
NF_CMD+=" --outdir \"${OUTDIR}\""
NF_CMD+=" --output_prefix \"${PREFIX}\""
NF_CMD+=" -profile ${PROFILE}"

if $SINGLE_END; then
    NF_CMD+=" --single_end true"
fi

if $SKIP_TRANSDECODER; then
    NF_CMD+=" --skip_transdecoder true"
fi

if [[ -n "$HOST_REF" ]]; then
    if [[ "$HOST_TYPE" == "proteome" ]]; then
        NF_CMD+=" --host_proteome \"${HOST_REF}\""
    else
        NF_CMD+=" --host_transcriptome \"${HOST_REF}\""
    fi
fi

if [[ -n "$TIERS" ]]; then
    NF_CMD+=" --expression_tiers \"${TIERS}\""
fi

if [[ "$CATEGORY_SET" == "custom" ]]; then
    NF_CMD+=" --custom_functional_categories \"${CUSTOM_FUNC}\""
    NF_CMD+=" --custom_cellular_categories \"${CUSTOM_CELL}\""
fi

if $RESUME; then
    NF_CMD+=" -resume"
fi

if [[ -n "$EXTRA_ARGS" ]]; then
    NF_CMD+=" ${EXTRA_ARGS}"
fi

# ── Display summary and confirm ──────────────────────────────────────────────

echo ""
echo -e "${CYAN}${BOLD}┌─────────────────────────────────────────────────────────────┐${NC}"
echo -e "${CYAN}${BOLD}│  Run Summary                                                │${NC}"
echo -e "${CYAN}${BOLD}└─────────────────────────────────────────────────────────────┘${NC}"
echo ""
echo -e "  ${BOLD}Reads:${NC}           ${READS}"
echo -e "  ${BOLD}Read layout:${NC}     $( $SINGLE_END && echo 'Single-end' || echo 'Paired-end' )"
echo -e "  ${BOLD}Translation:${NC}     $( $SKIP_TRANSDECODER && echo 'Direct CDS translation (TransDecoder skipped)' || echo 'TransDecoder ORF prediction' )"
echo -e "  ${BOLD}Transcriptome:${NC}   ${TRANSCRIPTS}"
echo -e "  ${BOLD}Category set:${NC}    ${CATEGORY_SET}"
if [[ -n "$HOST_REF" ]]; then
echo -e "  ${BOLD}Host reference:${NC}  ${HOST_REF} (${HOST_TYPE})"
fi
echo -e "  ${BOLD}Output:${NC}          ${OUTDIR}/"
echo -e "  ${BOLD}Prefix:${NC}          ${PREFIX}"
echo -e "  ${BOLD}Profile:${NC}         ${PROFILE}"
echo ""
echo -e "  ${DIM}Command:${NC}"
echo -e "  ${DIM}${NF_CMD}${NC}"
echo ""

if $INTERACTIVE; then
    if ! confirm "  Launch SCEPTR?"; then
        info "Cancelled. You can re-run with the command above."
        exit 0
    fi
fi

# ── Run ──────────────────────────────────────────────────────────────────────

echo ""
info "Launching SCEPTR..."
echo ""

eval $NF_CMD
exit_code=$?

echo ""
if (( exit_code == 0 )); then

    # Organise outputs: promote combined report, tidy intermediates
    if [[ -d "${OUTDIR}" ]]; then
        # Move combined report to top level of results
        if [[ -f "${OUTDIR}/enrichment_profiles/sceptr_report.html" ]]; then
            cp "${OUTDIR}/enrichment_profiles/sceptr_report.html" "${OUTDIR}/sceptr_report.html" 2>/dev/null
        fi
        for prefix_file in "${OUTDIR}/enrichment_profiles/"*"_report.html"; do
            [[ -f "$prefix_file" ]] || continue
            bn=$(basename "$prefix_file")
            [[ "$bn" == *"_BP_MF_report.html" || "$bn" == *"_CC_report.html" ]] || continue
            # These are already in the combined report
        done

        # Move redundant standalone reports to intermediate/ folder
        mkdir -p "${OUTDIR}/intermediate" 2>/dev/null
        for redundant in \
            "${OUTDIR}/enrichment_profiles/functional/"*"_BP_MF_report.html" \
            "${OUTDIR}/enrichment_profiles/cellular/"*"_CC_report.html" \
            "${OUTDIR}/functional/"*"_BP_MF_report.html" \
            "${OUTDIR}/landscape/"*"_landscape_report.html"; do
            if [[ -f "$redundant" ]]; then
                mv "$redundant" "${OUTDIR}/intermediate/" 2>/dev/null
            fi
        done

        # Also move the duplicate combined report from enrichment_profiles/
        # (the top-level copy is now the canonical one)
        if [[ -f "${OUTDIR}/sceptr_report.html" && -f "${OUTDIR}/enrichment_profiles/sceptr_report.html" ]]; then
            mv "${OUTDIR}/enrichment_profiles/sceptr_report.html" "${OUTDIR}/intermediate/" 2>/dev/null
        fi
    fi

    echo -e "${GREEN}${BOLD}┌─────────────────────────────────────────────────────────────┐${NC}"
    echo -e "${GREEN}${BOLD}│  SCEPTR completed successfully!                    │${NC}"
    echo -e "${GREEN}${BOLD}└─────────────────────────────────────────────────────────────┘${NC}"
    echo ""
    echo -e "  ${BOLD}Results:${NC} ${OUTDIR}/"
    echo ""
    echo -e "  ${BOLD}Report:${NC}"
    if [[ -f "${OUTDIR}/sceptr_report.html" ]]; then
        echo -e "    ${GREEN}→${NC} ${OUTDIR}/sceptr_report.html"
    fi
    if [[ -f "${OUTDIR}/go_enrichment/reports/"*".html" ]]; then
        for r in "${OUTDIR}"/go_enrichment/reports/*.html; do
            [[ -f "$r" ]] && echo -e "    ${GREEN}→${NC} ${r}"
        done
    fi
    echo ""
    echo -e "  ${DIM}Sláinte a chara!${NC}"
    echo ""
else
    echo -e "${RED}${BOLD}┌─────────────────────────────────────────────────────────────┐${NC}"
    echo -e "${RED}${BOLD}│  SCEPTR failed (exit code: ${exit_code})                          │${NC}"
    echo -e "${RED}${BOLD}└─────────────────────────────────────────────────────────────┘${NC}"
    echo ""
    echo "  Check .nextflow.log for details, or re-run with:"
    echo "  ${NF_CMD} -resume"
    echo ""
    exit $exit_code
fi
