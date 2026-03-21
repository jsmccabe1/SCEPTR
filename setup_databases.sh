#!/usr/bin/env bash
#
# SCEPTR Database Setup Script
# Downloads and prepares required databases for the SCEPTR pipeline.
#
# Databases:
#   1. UniProt Swiss-Prot (reviewed proteins) → DIAMOND database
#   2. Contaminant database (57 curated species across 5 categories) → DIAMOND database
#   3. Gene Ontology (GO) hierarchy file
#
# Total download: ~2 GB compressed → ~4 GB on disk after building
#
# Usage:
#   bash setup_databases.sh              # Download all databases
#   bash setup_databases.sh --check      # Check if databases are present
#   bash setup_databases.sh --go-only    # Only update GO hierarchy
#

set -euo pipefail

# ── Colours ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

info()    { echo -e "${BLUE}ℹ${NC}  $1"; }
success() { echo -e "${GREEN}✓${NC}  $1"; }
warn()    { echo -e "${YELLOW}⚠${NC}  $1"; }
error()   { echo -e "${RED}✗${NC}  $1"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"

# ── Check mode ───────────────────────────────────────────────────────────────
if [[ "${1:-}" == "--check" ]]; then
    echo -e "\n${CYAN}${BOLD}SCEPTR Database Status${NC}\n"
    all_ok=true

    if [[ -f "${DATA_DIR}/uniprot/uniprot.dmnd" ]]; then
        size=$(du -sh "${DATA_DIR}/uniprot/uniprot.dmnd" | cut -f1)
        success "UniProt DIAMOND database: ${size}"
    else
        warn "UniProt DIAMOND database: MISSING"
        all_ok=false
    fi

    if [[ -f "${DATA_DIR}/contaminants/contaminants_uniprot.dmnd" ]]; then
        size=$(du -sh "${DATA_DIR}/contaminants/contaminants_uniprot.dmnd" | cut -f1)
        success "Contaminant DIAMOND database: ${size}"
    else
        warn "Contaminant DIAMOND database: MISSING"
        all_ok=false
    fi

    if [[ -f "${DATA_DIR}/go/go-basic.obo" ]]; then
        size=$(du -sh "${DATA_DIR}/go/go-basic.obo" | cut -f1)
        success "GO hierarchy (go-basic.obo): ${size}"
    else
        warn "GO hierarchy: MISSING"
        all_ok=false
    fi

    echo ""
    if $all_ok; then
        success "All databases present. SCEPTR is ready to run."
    else
        warn "Some databases are missing. Run: bash setup_databases.sh"
    fi
    exit 0
fi

# ── Banner ───────────────────────────────────────────────────────────────────
echo ""
echo -e "${CYAN}${BOLD}┌─────────────────────────────────────────────────────────────┐${NC}"
echo -e "${CYAN}${BOLD}│   SCEPTR Database Setup                                     │${NC}"
echo -e "${CYAN}${BOLD}│   Downloading and building required databases                │${NC}"
echo -e "${CYAN}${BOLD}└─────────────────────────────────────────────────────────────┘${NC}"
echo ""

# ── Check dependencies ───────────────────────────────────────────────────────
USE_DOCKER_DIAMOND=false
if ! command -v wget &> /dev/null; then
    error "wget is required but not installed. Please install it first."
fi

if ! command -v diamond &> /dev/null; then
    if command -v docker &> /dev/null; then
        info "diamond not found locally - will use Docker (buchfink/diamond)"
        USE_DOCKER_DIAMOND=true
    else
        error "diamond is required but not installed. Install DIAMOND or Docker."
    fi
fi

# ── Create directory structure ───────────────────────────────────────────────
mkdir -p "${DATA_DIR}/uniprot"
mkdir -p "${DATA_DIR}/contaminants"
mkdir -p "${DATA_DIR}/go"

# ── 1. Gene Ontology ────────────────────────────────────────────────────────
if [[ "${1:-}" == "--go-only" ]] || [[ ! -f "${DATA_DIR}/go/go-basic.obo" ]]; then
    info "Downloading Gene Ontology hierarchy..."
    wget -q --show-progress -O "${DATA_DIR}/go/go-basic.obo" \
        "https://purl.obolibrary.org/obo/go/go-basic.obo"
    success "GO hierarchy downloaded ($(du -sh "${DATA_DIR}/go/go-basic.obo" | cut -f1))"
fi

if [[ "${1:-}" == "--go-only" ]]; then
    success "GO update complete."
    exit 0
fi

# ── 2. UniProt Swiss-Prot ───────────────────────────────────────────────────
if [[ -f "${DATA_DIR}/uniprot/uniprot.dmnd" ]]; then
    info "UniProt DIAMOND database already exists - skipping (delete to rebuild)"
else
    info "Downloading UniProt Swiss-Prot (reviewed proteins)..."
    wget -q --show-progress -O "${DATA_DIR}/uniprot/uniprot_sprot.fasta.gz" \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

    info "Decompressing..."
    gunzip -f "${DATA_DIR}/uniprot/uniprot_sprot.fasta.gz"

    info "Building DIAMOND database (this may take a few minutes)..."
    if $USE_DOCKER_DIAMOND; then
        docker run --rm -v "${DATA_DIR}:/data" \
            buchfink/diamond:latest \
            diamond makedb --in /data/uniprot/uniprot_sprot.fasta \
            --db /data/uniprot/uniprot -p 4
    else
        diamond makedb --in "${DATA_DIR}/uniprot/uniprot_sprot.fasta" \
            --db "${DATA_DIR}/uniprot/uniprot" -p 4
    fi

    success "UniProt DIAMOND database built ($(du -sh "${DATA_DIR}/uniprot/uniprot.dmnd" | cut -f1))"
fi

# ── 3. Contaminant Database ─────────────────────────────────────────────────
# Downloads complete proteomes (reviewed + unreviewed) from UniProt for
# 57 curated species across 5 categories: bacteria, viruses, mammals, fungi,
# plants. These species represent the most common sources of contamination
# in transcriptome sequencing experiments.
if [[ -f "${DATA_DIR}/contaminants/contaminants_uniprot.dmnd" ]]; then
    info "Contaminant database already exists - skipping (delete to rebuild)"
else
    info "Building contaminant database from curated species proteomes..."
    info "This downloads ~2 GB and may take 15-30 minutes depending on connection speed."
    echo ""

    CONTAM_DIR="${DATA_DIR}/contaminants"

    # Helper: build UniProt REST query from taxonomy IDs and download
    download_proteomes() {
        local output="$1"
        shift
        local ids=("$@")

        local query=""
        for id in "${ids[@]}"; do
            [[ -n "$query" ]] && query="${query}+OR+"
            query="${query}taxonomy_id%3A${id}"
        done

        wget -q --show-progress -O "$output" \
            "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28${query}%29"

        # Verify download is not empty
        if [[ ! -s "$output" ]]; then
            warn "Download appears empty: $output"
            return 1
        fi
    }

    # ── Bacteria (24 species) ──
    # E. coli, S. aureus, M. tuberculosis, P. aeruginosa, K. pneumoniae,
    # B. subtilis, H. pylori, C. jejuni, E. faecalis, E. faecium,
    # B. cereus, B. thuringiensis, L. pneumophila, S. marcescens,
    # M. pneumoniae, S. epidermidis, B. pseudomallei, B. cepacia,
    # S. maltophilia, B. diazoefficiens, B. japonicum, C. acnes,
    # B. fragilis, E. cloacae
    info "[1/5] Downloading bacterial proteomes (24 species)..."
    download_proteomes "${CONTAM_DIR}/bacteria.fasta.gz" \
        562 1280 1773 287 573 1423 210 197 1351 1352 \
        1396 1428 446 615 2104 1282 28450 292 \
        40324 114615 375 1496 817 550

    # ── Viruses (10 species) ──
    # HHV-6A, HHV-6B, HHV-7, MCMV, phiX174, AcMNPV, EBV, HIV-1, HIV-2, MLV
    info "[2/5] Downloading viral proteomes (10 species)..."
    download_proteomes "${CONTAM_DIR}/viruses.fasta.gz" \
        32603 32604 10372 10366 10847 46015 10376 11676 11709 11786

    # ── Mammals (9 species) ──
    # Human, mouse, rat, macaque, chimpanzee, cow, pig, dog, rabbit
    info "[3/5] Downloading mammalian proteomes (9 species)..."
    download_proteomes "${CONTAM_DIR}/mammals.fasta.gz" \
        9606 10090 10116 9544 9598 9913 9823 9615 9986

    # ── Fungi (8 species) ──
    # A. fumigatus, A. flavus, A. niger, P. chrysogenum, P. rubens,
    # C. albicans, S. cerevisiae, N. crassa
    info "[4/5] Downloading fungal proteomes (8 species)..."
    download_proteomes "${CONTAM_DIR}/fungi.fasta.gz" \
        5085 5059 5061 5076 500485 5476 4932 5141

    # ── Plants (6 species) ──
    # Arabidopsis, maize, rice, tomato, tobacco, wheat
    info "[5/5] Downloading plant proteomes (6 species)..."
    download_proteomes "${CONTAM_DIR}/plants.fasta.gz" \
        3702 4577 4530 4081 4097 4565

    info "Combining proteomes and building DIAMOND database..."
    zcat "${CONTAM_DIR}/bacteria.fasta.gz" \
         "${CONTAM_DIR}/viruses.fasta.gz" \
         "${CONTAM_DIR}/mammals.fasta.gz" \
         "${CONTAM_DIR}/fungi.fasta.gz" \
         "${CONTAM_DIR}/plants.fasta.gz" > "${CONTAM_DIR}/contaminants_combined.fasta"

    if $USE_DOCKER_DIAMOND; then
        docker run --rm -v "${DATA_DIR}:/data" \
            buchfink/diamond:latest \
            diamond makedb --in /data/contaminants/contaminants_combined.fasta \
            --db /data/contaminants/contaminants_uniprot -p 4
    else
        diamond makedb --in "${CONTAM_DIR}/contaminants_combined.fasta" \
            --db "${CONTAM_DIR}/contaminants_uniprot" -p 4
    fi

    # Record stats
    bact_count=$(zgrep -c "^>" "${CONTAM_DIR}/bacteria.fasta.gz" || true)
    virus_count=$(zgrep -c "^>" "${CONTAM_DIR}/viruses.fasta.gz" || true)
    mammal_count=$(zgrep -c "^>" "${CONTAM_DIR}/mammals.fasta.gz" || true)
    fungi_count=$(zgrep -c "^>" "${CONTAM_DIR}/fungi.fasta.gz" || true)
    plant_count=$(zgrep -c "^>" "${CONTAM_DIR}/plants.fasta.gz" || true)
    total_count=$(grep -c "^>" "${CONTAM_DIR}/contaminants_combined.fasta" || true)

    cat > "${CONTAM_DIR}/contaminant_db_stats.txt" << STATS
Total sequences: ${total_count}
Bacteria: ['Mycobacterium tuberculosis', 'Stenotrophomonas maltophilia', 'Bradyrhizobium diazoefficiens', 'Burkholderia pseudomallei', 'Klebsiella pneumoniae', 'Legionella pneumophila', 'Pseudomonas aeruginosa', 'Staphylococcus epidermidis', 'Staphylococcus aureus', 'Bacillus thuringiensis', 'Bradyrhizobium japonicum', 'Burkholderia cepacia', 'Helicobacter pylori', 'Serratia marcescens', 'Mycoplasma pneumoniae', 'Campylobacter jejuni', 'Enterococcus faecalis', 'Enterococcus faecium', 'Bacillus cereus', 'Bacillus subtilis', 'Escherichia coli', 'Cutibacterium acnes', 'Bacteroides fragilis', 'Enterobacter cloacae'] sequences
Viruses: ['Human betaherpesvirus 6A', 'Human betaherpesvirus 6B', 'Human betaherpesvirus 7', 'Murid betaherpesvirus 1', 'Escherichia virus phiX174', 'Autographa californica multiple nucleopolyhedrovirus', 'Human gammaherpesvirus 4', 'Human immunodeficiency virus 1', 'Human immunodeficiency virus 2', 'Murine leukemia virus'] sequences
Mammals: ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Macaca mulatta', 'Pan troglodytes', 'Bos taurus', 'Sus scrofa', 'Canis lupus familiaris', 'Oryctolagus cuniculus'] sequences
Fungi: ['Aspergillus fumigatus', 'Aspergillus flavus', 'Aspergillus niger', 'Penicillium chrysogenum', 'Penicillium rubens', 'Candida albicans', 'Saccharomyces cerevisiae', 'Neurospora crassa'] sequences
Plants: ['Arabidopsis thaliana', 'Zea mays', 'Oryza sativa', 'Solanum lycopersicum', 'Nicotiana tabacum', 'Triticum aestivum'] sequences
STATS

    echo ""
    info "Contaminant database composition:"
    echo "    Bacteria:  ${bact_count} proteins (24 species)"
    echo "    Viruses:   ${virus_count} proteins (10 species)"
    echo "    Mammals:   ${mammal_count} proteins (9 species)"
    echo "    Fungi:     ${fungi_count} proteins (8 species)"
    echo "    Plants:    ${plant_count} proteins (6 species)"
    echo "    Total:     ${total_count} proteins"

    # Clean up intermediate files
    rm -f "${CONTAM_DIR}/bacteria.fasta.gz" \
          "${CONTAM_DIR}/viruses.fasta.gz" \
          "${CONTAM_DIR}/mammals.fasta.gz" \
          "${CONTAM_DIR}/fungi.fasta.gz" \
          "${CONTAM_DIR}/plants.fasta.gz" \
          "${CONTAM_DIR}/contaminants_combined.fasta"

    success "Contaminant database built ($(du -sh "${CONTAM_DIR}/contaminants_uniprot.dmnd" | cut -f1))"
fi

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo -e "${GREEN}${BOLD}Database setup complete!${NC}"
echo ""
echo "  UniProt:       $(du -sh "${DATA_DIR}/uniprot/uniprot.dmnd" | cut -f1)"
echo "  Contaminants:  $(du -sh "${DATA_DIR}/contaminants/contaminants_uniprot.dmnd" | cut -f1)"
echo "  GO hierarchy:  $(du -sh "${DATA_DIR}/go/go-basic.obo" | cut -f1)"
echo "  Total:         $(du -sh "${DATA_DIR}" | cut -f1)"
echo ""
info "You can now run SCEPTR:"
echo -e "  ${CYAN}./run_sceptr.sh${NC}"
echo ""
