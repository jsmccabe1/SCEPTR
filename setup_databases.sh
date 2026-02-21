#!/usr/bin/env bash
#
# SCEPTR Database Setup Script
# Downloads and prepares required databases for the SCEPTR pipeline.
#
# Databases:
#   1. UniProt Swiss-Prot (reviewed proteins) → DIAMOND database
#   2. Contaminant database (common contaminants) → DIAMOND database
#   3. Gene Ontology (GO) hierarchy file
#
# Total download: ~500 MB compressed → ~3.5 GB on disk after building
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
        info "diamond not found locally — will use Docker (buchfink/diamond)"
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
    info "UniProt DIAMOND database already exists — skipping (delete to rebuild)"
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
if [[ -f "${DATA_DIR}/contaminants/contaminants_uniprot.dmnd" ]]; then
    info "Contaminant database already exists — skipping (delete to rebuild)"
else
    info "Building contaminant database from UniProt taxonomic groups..."
    info "Downloading bacterial proteins (Swiss-Prot, reviewed)..."
    wget -q --show-progress -O "${DATA_DIR}/contaminants/bacteria.fasta.gz" \
        "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28taxonomy_id%3A2%29+AND+%28reviewed%3Atrue%29"

    info "Downloading fungal proteins (Swiss-Prot, reviewed)..."
    wget -q --show-progress -O "${DATA_DIR}/contaminants/fungi.fasta.gz" \
        "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28taxonomy_id%3A4751%29+AND+%28reviewed%3Atrue%29"

    info "Downloading viral proteins (Swiss-Prot, reviewed)..."
    wget -q --show-progress -O "${DATA_DIR}/contaminants/viruses.fasta.gz" \
        "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28taxonomy_id%3A10239%29+AND+%28reviewed%3Atrue%29"

    info "Combining and building contaminant DIAMOND database..."
    zcat "${DATA_DIR}/contaminants/bacteria.fasta.gz" \
         "${DATA_DIR}/contaminants/fungi.fasta.gz" \
         "${DATA_DIR}/contaminants/viruses.fasta.gz" > "${DATA_DIR}/contaminants/contaminants_combined.fasta"

    if $USE_DOCKER_DIAMOND; then
        docker run --rm -v "${DATA_DIR}:/data" \
            buchfink/diamond:latest \
            diamond makedb --in /data/contaminants/contaminants_combined.fasta \
            --db /data/contaminants/contaminants_uniprot -p 4
    else
        diamond makedb --in "${DATA_DIR}/contaminants/contaminants_combined.fasta" \
            --db "${DATA_DIR}/contaminants/contaminants_uniprot" -p 4
    fi

    # Record stats
    bact_count=$(zgrep -c "^>" "${DATA_DIR}/contaminants/bacteria.fasta.gz" || true)
    fungi_count=$(zgrep -c "^>" "${DATA_DIR}/contaminants/fungi.fasta.gz" || true)
    virus_count=$(zgrep -c "^>" "${DATA_DIR}/contaminants/viruses.fasta.gz" || true)
    total_count=$(grep -c "^>" "${DATA_DIR}/contaminants/contaminants_combined.fasta" || true)

    cat > "${DATA_DIR}/contaminants/contaminant_db_stats.txt" << EOF
SCEPTR Contaminant Database Statistics
======================================
Date built: $(date -u +"%Y-%m-%d")
Source: UniProt Swiss-Prot (reviewed)

Bacteria (taxonomy:2):     ${bact_count} proteins
Fungi (taxonomy:4751):     ${fungi_count} proteins
Viruses (taxonomy:10239):  ${virus_count} proteins
Total:                     ${total_count} proteins
EOF

    # Clean up intermediate files
    rm -f "${DATA_DIR}/contaminants/bacteria.fasta.gz" \
          "${DATA_DIR}/contaminants/fungi.fasta.gz" \
          "${DATA_DIR}/contaminants/viruses.fasta.gz" \
          "${DATA_DIR}/contaminants/contaminants_combined.fasta"

    success "Contaminant database built ($(du -sh "${DATA_DIR}/contaminants/contaminants_uniprot.dmnd" | cut -f1))"
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
