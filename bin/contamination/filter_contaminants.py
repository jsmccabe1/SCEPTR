#!/usr/bin/env python3
"""
SCEPTR DeCon: Taxonomy-aware contamination filtering for predicted proteomes.

Filters contaminant sequences from DIAMOND BLAST results using taxonomy-aware
identity thresholds. Different organism groups require different thresholds
because sequence conservation varies by evolutionary distance.

Rationale for taxonomy-specific thresholds:
  - Bacterial/viral sequences at moderate identity (50-70%) against a eukaryotic
    query are very likely contaminants, since prokaryote-eukaryote divergence is
    ancient. Lower thresholds catch genuine contamination without risking
    false positives from conserved domains.
  - Fungal sequences require intermediate thresholds (70%) as fungi share more
    recent ancestry with other eukaryotes.
  - Mammalian/plant sequences at moderate identity are likely just conserved
    orthologs, not contaminants. Only near-identical hits (>=90%) indicate
    genuine cross-contamination from sample handling.

Coverage filtering ensures short domain-level alignments (e.g. a conserved
ATP-binding domain) don't trigger whole-sequence removal.
"""

import sys
import re
import csv
import json
import argparse
from Bio import SeqIO
from collections import defaultdict


# =============================================================================
# Taxonomy classification using word-boundary matching
# =============================================================================

# Each entry: (compiled_regex, is_case_sensitive)
# Word-boundary matching prevents substring false positives
# e.g. "coli" won't match "Brucella melitensis" and "rice" won't match "price"
TAXONOMY_PATTERNS = {
    "Viruses": re.compile(
        r'\b(?:virus|viral|viridae|phage|bacteriophage|retrovir|coronavir|herpesvir|'
        r'adenovir|parvovir|papillomavir|polyomavir|norovir|influenza|rhinovir|'
        r'enterovir|rotavir|flavivir|orbivir|bunyavir|arenavirus|filovir|'
        r'togavir|rhabdovir|paramyxovir|orthomyxovir|picornavir|calicivir|'
        r'astrovir|reovir|birnavir|hepadnavir|circo[v]ir|gemini[v]ir|'
        r'phiX174|lambda\s+phage)\b',
        re.IGNORECASE
    ),
    "Bacteria": re.compile(
        r'\b(?:Escherichia|Bacillus|Pseudomonas|Salmonella|Mycoplasma|'
        r'Staphylococcus|Streptococcus|Clostridium|Lactobacillus|'
        r'Mycobacterium|Listeria|Neisseria|Corynebacterium|Enterobacter|'
        r'Acinetobacter|Vibrio|Shigella|Chlamydia|Borrelia|Rickettsia|'
        r'Proteus|Yersinia|Klebsiella|Legionella|Helicobacter|Serratia|'
        r'Campylobacter|Enterococcus|Cutibacterium|Bacteroides|Thermus|'
        r'Deinococcus|Stenotrophomonas|Bradyrhizobium|Burkholderia|'
        r'Haemophilus|Bordetella|Treponema|Leptospira|Brucella|'
        r'Francisella|Bartonella|Coxiella|Agrobacterium|Rhizobium|'
        r'Caulobacter|Synechocystis|Synechococcus|Anabaena|Nostoc|'
        r'bacteria[l]?|actinobacteri|firmicutes|proteobacteri|'
        r'cyanobacteri|bacteroidetes|spirochaet)\b',
        re.IGNORECASE
    ),
    "Fungi": re.compile(
        r'\b(?:Saccharomyces|Aspergillus|Candida|Cryptococcus|Penicillium|'
        r'Neurospora|Fusarium|Alternaria|Cladosporium|Trichophyton|'
        r'Mucor|Rhizopus|Trichoderma|Malassezia|Rhodotorula|Geotrichum|'
        r'Schizosaccharomyces|Yarrowia|Pichia|Kluyveromyces|Botrytis|'
        r'Magnaporthe|Ustilago|Puccinia|Coprinopsis|Agaricus|'
        r'fung[iu]s|fung[ia]l|yeast|ascomycet|basidiomycet|zygomycet)\b',
        re.IGNORECASE
    ),
    "Mammals": re.compile(
        r'\b(?:Homo\s+sapiens|Mus\s+musculus|Rattus\s+norvegicus|'
        r'Macaca\s+mulatta|Pan\s+troglodytes|Bos\s+taurus|Sus\s+scrofa|'
        r'Canis\s+lupus|Oryctolagus\s+cuniculus|Ovis\s+aries|'
        r'Capra\s+hircus|Equus\s+caballus|Gorilla\s+gorilla|'
        r'Pongo\s+abelii|Felis\s+catus|Cavia\s+porcellus|'
        r'Cricetulus\s+griseus|Mesocricetus\s+auratus|'
        r'mammali[a]|primat[e]s|rodenti[a]|carnivora|'
        r'artiodactyl|perissodactyl|chiroptera|lagomorpha)\b',
        re.IGNORECASE
    ),
    "Plants": re.compile(
        r'\b(?:Arabidopsis|Nicotiana|Oryza\s+sativa|Zea\s+mays|'
        r'Solanum|Triticum|Hordeum|Brassica|Glycine\s+max|'
        r'Gossypium|Citrus|Vitis|Helianthus|Cucumis|Medicago|'
        r'Populus|Pinus|Eucalyptus|Sorghum|Avena|Secale|'
        r'Pisum|Phaseolus|Linum|Cannabis|Capsicum|Malus|Prunus|'
        r'Fragaria|Spinacia|Beta\s+vulgaris|Daucus|Lactuca|'
        r'viridiplantae|streptophyt[a]|embryophyt[a]|'
        r'angiosperm|gymnosperm|monocot|dicot|magnoliophyt)\b',
        re.IGNORECASE
    ),
}


def classify_taxonomy(taxonomy_str):
    """Classify a taxonomy/organism string into a contamination category.
    
    Uses word-boundary regex matching to avoid substring false positives.
    Returns the first matching category, checked in order of specificity:
    Viruses > Bacteria > Fungi > Mammals > Plants > Other.
    """
    if not taxonomy_str:
        return "Other"
    
    for category, pattern in TAXONOMY_PATTERNS.items():
        if pattern.search(taxonomy_str):
            return category
    
    return "Other"


# =============================================================================
# Threshold logic
# =============================================================================

# Default thresholds - overridden by command-line arguments
DEFAULT_THRESHOLDS = {
    "Bacteria":  70.0,
    "Viruses":   50.0,
    "Fungi":     70.0,
    "Mammals":   90.0,
    "Plants":    90.0,
    "Other":     50.0,
}

def get_threshold(category, default_identity):
    """Get identity threshold for a given taxonomy category."""
    return DEFAULT_THRESHOLDS.get(category, default_identity)


# =============================================================================
# Conserved protein family detection
# =============================================================================
# Ultra-conserved protein families show high cross-kingdom identity due to
# ancient functional constraints, not contamination. Ubiquitin is ~96% identical
# from yeast to human; histones H3/H4 are >90% across all eukaryotes; tubulins
# and actins are similarly conserved (Theodosius Dobzhansky, 1973; Ciccarelli
# et al., 2006, Science 311:1283-1287).
#
# For these families, only near-identical hits (>=99%) are flagged, since
# genuine cross-contamination would produce exact or near-exact matches,
# whereas orthologs plateau at their family-specific conservation ceiling.

CONSERVED_PROTEIN_PATTERN = re.compile(
    r'\b(?:'
    # Ubiquitin monomer/fusion only (NOT ubiquitin-conjugating, ubiquitin-ligase, etc.)
    r'[Pp]olyubiquitin|UBB\b|UBC\b|UBA52\b|RPS27A\b|'
    r'[Uu]biquitin-ribosomal|[Uu]biquitin-40S|[Uu]biquitin-60S|'
    # Histones (core and linker)
    r'[Hh]istone\s*H[1234]|H2A[._]|H2B[._]|H3[._]|H4[._]|'
    r'HIST[0-9]|CENPA|'
    # Tubulins
    r'[Tt]ubulin\s*(alpha|beta|gamma)|TUBA[0-9]|TUBB[0-9]|TUBG[0-9]|'
    # Actins (but not actinobacter, actinoporin, etc.)
    r'[Aa]ctin\s+(alpha|beta|gamma)|ACTB\b|ACTG\b|ACT[0-9]\b|'
    # Heat shock / chaperones
    r'[Hh]eat\s*shock\s*protein|HSP[0-9]|[Hh]sp[679]0|chaperonin|'
    r'GroEL|GroES|CCT[0-9]|TCP1|'
    # Ribosomal proteins (multiple naming conventions)
    r'[0-9]+[Ss]\s*ribosomal\s*protein|ribosomal\s*protein\s*[eSuL]+[0-9]|'
    r'ribosomal\s*subunit\s*protein|'
    # Elongation / translation factors
    r'[Ee]longation\s*factor\s*[12]|EF-?1\b|EF-?2\b|EEF[12]\b|ETF1\b|eIF-?[0-9]|'
    r'[Tt]ranslation\s*initiation\s*factor|'
    # Calmodulin / EF-hand
    r'[Cc]almodulin|EF-hand|'
    # Cyclophilins / peptidyl-prolyl isomerases
    r'[Cc]yclophilin|peptidyl-prolyl|PPIA\b|PPIB\b|'
    # ADP/ATP carrier
    r'ADP.ATP\s*carrier|adenine\s*nucleotide\s*transloca[a-z]*'
    r')\b',
    re.IGNORECASE
)

CONSERVED_THRESHOLD = 99.0  # Only flag if essentially identical


def is_conserved_protein(description):
    """Check if a BLAST hit description matches a known ultra-conserved protein family."""
    if not description:
        return False
    return bool(CONSERVED_PROTEIN_PATTERN.search(description))


# =============================================================================
# Main filtering logic
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="SCEPTR DeCon: Taxonomy-aware contamination filtering"
    )
    parser.add_argument("--identity", type=float, default=50.0,
                        help="Default minimum percent identity (default: 50.0)")
    parser.add_argument("--coverage", type=float, default=30.0,
                        help="Minimum query coverage percent (default: 30.0)")
    parser.add_argument("--evalue", type=float, default=1e-3,
                        help="Maximum e-value (default: 1e-3)")
    parser.add_argument("--proteome", required=True,
                        help="Input proteome FASTA file")
    parser.add_argument("--blast", default="contaminants_hits.tsv",
                        help="DIAMOND BLAST results file (outfmt 6)")
    parser.add_argument("--bacterial-threshold", type=float, default=70.0,
                        help="Identity threshold for bacterial hits (default: 70.0)")
    parser.add_argument("--viral-threshold", type=float, default=50.0,
                        help="Identity threshold for viral hits (default: 50.0)")
    parser.add_argument("--fungal-threshold", type=float, default=70.0,
                        help="Identity threshold for fungal hits (default: 70.0)")
    parser.add_argument("--eukaryotic-threshold", type=float, default=90.0,
                        help="Identity threshold for mammalian/plant hits (default: 90.0)")
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Update thresholds from command-line arguments
    DEFAULT_THRESHOLDS["Bacteria"] = args.bacterial_threshold
    DEFAULT_THRESHOLDS["Viruses"] = args.viral_threshold
    DEFAULT_THRESHOLDS["Fungi"] = args.fungal_threshold
    DEFAULT_THRESHOLDS["Mammals"] = args.eukaryotic_threshold
    DEFAULT_THRESHOLDS["Plants"] = args.eukaryotic_threshold
    DEFAULT_THRESHOLDS["Other"] = args.identity
    
    # -------------------------------------------------------------------------
    # Parse BLAST results
    # -------------------------------------------------------------------------
    # Expected DIAMOND outfmt:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen stitle
    
    all_hits = defaultdict(list)
    total_blast_hits = 0
    
    print("SCEPTR DeCon Contamination Filtering")
    print("=" * 40)
    print(f"\nParsing BLAST results from: {args.blast}")
    
    try:
        with open(args.blast) as f:
            for line in f:
                if not line.strip():
                    continue
                total_blast_hits += 1
                parts = line.strip().split("\t")
                
                if len(parts) < 14:
                    continue
                
                query_id = parts[0].split()[0]
                subject_id = parts[1]
                percent_identity = float(parts[2])
                alignment_length = int(parts[3])
                query_start = int(parts[6])
                query_end = int(parts[7])
                e_value = float(parts[10])
                bit_score = float(parts[11])
                query_length = int(parts[12])
                taxonomy = parts[13] if len(parts) > 13 else "Unknown"
                
                # Calculate query coverage
                alignment_span = abs(query_end - query_start) + 1
                query_coverage = (alignment_span / query_length * 100) if query_length > 0 else 0
                
                all_hits[query_id].append({
                    'subject_id': subject_id,
                    'percent_identity': percent_identity,
                    'alignment_length': alignment_length,
                    'query_coverage': query_coverage,
                    'query_length': query_length,
                    'e_value': e_value,
                    'bit_score': bit_score,
                    'taxonomy': taxonomy,
                })
    except FileNotFoundError:
        print(f"Warning: BLAST results file not found: {args.blast}")
        print("Assuming no contaminants detected.")
    
    print(f"Total BLAST hits parsed: {total_blast_hits}")
    print(f"Unique query sequences with hits: {len(all_hits)}")
    
    # -------------------------------------------------------------------------
    # Apply taxonomy-aware filtering
    # -------------------------------------------------------------------------
    # Check ALL hits per query (not just best hit) against their respective
    # category thresholds. A sequence is flagged as a contaminant if ANY hit
    # passes the threshold for its taxonomy category.
    
    final_contaminants = {}
    conserved_skipped = []  # Track conserved proteins that were spared
    
    for query_id, hits in all_hits.items():
        best_contaminant_hit = None
        best_contaminant_score = -1
        
        for hit in hits:
            # Classify this hit's taxonomy
            category = classify_taxonomy(hit['taxonomy'])
            threshold = get_threshold(category, args.identity)
            
            # Check if hit matches a conserved protein family
            conserved = is_conserved_protein(hit['taxonomy'])
            if conserved:
                threshold = CONSERVED_THRESHOLD
            
            # Apply all three filters: identity, coverage, e-value
            if (hit['percent_identity'] >= threshold and
                hit['query_coverage'] >= args.coverage and
                hit['e_value'] <= args.evalue):
                
                # Keep the strongest contaminant hit for reporting
                if hit['bit_score'] > best_contaminant_score:
                    best_contaminant_hit = dict(hit)
                    best_contaminant_hit['category'] = category
                    best_contaminant_hit['threshold_applied'] = threshold
                    best_contaminant_hit['conserved_protein'] = conserved
                    best_contaminant_score = hit['bit_score']
            
            elif conserved and hit['percent_identity'] >= get_threshold(category, args.identity):
                # Would have been flagged without conserved protein protection
                conserved_skipped.append({
                    'query_id': query_id,
                    'percent_identity': hit['percent_identity'],
                    'taxonomy': hit['taxonomy'][:80],
                    'category': category,
                })
        
        if best_contaminant_hit is not None:
            final_contaminants[query_id] = best_contaminant_hit
    
    print(f"\nContaminant sequences identified: {len(final_contaminants)}")
    if conserved_skipped:
        print(f"Conserved protein hits excluded: {len(conserved_skipped)} "
              f"(would have been flagged without conserved-family protection)")
    
    # -------------------------------------------------------------------------
    # Group by category for reporting
    # -------------------------------------------------------------------------
    category_groups = defaultdict(list)
    for seq_id, hit_info in final_contaminants.items():
        category_groups[hit_info['category']].append({
            'seq_id': seq_id,
            'percent_identity': hit_info['percent_identity'],
            'query_coverage': hit_info['query_coverage'],
            'e_value': hit_info['e_value'],
            'bit_score': hit_info['bit_score'],
            'taxonomy': hit_info['taxonomy'],
            'threshold_applied': hit_info['threshold_applied'],
        })
    
    # -------------------------------------------------------------------------
    # Write text report
    # -------------------------------------------------------------------------
    with open("contaminant_report.txt", "w") as report:
        report.write("SCEPTR DeCon Contamination Report\n")
        report.write("=" * 40 + "\n\n")
        
        report.write("Filtering Parameters:\n")
        report.write(f"  Default identity threshold: {args.identity}%\n")
        report.write(f"  Bacterial identity threshold: {args.bacterial_threshold}%\n")
        report.write(f"  Viral identity threshold: {args.viral_threshold}%\n")
        report.write(f"  Fungal identity threshold: {args.fungal_threshold}%\n")
        report.write(f"  Eukaryotic identity threshold: {args.eukaryotic_threshold}%\n")
        report.write(f"  Minimum query coverage: {args.coverage}%\n")
        report.write(f"  Maximum e-value: {args.evalue}\n\n")
        
        report.write(f"BLAST results: {total_blast_hits} hits for {len(all_hits)} sequences\n")
        report.write(f"Contaminant sequences identified: {len(final_contaminants)}\n\n")
        
        if category_groups:
            report.write("Category Summary:\n")
            for cat in sorted(category_groups.keys(), key=lambda c: len(category_groups[c]), reverse=True):
                seqs = category_groups[cat]
                report.write(f"  {cat}: {len(seqs)} sequences\n")
            
            report.write("\nDetailed Contaminant List:\n")
            report.write("-" * 40 + "\n")
            for cat in sorted(category_groups.keys()):
                report.write(f"\n{cat} (threshold: {get_threshold(cat, args.identity)}%):\n")
                for seq in sorted(category_groups[cat], key=lambda x: x['percent_identity'], reverse=True):
                    report.write(f"  {seq['seq_id']}: {seq['percent_identity']:.1f}% identity, "
                                f"{seq['query_coverage']:.1f}% coverage, "
                                f"E={seq['e_value']:.2e}\n")
                    report.write(f"    Hit: {seq['taxonomy'][:80]}\n")
        else:
            report.write("No contaminants identified.\n")
        
        # Report conserved proteins that were protected
        if conserved_skipped:
            report.write(f"\nConserved Protein Exclusions ({len(conserved_skipped)} hits):\n")
            report.write("-" * 40 + "\n")
            report.write("The following hits exceeded category thresholds but were excluded\n")
            report.write("because they match ultra-conserved protein families (threshold\n")
            report.write(f"elevated to {CONSERVED_THRESHOLD}% for these families):\n\n")
            for entry in sorted(conserved_skipped, key=lambda x: x['percent_identity'], reverse=True):
                report.write(f"  {entry['query_id']}: {entry['percent_identity']:.1f}% identity "
                            f"({entry['category']})\n")
                report.write(f"    {entry['taxonomy']}\n")
    
    # -------------------------------------------------------------------------
    # Write CSV for visualisation
    # -------------------------------------------------------------------------
    with open("contaminant_details.csv", "w", newline='') as csvfile:
        fieldnames = ['sequence_id', 'category', 'taxonomy', 'percent_identity',
                      'query_coverage', 'e_value', 'bit_score', 'threshold_applied',
                      'conserved_protein']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for seq_id, hit in final_contaminants.items():
            writer.writerow({
                'sequence_id': seq_id,
                'category': hit['category'],
                'taxonomy': hit['taxonomy'],
                'percent_identity': hit['percent_identity'],
                'query_coverage': hit['query_coverage'],
                'e_value': hit['e_value'],
                'bit_score': hit['bit_score'],
                'threshold_applied': hit['threshold_applied'],
                'conserved_protein': hit.get('conserved_protein', False),
            })
    
    # Write JSON with parameters for the visualisation script to read
    params_for_report = {
        'default_identity': args.identity,
        'bacterial_threshold': args.bacterial_threshold,
        'viral_threshold': args.viral_threshold,
        'fungal_threshold': args.fungal_threshold,
        'eukaryotic_threshold': args.eukaryotic_threshold,
        'conserved_threshold': CONSERVED_THRESHOLD,
        'coverage_threshold': args.coverage,
        'evalue_threshold': args.evalue,
        'total_blast_hits': total_blast_hits,
        'total_queries_with_hits': len(all_hits),
        'conserved_proteins_excluded': len(conserved_skipped),
    }
    with open("filtering_params.json", "w") as f:
        json.dump(params_for_report, f, indent=2)
    
    # -------------------------------------------------------------------------
    # Filter FASTA
    # -------------------------------------------------------------------------
    print("\nFiltering FASTA file...")
    kept_count = 0
    filtered_count = 0
    
    category_handles = {}

    try:
        with open(args.proteome) as infile, \
             open("filtered_proteome.fasta", "w") as filtered_out, \
             open("contaminant_sequences.fasta", "w") as contaminant_out:

            for record in SeqIO.parse(infile, "fasta"):
                seq_id = record.id.split()[0]

                if seq_id in final_contaminants:
                    category = final_contaminants[seq_id]['category']
                    SeqIO.write(record, contaminant_out, "fasta")

                    # Write to per-category file
                    cat_filename = f"{category.lower()}_contaminants.fasta"
                    if category not in category_handles:
                        category_handles[category] = open(cat_filename, "w")
                    SeqIO.write(record, category_handles[category], "fasta")

                    filtered_count += 1
                else:
                    SeqIO.write(record, filtered_out, "fasta")
                    kept_count += 1
    finally:
        for fh in category_handles.values():
            fh.close()
    
    # Append final counts to report
    with open("contaminant_report.txt", "a") as report:
        report.write(f"\nFinal Results:\n")
        report.write(f"  Sequences kept: {kept_count}\n")
        report.write(f"  Sequences filtered: {filtered_count}\n")
        report.write(f"  Removal rate: {filtered_count/(kept_count+filtered_count)*100:.2f}%\n")
    
    print(f"\nResults:")
    print(f"  Kept: {kept_count}")
    print(f"  Filtered: {filtered_count}")
    print(f"  Output: filtered_proteome.fasta")


if __name__ == "__main__":
    main()
