"""
Script 03: Codon Usage and RSCU Analysis
==========================================
Computes Relative Synonymous Codon Usage (RSCU), Effective Number of
Codons (ENC), and Codon Adaptation Index (CAI) for each cancer gene.

Author: Sudiksha
Project: GC Content Variability in Cancer Genes
Date: April 2025 - March 2026
"""

import os
import csv
import math
import numpy as np
from Bio import SeqIO
from collections import Counter

# ============================================================
# PATHS
# ============================================================
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
FASTA_FILE = os.path.join(DATA_DIR, "all_cancer_genes_cds.fasta")

# ============================================================
# GENETIC CODE TABLE
# Standard genetic code (NCBI translation table 1)
# ============================================================
CODON_TABLE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Amino acid to synonymous codons mapping
AA_TO_CODONS = {}
for codon, aa in CODON_TABLE.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(codon)

# Degeneracy classes
DEGENERACY = {aa: len(codons) for aa, codons in AA_TO_CODONS.items()}

# Human highly expressed gene reference codon usage (per 1000 codons)
# Source: Codon Usage Database (Kazusa) - Homo sapiens
HUMAN_REFERENCE_CODON_USAGE = {
    "TTT": 17.6, "TTC": 20.3, "TTA": 7.7, "TTG": 12.9,
    "CTT": 13.2, "CTC": 19.6, "CTA": 7.2, "CTG": 39.6,
    "ATT": 16.0, "ATC": 20.8, "ATA": 7.5, "ATG": 22.0,
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1,
    "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4,
    "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "TAT": 12.2, "TAC": 15.3, "TAA": 1.0, "TAG": 0.8,
    "CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2,
    "AAT": 17.0, "AAC": 19.1, "AAA": 24.4, "AAG": 31.9,
    "GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
    "TGT": 10.6, "TGC": 12.6, "TGA": 1.6, "TGG": 13.2,
    "CGT": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4,
    "AGT": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
}


def count_codons(seq):
    """Count all codons in a coding sequence."""
    seq = str(seq).upper()
    codons = [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    return Counter(codons)


def calculate_rscu(codon_counts):
    """
    Calculate Relative Synonymous Codon Usage (RSCU).
    RSCU(ij) = X(ij) / [(1/n_i) × Σ X(ij)]
    """
    rscu = {}

    for aa, codons in AA_TO_CODONS.items():
        # Skip non-degenerate amino acids (Met, Trp)
        if len(codons) == 1:
            rscu[codons[0]] = 1.0
            continue

        # Total count for this amino acid
        total = sum(codon_counts.get(c, 0) for c in codons)
        n_synonymous = len(codons)

        for codon in codons:
            observed = codon_counts.get(codon, 0)
            expected = total / n_synonymous if n_synonymous > 0 else 0
            rscu[codon] = observed / expected if expected > 0 else 0

    return rscu


def calculate_enc(codon_counts):
    """
    Calculate Effective Number of Codons (ENC).
    ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6

    where Fk = average homozygosity for k-fold degenerate amino acids.
    """
    # Group amino acids by degeneracy
    deg_groups = {}
    for aa, codons in AA_TO_CODONS.items():
        deg = len(codons)
        if deg > 1:  # Skip non-degenerate
            deg_groups.setdefault(deg, []).append(aa)

    def homozygosity(aa):
        """Calculate homozygosity (F) for an amino acid."""
        codons = AA_TO_CODONS[aa]
        counts = [codon_counts.get(c, 0) for c in codons]
        n = sum(counts)
        if n <= 1:
            return None  # Not enough data
        f = sum((c / n) ** 2 for c in counts)
        # Corrected F: F_hat = (n * F - 1) / (n - 1)
        f_corrected = (n * f - 1) / (n - 1)
        return max(f_corrected, 1.0 / len(codons))  # minimum is uniform usage

    # Calculate average F for each degeneracy class
    f_values = {}
    for deg, aas in deg_groups.items():
        f_list = [homozygosity(aa) for aa in aas if homozygosity(aa) is not None]
        if f_list:
            f_values[deg] = np.mean(f_list)
        else:
            f_values[deg] = 1.0 / deg  # Default to uniform

    # ENC formula
    enc = 2  # Met + Trp (non-degenerate)
    enc += 9 / f_values.get(2, 0.5)       # 9 two-fold degenerate AAs
    enc += 1 / f_values.get(3, 0.333)     # 1 three-fold (Ile)
    enc += 5 / f_values.get(4, 0.25)      # 5 four-fold degenerate AAs
    enc += 3 / f_values.get(6, 0.167)     # 3 six-fold degenerate AAs

    # Clamp to valid range [20, 61]
    return max(20, min(61, enc))


def calculate_enc_expected(gc3s):
    """
    Calculate expected ENC under pure mutational pressure.
    ENC_expected = 2 + s + 29 / (s^2 + (1-s)^2)
    where s = GC3s (as fraction, not percentage)
    """
    s = gc3s / 100.0
    denom = s ** 2 + (1 - s) ** 2
    if denom == 0:
        return 61
    return 2 + s + 29 / denom


def calculate_cai(codon_counts, reference=HUMAN_REFERENCE_CODON_USAGE):
    """
    Calculate Codon Adaptation Index (CAI).
    CAI = geometric mean of w(i) for all codons in the gene.
    w(i) = RSCU(i) / RSCU(max) for each amino acid.
    """
    # Calculate reference RSCU
    ref_rscu = {}
    for aa, codons in AA_TO_CODONS.items():
        if len(codons) == 1:
            continue
        total = sum(reference.get(c, 0) for c in codons)
        for codon in codons:
            ref_rscu[codon] = reference.get(codon, 0) / (total / len(codons)) if total > 0 else 0

    # Calculate w values (relative adaptiveness)
    w_values = {}
    for aa, codons in AA_TO_CODONS.items():
        if len(codons) == 1:
            continue
        max_rscu = max(ref_rscu.get(c, 0) for c in codons)
        for codon in codons:
            w_values[codon] = ref_rscu.get(codon, 0) / max_rscu if max_rscu > 0 else 0

    # Calculate CAI
    log_sum = 0
    count = 0
    for codon, cnt in codon_counts.items():
        if codon in w_values and w_values[codon] > 0 and CODON_TABLE.get(codon) not in ("*", "M", "W"):
            log_sum += cnt * math.log(w_values[codon])
            count += cnt

    if count == 0:
        return 0
    return math.exp(log_sum / count)


def analyze_gene_codons(record):
    """Perform complete codon usage analysis for a single gene."""
    seq = str(record.seq).upper()
    gene_name = record.id

    # Count codons (exclude stop codon at end)
    cds = seq[:-3] if seq[-3:] in ("TAA", "TAG", "TGA") else seq
    codon_counts = count_codons(cds)

    # Calculate metrics
    rscu = calculate_rscu(codon_counts)
    enc = calculate_enc(codon_counts)
    cai = calculate_cai(codon_counts)

    # Calculate GC3s for ENC expected
    codons_list = [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]
    exclude = {"ATG", "TGG"}
    syn_codons = [c for c in codons_list if c not in exclude and len(c) == 3]
    gc3s = sum(1 for c in syn_codons if c[2] in "GC") / len(syn_codons) * 100 if syn_codons else 50

    enc_expected = calculate_enc_expected(gc3s)

    return {
        "gene_name": gene_name,
        "codon_counts": codon_counts,
        "rscu": rscu,
        "enc": enc,
        "enc_expected": enc_expected,
        "cai": cai,
        "gc3s": gc3s,
        "total_codons": sum(codon_counts.values()),
    }


def save_rscu_table(all_results, output_file):
    """Save RSCU values for all genes to a TSV file."""
    # Get all codons (sorted)
    sense_codons = sorted([c for c, aa in CODON_TABLE.items() if aa != "*"])

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Header
        header = ["Codon", "Amino_Acid"] + [r["gene_name"] for r in all_results] + ["Mean", "SD"]
        writer.writerow(header)

        for codon in sense_codons:
            aa = CODON_TABLE[codon]
            rscu_values = [r["rscu"].get(codon, 0) for r in all_results]
            row = [codon, aa] + [f"{v:.3f}" for v in rscu_values]
            row += [f"{np.mean(rscu_values):.3f}", f"{np.std(rscu_values):.3f}"]
            writer.writerow(row)


def save_enc_data(all_results, output_file):
    """Save ENC analysis data."""
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Gene", "ENC", "ENC_Expected", "GC3s", "CAI", "Total_Codons",
                         "ENC_Deviation", "Below_Curve"])

        for r in all_results:
            deviation = r["enc_expected"] - r["enc"]
            below = "Yes" if r["enc"] < r["enc_expected"] * 0.90 else "No"
            writer.writerow([
                r["gene_name"],
                f"{r['enc']:.2f}",
                f"{r['enc_expected']:.2f}",
                f"{r['gc3s']:.2f}",
                f"{r['cai']:.4f}",
                r["total_codons"],
                f"{deviation:.2f}",
                below,
            ])


def identify_preferred_avoided_codons(all_results):
    """Identify overall preferred and avoided codons across all genes."""
    print("\n" + "=" * 70)
    print("PREFERRED AND AVOIDED CODONS (across all 50 cancer genes)")
    print("=" * 70)

    # Calculate mean RSCU for each codon
    sense_codons = [c for c, aa in CODON_TABLE.items() if aa not in ("*", "M", "W")]
    mean_rscu = {}
    for codon in sense_codons:
        values = [r["rscu"].get(codon, 0) for r in all_results]
        mean_rscu[codon] = np.mean(values)

    # Preferred codons (RSCU > 1.6)
    preferred = [(c, mean_rscu[c], CODON_TABLE[c]) for c in sense_codons if mean_rscu[c] > 1.6]
    preferred.sort(key=lambda x: x[1], reverse=True)

    print("\nOverrepresented codons (RSCU > 1.6):")
    print(f"  {'Codon':<8} {'AA':<6} {'RSCU':<8} {'Type'}")
    print(f"  {'-'*35}")
    for codon, rscu, aa in preferred:
        ending = "GC-ending" if codon[2] in "GC" else "AT-ending"
        print(f"  {codon:<8} {aa:<6} {rscu:<8.3f} {ending}")

    # Avoided codons (RSCU < 0.6)
    avoided = [(c, mean_rscu[c], CODON_TABLE[c]) for c in sense_codons if mean_rscu[c] < 0.6]
    avoided.sort(key=lambda x: x[1])

    print(f"\nUnderrepresented codons (RSCU < 0.6):")
    print(f"  {'Codon':<8} {'AA':<6} {'RSCU':<8} {'Type'}")
    print(f"  {'-'*35}")
    for codon, rscu, aa in avoided:
        ending = "GC-ending" if codon[2] in "GC" else "AT-ending"
        print(f"  {codon:<8} {aa:<6} {rscu:<8.3f} {ending}")


def print_enc_summary(all_results):
    """Print ENC summary by gene category."""
    print("\n" + "=" * 70)
    print("ENC ANALYSIS SUMMARY")
    print("=" * 70)

    for r in sorted(all_results, key=lambda x: x["enc"]):
        status = "STRONG BIAS" if r["enc"] < 35 else "MODERATE" if r["enc"] < 45 else "LOW BIAS"
        print(f"  {r['gene_name']:<10} ENC={r['enc']:5.1f}  GC3s={r['gc3s']:5.1f}%  "
              f"CAI={r['cai']:.3f}  [{status}]")


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("CODON USAGE ANALYSIS (RSCU, ENC, CAI)")
    print("=" * 70)

    if not os.path.exists(FASTA_FILE):
        print(f"ERROR: Input file not found: {FASTA_FILE}")
        print("Please run 01_data_acquisition.py first.")
        return

    # Analyze each gene
    all_results = []
    for record in SeqIO.parse(FASTA_FILE, "fasta"):
        print(f"  Analyzing {record.id}...")
        result = analyze_gene_codons(record)
        all_results.append(result)

    # Save results
    rscu_file = os.path.join(DATA_DIR, "rscu_all_genes.tsv")
    save_rscu_table(all_results, rscu_file)
    print(f"\nRSCU table saved: {rscu_file}")

    enc_file = os.path.join(DATA_DIR, "enc_analysis.tsv")
    save_enc_data(all_results, enc_file)
    print(f"ENC data saved: {enc_file}")

    # Print summaries
    identify_preferred_avoided_codons(all_results)
    print_enc_summary(all_results)

    print("\nCodon usage analysis complete!")


if __name__ == "__main__":
    main()
