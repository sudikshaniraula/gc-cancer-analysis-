"""
Script 02: GC Content and Nucleotide Composition Analysis
==========================================================
Computes overall GC content, positional GC content (GC1, GC2, GC3),
individual nucleotide frequencies, and CpG O/E ratios for each gene.

Author: Sudiksha
Project: GC Content Variability in Cancer Genes
Date: April 2025 - March 2026
"""

import os
import csv
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from collections import Counter

# ============================================================
# PATHS
# ============================================================
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "data")
FASTA_FILE = os.path.join(DATA_DIR, "all_cancer_genes_cds.fasta")


def calculate_nucleotide_frequencies(seq):
    """Calculate individual nucleotide frequencies."""
    seq = str(seq).upper()
    length = len(seq)
    counts = Counter(seq)
    return {
        "A": counts.get("A", 0),
        "T": counts.get("T", 0),
        "G": counts.get("G", 0),
        "C": counts.get("C", 0),
        "A%": counts.get("A", 0) / length * 100,
        "T%": counts.get("T", 0) / length * 100,
        "G%": counts.get("G", 0) / length * 100,
        "C%": counts.get("C", 0) / length * 100,
        "length": length,
    }


def calculate_gc_content(seq):
    """Calculate overall GC content percentage."""
    seq = str(seq).upper()
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100


def calculate_positional_gc(seq):
    """
    Calculate GC content at each codon position (GC1, GC2, GC3).
    Assumes the sequence starts at the first codon position.
    """
    seq = str(seq).upper()
    codons = [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    n_codons = len(codons)

    if n_codons == 0:
        return {"GC1": 0, "GC2": 0, "GC3": 0}

    gc1 = sum(1 for c in codons if c[0] in "GC") / n_codons * 100
    gc2 = sum(1 for c in codons if c[1] in "GC") / n_codons * 100
    gc3 = sum(1 for c in codons if c[2] in "GC") / n_codons * 100

    return {"GC1": gc1, "GC2": gc2, "GC3": gc3, "GC12": (gc1 + gc2) / 2}


def calculate_gc3s(seq):
    """
    Calculate GC3s (GC at third synonymous position).
    Excludes Met (ATG) and Trp (TGG) codons which have no synonymous alternatives.
    Also excludes stop codons.
    """
    seq = str(seq).upper()
    codons = [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]

    # Non-degenerate codons to exclude
    exclude = {"ATG", "TGG", "TAA", "TAG", "TGA"}

    synonymous_codons = [c for c in codons if c not in exclude and len(c) == 3]
    if not synonymous_codons:
        return 0

    gc3s = sum(1 for c in synonymous_codons if c[2] in "GC") / len(synonymous_codons) * 100
    return gc3s


def calculate_nucleotide_skew(seq):
    """Calculate GC-skew and AT-skew."""
    seq = str(seq).upper()
    counts = Counter(seq)
    g, c = counts.get("G", 0), counts.get("C", 0)
    a, t = counts.get("A", 0), counts.get("T", 0)

    gc_skew = (g - c) / (g + c) if (g + c) > 0 else 0
    at_skew = (a - t) / (a + t) if (a + t) > 0 else 0

    return {"GC_skew": gc_skew, "AT_skew": at_skew}


def calculate_cpg_oe(seq):
    """
    Calculate CpG observed/expected ratio.
    CpG O/E = (f(CpG) × N) / (f(C) × f(G))
    """
    seq = str(seq).upper()
    length = len(seq)
    counts = Counter(seq)

    # Count CpG dinucleotides
    cpg_count = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "CG")

    c_count = counts.get("C", 0)
    g_count = counts.get("G", 0)

    if c_count == 0 or g_count == 0:
        return 0

    cpg_observed = cpg_count / (length - 1)
    cpg_expected = (c_count / length) * (g_count / length)

    return cpg_observed / cpg_expected if cpg_expected > 0 else 0


def calculate_dinucleotide_frequencies(seq):
    """Calculate all 16 dinucleotide frequencies."""
    seq = str(seq).upper()
    dinucs = [seq[i:i+2] for i in range(len(seq) - 1)]
    total = len(dinucs)
    counts = Counter(dinucs)

    bases = "ATGC"
    freqs = {}
    for b1 in bases:
        for b2 in bases:
            dinuc = b1 + b2
            freqs[dinuc] = counts.get(dinuc, 0) / total * 100 if total > 0 else 0

    return freqs


def sliding_window_gc(seq, window_size=100, step_size=25):
    """
    Compute GC content across the sequence using a sliding window.
    Returns list of (position, gc_content) tuples.
    """
    seq = str(seq).upper()
    results = []

    for i in range(0, len(seq) - window_size + 1, step_size):
        window = seq[i:i + window_size]
        gc = (window.count("G") + window.count("C")) / window_size * 100
        position = i + window_size // 2  # center of window
        results.append((position, gc))

    return results


def calculate_pr2_bias(seq):
    """
    Calculate PR2 (Parity Rule 2) bias values at third codon position.
    Returns A3/(A3+T3) and G3/(G3+C3).
    """
    seq = str(seq).upper()
    codons = [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]

    # Count bases at third position
    third_bases = [c[2] for c in codons if len(c) == 3]
    counts = Counter(third_bases)

    a3, t3 = counts.get("A", 0), counts.get("T", 0)
    g3, c3 = counts.get("G", 0), counts.get("C", 0)

    at_bias = a3 / (a3 + t3) if (a3 + t3) > 0 else 0.5
    gc_bias = g3 / (g3 + c3) if (g3 + c3) > 0 else 0.5

    return {"A3/(A3+T3)": at_bias, "G3/(G3+C3)": gc_bias}


def analyze_all_genes(fasta_file):
    """Run all nucleotide composition analyses on all genes."""
    results = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_name = record.id
        seq = record.seq

        # Parse category from description
        desc = record.description
        if "ONCOGENE" in desc.upper():
            category = "Oncogene"
        elif "TSG" in desc.upper():
            category = "TSG"
        else:
            category = "Dual-Role"

        # Run all analyses
        nuc_freq = calculate_nucleotide_frequencies(seq)
        gc_overall = calculate_gc_content(seq)
        gc_pos = calculate_positional_gc(seq)
        gc3s = calculate_gc3s(seq)
        skew = calculate_nucleotide_skew(seq)
        cpg_oe = calculate_cpg_oe(seq)
        pr2 = calculate_pr2_bias(seq)

        result = {
            "Gene": gene_name,
            "Category": category,
            "CDS_Length": nuc_freq["length"],
            "A%": round(nuc_freq["A%"], 2),
            "T%": round(nuc_freq["T%"], 2),
            "G%": round(nuc_freq["G%"], 2),
            "C%": round(nuc_freq["C%"], 2),
            "GC%": round(gc_overall, 2),
            "AT%": round(100 - gc_overall, 2),
            "GC1%": round(gc_pos["GC1"], 2),
            "GC2%": round(gc_pos["GC2"], 2),
            "GC3%": round(gc_pos["GC3"], 2),
            "GC12%": round(gc_pos["GC12"], 2),
            "GC3s%": round(gc3s, 2),
            "GC_skew": round(skew["GC_skew"], 4),
            "AT_skew": round(skew["AT_skew"], 4),
            "CpG_OE": round(cpg_oe, 4),
            "A3/(A3+T3)": round(pr2["A3/(A3+T3)"], 4),
            "G3/(G3+C3)": round(pr2["G3/(G3+C3)"], 4),
            "Purine%": round(nuc_freq["A%"] + nuc_freq["G%"], 2),
            "Pyrimidine%": round(nuc_freq["T%"] + nuc_freq["C%"], 2),
        }

        results.append(result)
        print(f"  {gene_name}: GC={gc_overall:.1f}%, GC3={gc_pos['GC3']:.1f}%, "
              f"CpG O/E={cpg_oe:.3f}")

    return results


def save_results(results, output_file):
    """Save results to TSV file."""
    if not results:
        print("No results to save!")
        return

    fieldnames = results[0].keys()
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"\nResults saved to: {output_file}")


def print_summary_statistics(results):
    """Print summary statistics by gene category."""
    categories = {"Oncogene": [], "TSG": [], "Dual-Role": []}

    for r in results:
        cat = r["Category"]
        if cat in categories:
            categories[cat].append(r)

    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)

    for cat, genes in categories.items():
        if not genes:
            continue

        gc_values = [g["GC%"] for g in genes]
        gc3_values = [g["GC3%"] for g in genes]
        cpg_values = [g["CpG_OE"] for g in genes]
        len_values = [g["CDS_Length"] for g in genes]

        print(f"\n{cat} (n={len(genes)}):")
        print(f"  GC%:      {np.mean(gc_values):.1f} ± {np.std(gc_values):.1f} "
              f"(range: {min(gc_values):.1f} - {max(gc_values):.1f})")
        print(f"  GC3%:     {np.mean(gc3_values):.1f} ± {np.std(gc3_values):.1f}")
        print(f"  CpG O/E:  {np.mean(cpg_values):.3f} ± {np.std(cpg_values):.3f}")
        print(f"  CDS Len:  {np.mean(len_values):.0f} ± {np.std(len_values):.0f} bp")


def generate_sliding_window_data(fasta_file, output_dir):
    """Generate sliding window GC data for key genes."""
    key_genes = ["TP53", "EGFR", "BRCA1", "KRAS", "MYC", "APC", "BRAF", "PTEN"]

    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in key_genes:
            sw_data = sliding_window_gc(record.seq, window_size=100, step_size=25)
            outfile = os.path.join(output_dir, f"{record.id}_sliding_window.tsv")
            with open(outfile, "w") as f:
                f.write("Position\tGC_Content\n")
                for pos, gc in sw_data:
                    f.write(f"{pos}\t{gc:.2f}\n")
            print(f"  Sliding window data saved for {record.id}")


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("GC CONTENT & NUCLEOTIDE COMPOSITION ANALYSIS")
    print("=" * 70)

    if not os.path.exists(FASTA_FILE):
        print(f"ERROR: Input file not found: {FASTA_FILE}")
        print("Please run 01_data_acquisition.py first.")
        return

    # Run analysis
    print("\nAnalyzing nucleotide composition...")
    results = analyze_all_genes(FASTA_FILE)

    # Save results
    output_file = os.path.join(OUTPUT_DIR, "gc_content_results.tsv")
    save_results(results, output_file)

    # Print summary
    print_summary_statistics(results)

    # Generate sliding window data
    print("\nGenerating sliding window data for key genes...")
    generate_sliding_window_data(FASTA_FILE, OUTPUT_DIR)

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
