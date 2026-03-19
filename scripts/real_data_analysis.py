#!/usr/bin/env python3
"""
=============================================================================
REAL DATA ANALYSIS — NCBI Sequences
GC Content Variability & Nucleotide Composition in Cancer-Associated Genes
=============================================================================
Downloads REAL coding sequences from NCBI RefSeq, extracts CDS regions,
and performs the complete compositional analysis with biological data.

Author: Sudiksha
Date: March 2026
=============================================================================
"""

import os, sys, time, math, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from collections import Counter
from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction

warnings.filterwarnings('ignore')
Entrez.email = "research@university.edu"

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG_DIR = os.path.join(BASE_DIR, "figures")
DATA_DIR = os.path.join(BASE_DIR, "data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 150, "font.size": 11, "font.family": "serif",
    "axes.labelsize": 13, "axes.titlesize": 14,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 10, "figure.facecolor": "white",
})
COLORS = {"Oncogene": "#E74C3C", "TSG": "#3498DB", "Dual-Role": "#2ECC71"}

# ============================================================
# GENE LIST WITH REFSEQ ACCESSIONS
# ============================================================
GENES = {
    # ONCOGENES
    "KRAS":   {"acc": "NM_004985.5",     "cat": "Oncogene"},
    "BRAF":   {"acc": "NM_004333.6",     "cat": "Oncogene"},
    "MYC":    {"acc": "NM_002467.6",     "cat": "Oncogene"},
    "EGFR":   {"acc": "NM_005228.5",     "cat": "Oncogene"},
    "ERBB2":  {"acc": "NM_004448.4",     "cat": "Oncogene"},
    "PIK3CA": {"acc": "NM_006218.4",     "cat": "Oncogene"},
    "ABL1":   {"acc": "NM_005157.6",     "cat": "Oncogene"},
    "RET":    {"acc": "NM_020975.6",     "cat": "Oncogene"},
    "KIT":    {"acc": "NM_000222.3",     "cat": "Oncogene"},
    "MET":    {"acc": "NM_000245.4",     "cat": "Oncogene"},
    "ALK":    {"acc": "NM_004304.5",     "cat": "Oncogene"},
    "NRAS":   {"acc": "NM_002524.5",     "cat": "Oncogene"},
    "HRAS":   {"acc": "NM_005343.4",     "cat": "Oncogene"},
    "FLT3":   {"acc": "NM_004119.3",     "cat": "Oncogene"},
    "JAK2":   {"acc": "NM_004972.4",     "cat": "Oncogene"},
    "FGFR1":  {"acc": "NM_023110.3",     "cat": "Oncogene"},
    "FGFR2":  {"acc": "NM_000141.5",     "cat": "Oncogene"},
    "FGFR3":  {"acc": "NM_000142.5",     "cat": "Oncogene"},
    "CDK4":   {"acc": "NM_000075.4",     "cat": "Oncogene"},
    "CCND1":  {"acc": "NM_053056.3",     "cat": "Oncogene"},
    # TUMOR SUPPRESSORS
    "TP53":   {"acc": "NM_000546.6",     "cat": "TSG"},
    "BRCA1":  {"acc": "NM_007294.4",     "cat": "TSG"},
    "BRCA2":  {"acc": "NM_000059.4",     "cat": "TSG"},
    "RB1":    {"acc": "NM_000321.3",     "cat": "TSG"},
    "APC":    {"acc": "NM_000038.6",     "cat": "TSG"},
    "PTEN":   {"acc": "NM_000314.8",     "cat": "TSG"},
    "VHL":    {"acc": "NM_000551.4",     "cat": "TSG"},
    "WT1":    {"acc": "NM_024426.6",     "cat": "TSG"},
    "NF1":    {"acc": "NM_001042492.3",  "cat": "TSG"},
    "NF2":    {"acc": "NM_000268.4",     "cat": "TSG"},
    "CDKN2A": {"acc": "NM_000077.5",     "cat": "TSG"},
    "SMAD4":  {"acc": "NM_005359.6",     "cat": "TSG"},
    "STK11":  {"acc": "NM_000455.5",     "cat": "TSG"},
    "MLH1":   {"acc": "NM_000249.4",     "cat": "TSG"},
    "MSH2":   {"acc": "NM_000251.3",     "cat": "TSG"},
    "ATM":    {"acc": "NM_000051.4",     "cat": "TSG"},
    "BAP1":   {"acc": "NM_004656.4",     "cat": "TSG"},
    "CHEK2":  {"acc": "NM_007194.4",     "cat": "TSG"},
    "CDH1":   {"acc": "NM_004360.5",     "cat": "TSG"},
    "FBXW7":  {"acc": "NM_033632.3",     "cat": "TSG"},
    # DUAL-ROLE
    "NOTCH1": {"acc": "NM_017617.5",     "cat": "Dual-Role"},
    "IDH1":   {"acc": "NM_005896.4",     "cat": "Dual-Role"},
    "IDH2":   {"acc": "NM_002168.4",     "cat": "Dual-Role"},
    "EZH2":   {"acc": "NM_004456.5",     "cat": "Dual-Role"},
    "DNMT3A": {"acc": "NM_175629.2",     "cat": "Dual-Role"},
    "TET2":   {"acc": "NM_001127208.3",  "cat": "Dual-Role"},
    "SF3B1":  {"acc": "NM_012433.4",     "cat": "Dual-Role"},
    "NPM1":   {"acc": "NM_002520.7",     "cat": "Dual-Role"},
    "CTNNB1": {"acc": "NM_001904.4",     "cat": "Dual-Role"},
    "MDM2":   {"acc": "NM_002392.6",     "cat": "Dual-Role"},
}

# ============================================================
# GENETIC CODE
# ============================================================
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
AA_CODONS = {}
for c, a in CODON_TABLE.items():
    if a != "*": AA_CODONS.setdefault(a, []).append(c)


# ============================================================
# STEP 1: DOWNLOAD REAL CDS FROM NCBI
# ============================================================
def fetch_cds(gene_name, accession):
    """Fetch CDS from NCBI GenBank record."""
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Extract CDS feature
        for feature in record.features:
            if feature.type == "CDS":
                cds_seq = str(feature.extract(record.seq)).upper()
                # Validate
                if cds_seq[:3] == "ATG" and len(cds_seq) % 3 == 0:
                    return cds_seq
                elif len(cds_seq) % 3 == 0:
                    return cds_seq  # Some CDS may not start with ATG in record

        # Fallback: try fasta_cds_na
        handle2 = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta_cds_na", retmode="text")
        records = list(SeqIO.parse(handle2, "fasta"))
        handle2.close()
        if records:
            seq = str(max(records, key=lambda r: len(r.seq)).seq).upper()
            return seq

        return None
    except Exception as e:
        print(f"    ERROR: {e}")
        return None


# ============================================================
# ANALYSIS FUNCTIONS
# ============================================================
def calc_positional_gc(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq)-len(seq)%3, 3)]
    n = len(codons)
    if n == 0: return 0,0,0,0
    gc1 = sum(1 for c in codons if c[0] in 'GC') / n * 100
    gc2 = sum(1 for c in codons if c[1] in 'GC') / n * 100
    gc3 = sum(1 for c in codons if c[2] in 'GC') / n * 100
    return gc1, gc2, gc3, (gc1+gc2)/2

def calc_gc3s(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq)-len(seq)%3, 3)]
    exclude = {"ATG","TGG","TAA","TAG","TGA"}
    syn = [c for c in codons if c not in exclude and len(c)==3]
    if not syn: return 50.0
    return sum(1 for c in syn if c[2] in 'GC') / len(syn) * 100

def calc_skew(seq):
    cnt = Counter(seq)
    g, c = cnt.get('G',0), cnt.get('C',0)
    a, t = cnt.get('A',0), cnt.get('T',0)
    return ((g-c)/(g+c) if (g+c)>0 else 0, (a-t)/(a+t) if (a+t)>0 else 0)

def calc_cpg_oe(seq):
    n = len(seq)
    cnt = Counter(seq)
    cpg = sum(1 for i in range(n-1) if seq[i:i+2]=='CG')
    c_n, g_n = cnt.get('C',0), cnt.get('G',0)
    if c_n==0 or g_n==0: return 0
    return (cpg/(n-1)) / ((c_n/n)*(g_n/n))

def calc_pr2(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq)-len(seq)%3, 3)]
    third = [c[2] for c in codons if len(c)==3]
    cnt = Counter(third)
    a3,t3 = cnt.get('A',0), cnt.get('T',0)
    g3,c3 = cnt.get('G',0), cnt.get('C',0)
    return (a3/(a3+t3) if (a3+t3)>0 else 0.5, g3/(g3+c3) if (g3+c3)>0 else 0.5)

def calc_rscu(seq):
    cds = seq[:-3] if seq[-3:] in ("TAA","TAG","TGA") else seq
    codons = [cds[i:i+3] for i in range(0, len(cds)-len(cds)%3, 3)]
    counts = Counter(codons)
    rscu = {}
    for aa, aa_codons in AA_CODONS.items():
        total = sum(counts.get(c,0) for c in aa_codons)
        n_syn = len(aa_codons)
        for c in aa_codons:
            exp = total/n_syn if n_syn>0 else 0
            rscu[c] = counts.get(c,0)/exp if exp>0 else 0
    return rscu

def calc_enc(seq):
    cds = seq[:-3] if seq[-3:] in ("TAA","TAG","TGA") else seq
    codons = [cds[i:i+3] for i in range(0, len(cds)-len(cds)%3, 3)]
    counts = Counter(codons)
    deg_groups = {}
    for aa, aa_codons in AA_CODONS.items():
        d = len(aa_codons)
        if d > 1: deg_groups.setdefault(d, []).append(aa)

    def homozygosity(aa):
        aa_c = AA_CODONS[aa]
        cnts = [counts.get(c,0) for c in aa_c]
        n = sum(cnts)
        if n <= 1: return None
        f = sum((x/n)**2 for x in cnts)
        return max((n*f-1)/(n-1), 1/len(aa_c))

    f_vals = {}
    for d, aas in deg_groups.items():
        fl = [homozygosity(a) for a in aas if homozygosity(a) is not None]
        f_vals[d] = np.mean(fl) if fl else 1/d

    enc = 2
    enc += 9/f_vals.get(2, 0.5)
    enc += 1/f_vals.get(3, 0.333)
    enc += 5/f_vals.get(4, 0.25)
    enc += 3/f_vals.get(6, 0.167)
    return max(20, min(61, enc))

def enc_expected(gc3s):
    s = gc3s / 100
    d = s**2 + (1-s)**2
    return 2 + s + 29/d if d>0 else 61

def sliding_window_gc(seq, window=100, step=25):
    results = []
    for i in range(0, len(seq)-window+1, step):
        w = seq[i:i+window]
        gc = (w.count('G')+w.count('C'))/window*100
        results.append((i+window//2, gc))
    return results


# ============================================================
# MAIN EXECUTION
# ============================================================
print("=" * 75)
print("  REAL DATA ANALYSIS — NCBI RefSeq Sequences")
print("  GC Content & Nucleotide Composition in Cancer Genes")
print("=" * 75)

# --- Download sequences ---
print("\nSTEP 1: Downloading CDS from NCBI (50 genes)...")
print("-" * 75)

sequences = {}
failed = []

for i, (gene, info) in enumerate(GENES.items(), 1):
    print(f"  [{i:02d}/50] {gene} ({info['acc']})...", end=" ", flush=True)
    cds = fetch_cds(gene, info["acc"])
    if cds:
        sequences[gene] = cds
        gc = (cds.count('G')+cds.count('C'))/len(cds)*100
        print(f"OK — {len(cds)} bp, GC={gc:.1f}%")
    else:
        failed.append(gene)
        print("FAILED")
    time.sleep(0.35)  # NCBI rate limit

print(f"\nDownloaded: {len(sequences)}/50")
if failed:
    print(f"Failed: {', '.join(failed)}")

# Save FASTA
fasta_file = os.path.join(DATA_DIR, "real_cancer_genes_cds.fasta")
with open(fasta_file, "w") as f:
    for gene, seq in sequences.items():
        f.write(f">{gene} | {GENES[gene]['cat']} | {GENES[gene]['acc']} | {len(seq)}bp\n")
        for j in range(0, len(seq), 70):
            f.write(seq[j:j+70] + "\n")
print(f"Saved: {fasta_file}")

# --- Analyze all genes ---
print(f"\n\nSTEP 2: Computing nucleotide composition on REAL sequences...")
print("-" * 75)

results = []
all_rscu = {}

for gene, seq in sequences.items():
    cat = GENES[gene]["cat"]
    gc_overall = (seq.count('G')+seq.count('C'))/len(seq)*100
    gc1, gc2, gc3, gc12 = calc_positional_gc(seq)
    gc3s = calc_gc3s(seq)
    gc_skew, at_skew = calc_skew(seq)
    cpg_oe = calc_cpg_oe(seq)
    at_bias, gc_bias = calc_pr2(seq)
    rscu = calc_rscu(seq)
    enc_val = calc_enc(seq)
    enc_exp = enc_expected(gc3s)
    all_rscu[gene] = rscu

    cnt = Counter(seq)
    total = len(seq)

    results.append({
        "Gene": gene, "Category": cat, "CDS_Length": len(seq),
        "A%": cnt['A']/total*100, "T%": cnt['T']/total*100,
        "G%": cnt['G']/total*100, "C%": cnt['C']/total*100,
        "GC%": gc_overall, "AT%": 100-gc_overall,
        "GC1%": gc1, "GC2%": gc2, "GC3%": gc3, "GC12%": gc12,
        "GC3s%": gc3s, "GC_skew": gc_skew, "AT_skew": at_skew,
        "CpG_OE": cpg_oe, "A3/(A3+T3)": at_bias, "G3/(G3+C3)": gc_bias,
        "ENC": enc_val, "ENC_expected": enc_exp,
        "Sequence": seq,
    })

df = pd.DataFrame(results)

# Print full gene table
print(f"\n{'Gene':<10} {'Cat':<10} {'Len':>6} {'GC%':>7} {'GC1%':>7} {'GC2%':>7} {'GC3%':>7} {'GC3s%':>7} {'ENC':>6} {'CpG':>7}")
print("-" * 90)
for _, r in df.iterrows():
    print(f"{r['Gene']:<10} {r['Category']:<10} {r['CDS_Length']:>6} {r['GC%']:>7.2f} "
          f"{r['GC1%']:>7.2f} {r['GC2%']:>7.2f} {r['GC3%']:>7.2f} {r['GC3s%']:>7.2f} "
          f"{r['ENC']:>6.1f} {r['CpG_OE']:>7.4f}")

# --- Summary Statistics ---
print("\n\n" + "=" * 75)
print("  STEP 3: SUMMARY STATISTICS BY CATEGORY")
print("=" * 75)

onc = df[df["Category"]=="Oncogene"]
tsg = df[df["Category"]=="TSG"]
dual = df[df["Category"]=="Dual-Role"]

for cat_name, sub in [("Oncogene", onc), ("TSG", tsg), ("Dual-Role", dual)]:
    print(f"\n--- {cat_name} (n={len(sub)}) ---")
    for p in ["GC%","GC1%","GC2%","GC3%","GC3s%","ENC","CpG_OE","CDS_Length"]:
        d = sub[p]
        print(f"  {p:<12} Mean={d.mean():>8.2f}  SD={d.std():>7.2f}  Range=[{d.min():.2f} - {d.max():.2f}]")

# --- Statistical Tests ---
print("\n\n" + "=" * 75)
print("  STEP 4: MANN-WHITNEY U TESTS (Oncogenes vs TSGs)")
print("=" * 75)

print(f"\n{'Parameter':<12} {'Oncogene':>14} {'TSG':>14} {'U-stat':>10} {'p-value':>12} {'Sig':>5}")
print("-" * 70)
for p in ["GC%","GC1%","GC2%","GC3%","GC3s%","ENC","CpG_OE","CDS_Length","A%","T%","G%","C%"]:
    u, pv = stats.mannwhitneyu(onc[p], tsg[p], alternative="two-sided")
    sig = "***" if pv<0.001 else "**" if pv<0.01 else "*" if pv<0.05 else "ns"
    print(f"{p:<12} {onc[p].mean():>14.2f} {tsg[p].mean():>14.2f} {u:>10.1f} {pv:>12.6f} {sig:>5}")

# Kruskal-Wallis
print("\n\nKruskal-Wallis (3 categories):")
for p in ["GC%","GC3%","ENC","CpG_OE"]:
    h, pv = stats.kruskal(onc[p], tsg[p], dual[p])
    sig = "***" if pv<0.001 else "**" if pv<0.01 else "*" if pv<0.05 else "ns"
    print(f"  {p:<12} H={h:>8.3f}  p={pv:.6f}  {sig}")

# --- Correlation ---
print("\n\n" + "=" * 75)
print("  STEP 5: PEARSON CORRELATION ANALYSIS")
print("=" * 75)
pairs = [("GC%","GC3%"),("GC%","ENC"),("GC%","CpG_OE"),("GC3%","GC12%"),
         ("GC%","CDS_Length"),("GC3%","ENC"),("GC_skew","AT_skew")]
print(f"\n{'Pair':<25} {'r':>8} {'r²':>8} {'p-value':>12} {'Strength':<20}")
print("-" * 75)
for p1,p2 in pairs:
    r, pv = stats.pearsonr(df[p1], df[p2])
    s = "Strong" if abs(r)>=0.7 else "Moderate" if abs(r)>=0.4 else "Weak"
    d = "positive" if r>0 else "negative"
    print(f"{p1+' vs '+p2:<25} {r:>8.4f} {r**2:>8.4f} {pv:>12.2e} {s+' '+d:<20}")

# --- Neutrality Plot ---
print("\n\n" + "=" * 75)
print("  STEP 6: NEUTRALITY PLOT REGRESSION (GC12 vs GC3)")
print("=" * 75)
slope, intercept, r_val, p_val, se = stats.linregress(df["GC3%"], df["GC12%"])
print(f"\n  ALL GENES (n={len(df)}):")
print(f"  Equation: GC12 = {slope:.4f} × GC3 + {intercept:.2f}")
print(f"  R² = {r_val**2:.4f}, p = {p_val:.2e}")
print(f"  ┌───────────────────────────────────────────────┐")
print(f"  │  Mutation pressure contribution:   {slope*100:>5.1f}%      │")
print(f"  │  Natural selection contribution:   {(1-slope)*100:>5.1f}%      │")
print(f"  └───────────────────────────────────────────────┘")

for cat_name, sub in [("Oncogene", onc), ("TSG", tsg), ("Dual-Role", dual)]:
    if len(sub) >= 3:
        s, i, r, p, _ = stats.linregress(sub["GC3%"], sub["GC12%"])
        print(f"  {cat_name} (n={len(sub)}): slope={s:.4f}, R²={r**2:.4f}, "
              f"mutation={s*100:.1f}%, selection={(1-s)*100:.1f}%")

# --- ENC-GC3 Analysis ---
print("\n\n" + "=" * 75)
print("  STEP 7: ENC-GC3s PLOT ANALYSIS")
print("=" * 75)
below = sum(1 for _,r in df.iterrows() if r["ENC"] < enc_expected(r["GC3s%"]) * 0.95)
on = len(df) - below
print(f"  Below expected curve (selection): {below}/{len(df)} ({below/len(df)*100:.0f}%)")
print(f"  On/near curve (mutation only):    {on}/{len(df)} ({on/len(df)*100:.0f}%)")

# --- RSCU ---
print("\n\n" + "=" * 75)
print("  STEP 8: CODON USAGE — RSCU ANALYSIS")
print("=" * 75)
sense = [c for c, a in CODON_TABLE.items() if a not in ("*","M","W")]
mean_rscu = {}
for codon in sense:
    vals = [all_rscu[g].get(codon, 0) for g in sequences.keys()]
    mean_rscu[codon] = np.mean(vals)

preferred = sorted([(c,mean_rscu[c],CODON_TABLE[c]) for c in sense if mean_rscu[c]>1.6],
                   key=lambda x: x[1], reverse=True)
avoided = sorted([(c,mean_rscu[c],CODON_TABLE[c]) for c in sense if mean_rscu[c]<0.6],
                 key=lambda x: x[1])

print(f"\nOVERREPRESENTED (RSCU > 1.6):")
print(f"  {'Codon':<8} {'AA':<5} {'RSCU':>8} {'Ending':<10}")
print(f"  {'-'*33}")
for c,r,a in preferred:
    print(f"  {c:<8} {a:<5} {r:>8.3f} {'GC-end' if c[2] in 'GC' else 'AT-end':<10}")

print(f"\nUNDERREPRESENTED (RSCU < 0.6):")
print(f"  {'Codon':<8} {'AA':<5} {'RSCU':>8} {'Ending':<10}")
print(f"  {'-'*33}")
for c,r,a in avoided:
    print(f"  {c:<8} {a:<5} {r:>8.3f} {'GC-end' if c[2] in 'GC' else 'AT-end':<10}")

# --- CpG ---
print("\n\n" + "=" * 75)
print("  STEP 9: CpG ISLAND ANALYSIS")
print("=" * 75)
df["CpG_Island"] = (df["CpG_OE"] >= 0.6) & (df["GC%"] >= 50)
onc = df[df["Category"]=="Oncogene"]
tsg = df[df["Category"]=="TSG"]
for cat_name in ["Oncogene","TSG","Dual-Role"]:
    sub = df[df["Category"]==cat_name]
    n = sub["CpG_Island"].sum()
    print(f"  {cat_name:<12}: {n}/{len(sub)} ({n/len(sub)*100:.0f}%) meet CpG island criteria")

a = onc["CpG_Island"].sum(); b = len(onc)-a
c = tsg["CpG_Island"].sum(); d = len(tsg)-c
odds, pf = stats.fisher_exact([[a,b],[c,d]])
print(f"  Fisher's exact (Onc vs TSG): OR={odds:.3f}, p={pf:.6f}")

# --- PR2 ---
print("\n\n" + "=" * 75)
print("  STEP 10: PR2 BIAS ANALYSIS")
print("=" * 75)
print(f"  Mean A3/(A3+T3) = {df['A3/(A3+T3)'].mean():.4f} ± {df['A3/(A3+T3)'].std():.4f}")
print(f"  Mean G3/(G3+C3) = {df['G3/(G3+C3)'].mean():.4f} ± {df['G3/(G3+C3)'].std():.4f}")
dev = sum(1 for _,r in df.iterrows() if abs(r['A3/(A3+T3)']-0.5)>0.02 or abs(r['G3/(G3+C3)']-0.5)>0.02)
print(f"  Genes deviating from parity: {dev}/{len(df)} ({dev/len(df)*100:.0f}%)")


# ==============================================================
# FIGURE GENERATION
# ==============================================================
print("\n\n" + "=" * 75)
print("  GENERATING FIGURES WITH REAL DATA...")
print("=" * 75)

# Fig 1: GC content bars
fig, ax = plt.subplots(figsize=(16,7))
df_s = df.sort_values("GC%").reset_index(drop=True)
ax.bar(range(len(df_s)), df_s["GC%"], color=[COLORS[c] for c in df_s["Category"]],
       edgecolor="black", lw=0.4, alpha=0.85)
ax.axhline(41, color="gray", ls="--", lw=1.5, label="Human genome avg (41%)")
ax.axhline(50, color="black", ls=":", lw=0.8, alpha=0.4)
ax.set_xticks(range(len(df_s))); ax.set_xticklabels(df_s["Gene"], rotation=90, fontsize=7)
ax.set_ylabel("GC Content (%)"); ax.set_ylim(30,75)
ax.set_title("Figure 1: GC Content of 50 Cancer Genes (REAL NCBI Data)", fontweight="bold")
lp = [mpatches.Patch(color=COLORS[c], label=f"{c} (n={len(df[df['Category']==c])})") for c in COLORS]
lp.append(plt.Line2D([0],[0], color="gray", ls="--", label="Genome avg (41%)"))
ax.legend(handles=lp, loc="upper left")
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig1_gc_content_bars.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig1_gc_content_bars.png")

# Fig 2: Box plots
fig, axes = plt.subplots(1,4, figsize=(18,5))
for ax, m, t in zip(axes, ["GC%","GC1%","GC2%","GC3%"], ["Overall GC","GC 1st Pos","GC 2nd Pos","GC 3rd Pos"]):
    order = ["Oncogene","TSG","Dual-Role"]; pal = [COLORS[c] for c in order]
    sns.boxplot(data=df, x="Category", y=m, order=order, palette=pal, ax=ax, width=0.5)
    sns.stripplot(data=df, x="Category", y=m, order=order, color="black", size=4, alpha=0.3, ax=ax, jitter=True)
    _, pv = stats.mannwhitneyu(onc[m], tsg[m]); sig = "***" if pv<0.001 else "**" if pv<0.01 else "*" if pv<0.05 else "ns"
    ax.set_title(f"{t}\n(p={pv:.4f} {sig})", fontsize=11, fontweight="bold"); ax.set_xlabel("")
plt.suptitle("Figure 2: GC Content by Category (REAL Data)", fontsize=14, fontweight="bold", y=1.03)
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig2_gc_boxplots.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig2_gc_boxplots.png")

# Fig 3: ENC-GC3s
fig, ax = plt.subplots(figsize=(10,7))
gc3r = np.linspace(20,85,200)
ax.plot(gc3r, [enc_expected(g) for g in gc3r], "k-", lw=2, label="Expected (mutation only)")
for cat in COLORS:
    sub = df[df["Category"]==cat]
    ax.scatter(sub["GC3s%"], sub["ENC"], c=COLORS[cat], s=70, edgecolors="black", lw=0.5, label=f"{cat}", alpha=0.85, zorder=5)
    for _,r in sub.iterrows():
        ax.annotate(r["Gene"], (r["GC3s%"],r["ENC"]), fontsize=5.5, ha="center", va="bottom", alpha=0.65)
ax.set_xlabel("GC3s (%)"); ax.set_ylabel("ENC"); ax.set_xlim(20,85); ax.set_ylim(25,65)
ax.set_title("Figure 3: ENC-GC3s Plot (REAL Data)", fontweight="bold")
ax.legend(loc="lower left"); ax.grid(True, alpha=0.2)
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig3_enc_gc3_plot.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig3_enc_gc3_plot.png")

# Fig 4: Neutrality plot
fig, ax = plt.subplots(figsize=(9,7))
for cat in COLORS:
    sub = df[df["Category"]==cat]
    ax.scatter(sub["GC3%"], sub["GC12%"], c=COLORS[cat], s=70, edgecolors="black", lw=0.5, label=cat, alpha=0.85, zorder=5)
sl,ic,rv,pv,_ = stats.linregress(df["GC3%"], df["GC12%"])
x_ln = np.linspace(df["GC3%"].min()-5, df["GC3%"].max()+5, 100)
ax.plot(x_ln, sl*x_ln+ic, "r--", lw=2, label=f"Regression (slope={sl:.3f})")
ax.plot([20,85],[20,85], "k:", lw=1, alpha=0.4, label="Complete neutrality")
ax.set_xlabel("GC3 (%)"); ax.set_ylabel("GC12 (%)")
ax.set_title("Figure 4: Neutrality Plot (REAL Data)", fontweight="bold")
ax.legend(loc="upper left"); ax.grid(True, alpha=0.2)
txt = f"Slope = {sl:.3f} (R² = {rv**2:.3f})\nMutation: {sl*100:.1f}%\nSelection: {(1-sl)*100:.1f}%"
ax.text(0.95, 0.05, txt, transform=ax.transAxes, fontsize=10, ha="right", va="bottom",
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig4_neutrality_plot.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig4_neutrality_plot.png")

# Fig 5: PR2 bias
fig, ax = plt.subplots(figsize=(8,8))
for cat in COLORS:
    sub = df[df["Category"]==cat]
    ax.scatter(sub["G3/(G3+C3)"], sub["A3/(A3+T3)"], c=COLORS[cat], s=70, edgecolors="black", lw=0.5, label=cat, alpha=0.85, zorder=5)
ax.axhline(0.5, color="gray", ls="--", lw=1, alpha=0.6); ax.axvline(0.5, color="gray", ls="--", lw=1, alpha=0.6)
ax.plot(0.5, 0.5, "k+", ms=15, mew=2, zorder=10)
ax.set_xlabel("G3/(G3+C3)"); ax.set_ylabel("A3/(A3+T3)")
ax.set_title("Figure 5: PR2 Bias Plot (REAL Data)", fontweight="bold")
ax.set_xlim(0.3,0.7); ax.set_ylim(0.3,0.7); ax.set_aspect("equal")
ax.legend(loc="upper right"); ax.grid(True, alpha=0.15)
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig5_pr2_bias.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig5_pr2_bias.png")

# Fig 6: CpG
fig, (ax1,ax2) = plt.subplots(1,2, figsize=(14,6))
df_cs = df.sort_values("CpG_OE").reset_index(drop=True)
ax1.barh(range(len(df_cs)), df_cs["CpG_OE"], color=[COLORS[c] for c in df_cs["Category"]], edgecolor="black", lw=0.3, alpha=0.85)
ax1.axvline(0.6, color="red", ls="--", lw=1.5, label="CpG island threshold"); ax1.set_yticks(range(len(df_cs)))
ax1.set_yticklabels(df_cs["Gene"], fontsize=6); ax1.set_xlabel("CpG O/E"); ax1.set_title("CpG O/E by Gene", fontweight="bold"); ax1.legend(fontsize=8)
order = ["Oncogene","TSG","Dual-Role"]
sns.boxplot(data=df, x="Category", y="CpG_OE", order=order, palette=[COLORS[c] for c in order], ax=ax2, width=0.5)
sns.stripplot(data=df, x="Category", y="CpG_OE", order=order, color="black", size=5, alpha=0.3, ax=ax2, jitter=True)
ax2.axhline(0.6, color="red", ls="--", lw=1.5, alpha=0.6); ax2.set_title("CpG O/E by Category", fontweight="bold")
plt.suptitle("Figure 6: CpG Analysis (REAL Data)", fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig6_cpg_analysis.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig6_cpg_analysis.png")

# Fig 7: Nucleotide composition
fig, ax = plt.subplots(figsize=(16,6))
df_ns = df.sort_values("GC%").reset_index(drop=True); x = range(len(df_ns))
ax.bar(x, df_ns["A%"], label="A", color="#2196F3", alpha=0.85)
ax.bar(x, df_ns["T%"], bottom=df_ns["A%"], label="T", color="#FF9800", alpha=0.85)
ax.bar(x, df_ns["G%"], bottom=df_ns["A%"]+df_ns["T%"], label="G", color="#4CAF50", alpha=0.85)
ax.bar(x, df_ns["C%"], bottom=df_ns["A%"]+df_ns["T%"]+df_ns["G%"], label="C", color="#F44336", alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(df_ns["Gene"], rotation=90, fontsize=7)
ax.set_ylabel("Nucleotide (%)"); ax.set_title("Figure 7: Nucleotide Composition (REAL Data)", fontweight="bold")
ax.legend(ncol=4, loc="upper left"); ax.set_ylim(0,105)
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig7_nucleotide_composition.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig7_nucleotide_composition.png")

# Fig 8: Correlation matrix
fig, ax = plt.subplots(figsize=(10,8))
cols = ["GC%","GC1%","GC2%","GC3%","ENC","CpG_OE","CDS_Length","GC_skew","AT_skew"]
corr = df[cols].corr()
mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
sns.heatmap(corr, mask=mask, annot=True, fmt=".2f", cmap="RdBu_r", center=0, vmin=-1, vmax=1,
            square=True, lw=1, ax=ax, annot_kws={"size":9})
ax.set_title("Figure 8: Correlation Matrix (REAL Data)", fontweight="bold", pad=15)
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig8_correlation_matrix.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig8_correlation_matrix.png")

# Fig 9: GC vs Length
fig, ax = plt.subplots(figsize=(9,7))
for cat in COLORS:
    sub = df[df["Category"]==cat]
    ax.scatter(sub["CDS_Length"], sub["GC%"], c=COLORS[cat], s=70, edgecolors="black", lw=0.5, label=cat, alpha=0.85, zorder=5)
    for _,r in sub.iterrows():
        ax.annotate(r["Gene"], (r["CDS_Length"],r["GC%"]), fontsize=5.5, ha="left", va="bottom", alpha=0.6)
sl2,ic2,rv2,pv2,_ = stats.linregress(df["CDS_Length"], df["GC%"])
x2 = np.linspace(0, df["CDS_Length"].max()*1.05, 100)
ax.plot(x2, sl2*x2+ic2, "k--", lw=1.5, alpha=0.5)
ax.set_xlabel("CDS Length (bp)"); ax.set_ylabel("GC (%)"); ax.set_title("Figure 9: GC vs Length (REAL Data)", fontweight="bold")
ax.legend(); ax.grid(True, alpha=0.2)
ax.text(0.95, 0.95, f"r = {rv2:.3f}\np = {pv2:.4f}", transform=ax.transAxes, fontsize=10, ha="right", va="top",
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig9_gc_vs_length.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig9_gc_vs_length.png")

# Fig 10: Sliding window
fig, axes = plt.subplots(2,3, figsize=(16,9))
key_genes = ["TP53","EGFR","BRCA1","KRAS","MYC","APC"]
for ax, gene in zip(axes.flat, key_genes):
    if gene in sequences:
        seq = sequences[gene]; cat = GENES[gene]["cat"]
        gc_avg = (seq.count('G')+seq.count('C'))/len(seq)*100
        sw = sliding_window_gc(seq, window=min(100, len(seq)//5), step=max(10, len(seq)//100))
        if sw:
            pos, gv = zip(*sw)
            ax.plot(pos, gv, color=COLORS[cat], lw=1.2, alpha=0.8)
            ax.fill_between(pos, gv, alpha=0.12, color=COLORS[cat])
            ax.axhline(gc_avg, color="black", ls="--", lw=1, alpha=0.5)
        ax.set_title(f"{gene} ({cat})\nGC={gc_avg:.1f}%", fontweight="bold", fontsize=10)
        ax.set_xlabel("Position (bp)", fontsize=8); ax.set_ylabel("GC%", fontsize=8)
        ax.set_ylim(15,85); ax.grid(True, alpha=0.2)
plt.suptitle("Figure 10: Sliding Window GC Content (REAL Data)", fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout(); plt.savefig(os.path.join(FIG_DIR,"fig10_sliding_window.png"), dpi=200, bbox_inches="tight"); plt.close()
print("  [SAVED] fig10_sliding_window.png")

# --- Save data ---
df_save = df.drop(columns=["Sequence"])
df_save.to_csv(os.path.join(DATA_DIR, "real_gc_content_results.tsv"), sep="\t", index=False)
print(f"\n  [SAVED] data/real_gc_content_results.tsv")

# === FINAL SUMMARY ===
print("\n\n" + "=" * 75)
print("  REAL DATA ANALYSIS — FINAL SUMMARY")
print("=" * 75)
print(f"""
  Genes analyzed:           {len(df)} (from NCBI RefSeq)
  Oncogenes:                {len(onc)} (mean GC = {onc['GC%'].mean():.2f} ± {onc['GC%'].std():.2f}%)
  Tumor Suppressors:        {len(tsg)} (mean GC = {tsg['GC%'].mean():.2f} ± {tsg['GC%'].std():.2f}%)
  Dual-Role:                {len(dual)} (mean GC = {dual['GC%'].mean():.2f} ± {dual['GC%'].std():.2f}%)

  GC% Onc vs TSG:           p = {stats.mannwhitneyu(onc['GC%'],tsg['GC%']).pvalue:.6f}
  GC3% Onc vs TSG:          p = {stats.mannwhitneyu(onc['GC3%'],tsg['GC3%']).pvalue:.6f}
  ENC Onc vs TSG:            p = {stats.mannwhitneyu(onc['ENC'],tsg['ENC']).pvalue:.6f}

  Neutrality slope:          {slope:.4f}
  Mutation pressure:         {slope*100:.1f}%
  Natural selection:         {(1-slope)*100:.1f}%

  Preferred codons (>1.6):   {len(preferred)}
  Avoided codons (<0.6):     {len(avoided)}
  Genes below ENC curve:     {below}/{len(df)} ({below/len(df)*100:.0f}%)

  Figures generated:         10 (all from REAL data)
  Data saved:                real_gc_content_results.tsv
""")
print("=" * 75)
