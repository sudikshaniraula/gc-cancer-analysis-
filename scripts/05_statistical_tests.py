"""
Script 05: Statistical Analysis
=================================
Performs all statistical tests comparing compositional features
between oncogenes and tumor suppressor genes.

Tests performed:
- Shapiro-Wilk normality test
- Mann-Whitney U test (Oncogenes vs TSGs)
- Kruskal-Wallis test (all 3 categories)
- Pearson correlation analysis
- Fisher's exact test (CpG island classification)
- Linear regression (Neutrality plot)

Author: Sudiksha
Project: GC Content Variability in Cancer Genes
Date: April 2025 - March 2026
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations

# ============================================================
# PATHS
# ============================================================
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_FILE = os.path.join(DATA_DIR, "statistical_results.txt")


def load_data():
    """Load analysis results."""
    gc_file = os.path.join(DATA_DIR, "gc_content_results.tsv")
    if os.path.exists(gc_file):
        return pd.read_csv(gc_file, sep="\t")
    else:
        print("WARNING: Data file not found. Run analysis scripts first.")
        return None


def normality_tests(df, output):
    """Shapiro-Wilk normality test for key parameters."""
    output.append("\n" + "=" * 70)
    output.append("1. NORMALITY TESTS (Shapiro-Wilk)")
    output.append("=" * 70)
    output.append(f"{'Parameter':<15} {'Category':<12} {'W-stat':>10} {'p-value':>12} {'Normal?':>10}")
    output.append("-" * 60)

    parameters = ["GC%", "GC3%", "CpG_OE"]
    categories = ["Oncogene", "TSG", "Dual-Role"]

    all_normal = True
    for param in parameters:
        for cat in categories:
            data = df[df["Category"] == cat][param].dropna()
            if len(data) >= 3:
                w_stat, p_value = stats.shapiro(data)
                is_normal = "Yes" if p_value > 0.05 else "No"
                if p_value <= 0.05:
                    all_normal = False
                output.append(f"{param:<15} {cat:<12} {w_stat:>10.4f} {p_value:>12.6f} {is_normal:>10}")

    output.append(f"\nConclusion: {'All data normally distributed' if all_normal else 'Some parameters non-normal → use non-parametric tests'}")


def mann_whitney_tests(df, output):
    """Mann-Whitney U test comparing Oncogenes vs TSGs."""
    output.append("\n" + "=" * 70)
    output.append("2. MANN-WHITNEY U TEST (Oncogenes vs Tumor Suppressor Genes)")
    output.append("=" * 70)

    onc = df[df["Category"] == "Oncogene"]
    tsg = df[df["Category"] == "TSG"]

    parameters = ["GC%", "GC1%", "GC2%", "GC3%", "GC3s%", "CpG_OE",
                  "A%", "T%", "G%", "C%", "GC_skew", "AT_skew",
                  "CDS_Length", "Purine%", "Pyrimidine%"]

    output.append(f"\n{'Parameter':<15} {'Onc Mean±SD':<20} {'TSG Mean±SD':<20} {'U-stat':>10} {'p-value':>12} {'Sig':>6}")
    output.append("-" * 85)

    significant_count = 0
    bonferroni_alpha = 0.05 / len(parameters)

    for param in parameters:
        onc_data = onc[param].dropna()
        tsg_data = tsg[param].dropna()

        if len(onc_data) > 0 and len(tsg_data) > 0:
            u_stat, p_value = stats.mannwhitneyu(onc_data, tsg_data, alternative="two-sided")

            onc_str = f"{onc_data.mean():.2f}±{onc_data.std():.2f}"
            tsg_str = f"{tsg_data.mean():.2f}±{tsg_data.std():.2f}"

            if p_value < 0.001:
                sig = "***"
            elif p_value < 0.01:
                sig = "**"
            elif p_value < 0.05:
                sig = "*"
            else:
                sig = "ns"

            if p_value < 0.05:
                significant_count += 1

            output.append(f"{param:<15} {onc_str:<20} {tsg_str:<20} {u_stat:>10.1f} {p_value:>12.6f} {sig:>6}")

    output.append(f"\nSignificance levels: * p<0.05, ** p<0.01, *** p<0.001")
    output.append(f"Bonferroni-corrected alpha: {bonferroni_alpha:.4f}")
    output.append(f"Significant parameters: {significant_count}/{len(parameters)}")


def kruskal_wallis_test(df, output):
    """Kruskal-Wallis test comparing all three gene categories."""
    output.append("\n" + "=" * 70)
    output.append("3. KRUSKAL-WALLIS TEST (All Three Categories)")
    output.append("=" * 70)

    parameters = ["GC%", "GC3%", "CpG_OE", "CDS_Length"]

    output.append(f"\n{'Parameter':<15} {'H-stat':>10} {'p-value':>12} {'Sig':>6}")
    output.append("-" * 45)

    for param in parameters:
        groups = [df[df["Category"] == cat][param].dropna() for cat in ["Oncogene", "TSG", "Dual-Role"]]
        groups = [g for g in groups if len(g) > 0]

        if len(groups) >= 2:
            h_stat, p_value = stats.kruskal(*groups)
            sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
            output.append(f"{param:<15} {h_stat:>10.4f} {p_value:>12.6f} {sig:>6}")


def correlation_analysis(df, output):
    """Pearson correlation analysis between compositional parameters."""
    output.append("\n" + "=" * 70)
    output.append("4. PEARSON CORRELATION ANALYSIS")
    output.append("=" * 70)

    param_pairs = [
        ("GC%", "GC3%"),
        ("GC%", "CDS_Length"),
        ("GC%", "CpG_OE"),
        ("GC3%", "GC12%"),
        ("GC3%", "GC1%"),
        ("GC3%", "GC2%"),
        ("GC_skew", "AT_skew"),
        ("CDS_Length", "CpG_OE"),
        ("A%", "T%"),
        ("G%", "C%"),
    ]

    output.append(f"\n{'Parameter Pair':<25} {'r':>8} {'r²':>8} {'p-value':>12} {'Sig':>6} {'Interpretation':>20}")
    output.append("-" * 85)

    for p1, p2 in param_pairs:
        if p1 in df.columns and p2 in df.columns:
            data1 = df[p1].dropna()
            data2 = df[p2].dropna()
            # Align indices
            common = data1.index.intersection(data2.index)
            if len(common) >= 3:
                r, p_value = stats.pearsonr(data1[common], data2[common])
                sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"

                if abs(r) >= 0.7:
                    interp = "Strong"
                elif abs(r) >= 0.4:
                    interp = "Moderate"
                elif abs(r) >= 0.2:
                    interp = "Weak"
                else:
                    interp = "Negligible"

                direction = "positive" if r > 0 else "negative"

                output.append(f"{p1+' vs '+p2:<25} {r:>8.4f} {r**2:>8.4f} {p_value:>12.6f} {sig:>6} {interp+' '+direction:>20}")


def neutrality_plot_regression(df, output):
    """Linear regression for neutrality plot (GC12 vs GC3)."""
    output.append("\n" + "=" * 70)
    output.append("5. NEUTRALITY PLOT REGRESSION ANALYSIS")
    output.append("=" * 70)

    # All genes
    slope, intercept, r_value, p_value, std_err = stats.linregress(df["GC3%"], df["GC12%"])
    output.append(f"\nAll genes (n={len(df)}):")
    output.append(f"  Regression: GC12 = {slope:.4f} × GC3 + {intercept:.2f}")
    output.append(f"  R² = {r_value**2:.4f}, p = {p_value:.6f}")
    output.append(f"  Mutation pressure contribution: {slope*100:.1f}%")
    output.append(f"  Natural selection contribution: {(1-slope)*100:.1f}%")

    # By category
    for cat in ["Oncogene", "TSG", "Dual-Role"]:
        subset = df[df["Category"] == cat]
        if len(subset) >= 3:
            s, i, r, p, se = stats.linregress(subset["GC3%"], subset["GC12%"])
            output.append(f"\n{cat} (n={len(subset)}):")
            output.append(f"  Regression: GC12 = {s:.4f} × GC3 + {i:.2f}")
            output.append(f"  R² = {r**2:.4f}, p = {p:.6f}")
            output.append(f"  Mutation pressure: {s*100:.1f}% | Selection: {(1-s)*100:.1f}%")


def fisher_cpg_test(df, output):
    """Fisher's exact test for CpG island classification."""
    output.append("\n" + "=" * 70)
    output.append("6. FISHER'S EXACT TEST (CpG Island Classification)")
    output.append("=" * 70)

    # CpG island criteria: CpG O/E >= 0.6 AND GC% >= 50%
    df["CpG_Island"] = (df["CpG_OE"] >= 0.6) & (df["GC%"] >= 50)

    onc = df[df["Category"] == "Oncogene"]
    tsg = df[df["Category"] == "TSG"]

    # Contingency table
    a = onc["CpG_Island"].sum()  # Oncogene + CpG island
    b = len(onc) - a              # Oncogene + no CpG island
    c = tsg["CpG_Island"].sum()  # TSG + CpG island
    d = len(tsg) - c              # TSG + no CpG island

    table = np.array([[a, b], [c, d]])
    odds_ratio, p_value = stats.fisher_exact(table)

    output.append(f"\nContingency Table:")
    output.append(f"{'':>20} {'CpG Island':>12} {'Non-CpG':>12} {'Total':>8}")
    output.append(f"{'Oncogenes':>20} {int(a):>12} {int(b):>12} {int(a+b):>8}")
    output.append(f"{'TSGs':>20} {int(c):>12} {int(d):>12} {int(c+d):>8}")
    output.append(f"\nOdds Ratio: {odds_ratio:.4f}")
    output.append(f"p-value (Fisher's exact): {p_value:.6f}")
    sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
    output.append(f"Significance: {sig}")
    output.append(f"\nConclusion: {'Significant difference in CpG island frequency between oncogenes and TSGs' if p_value < 0.05 else 'No significant difference'}")


def descriptive_statistics(df, output):
    """Comprehensive descriptive statistics."""
    output.append("\n" + "=" * 70)
    output.append("7. DESCRIPTIVE STATISTICS BY CATEGORY")
    output.append("=" * 70)

    parameters = ["GC%", "GC1%", "GC2%", "GC3%", "GC3s%",
                  "A%", "T%", "G%", "C%",
                  "CpG_OE", "GC_skew", "AT_skew", "CDS_Length"]

    for cat in ["Oncogene", "TSG", "Dual-Role", "All"]:
        if cat == "All":
            subset = df
        else:
            subset = df[df["Category"] == cat]

        output.append(f"\n--- {cat} (n={len(subset)}) ---")
        output.append(f"{'Parameter':<15} {'Mean':>10} {'SD':>10} {'Min':>10} {'Max':>10} {'Median':>10}")
        output.append("-" * 65)

        for param in parameters:
            if param in subset.columns:
                data = subset[param].dropna()
                if len(data) > 0:
                    output.append(f"{param:<15} {data.mean():>10.2f} {data.std():>10.2f} "
                                  f"{data.min():>10.2f} {data.max():>10.2f} {data.median():>10.2f}")


def main():
    """Run all statistical analyses."""
    print("=" * 70)
    print("STATISTICAL ANALYSIS")
    print("=" * 70)

    df = load_data()
    if df is None:
        return

    output = []
    output.append("STATISTICAL ANALYSIS REPORT")
    output.append(f"Generated for: GC Content Variability in Cancer Genes")
    output.append(f"Number of genes: {len(df)}")
    output.append(f"Categories: {df['Category'].value_counts().to_dict()}")

    # Run all tests
    descriptive_statistics(df, output)
    normality_tests(df, output)
    mann_whitney_tests(df, output)
    kruskal_wallis_test(df, output)
    correlation_analysis(df, output)
    neutrality_plot_regression(df, output)
    fisher_cpg_test(df, output)

    # Save report
    report_text = "\n".join(output)
    with open(OUTPUT_FILE, "w") as f:
        f.write(report_text)

    print(report_text)
    print(f"\nFull report saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
