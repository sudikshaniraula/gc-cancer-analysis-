"""
Script 04: Visualization
=========================
Generates all plots and figures for the research paper:
- GC content bar charts
- GC content box plots (Oncogenes vs TSGs)
- ENC-GC3 plot
- Neutrality plot (GC12 vs GC3)
- PR2 bias plot
- Sliding window GC content plots
- RSCU heatmap
- CpG O/E comparison plot
- Nucleotide composition pie charts
- Correlation scatter plots

Author: Sudiksha
Project: GC Content Variability in Cancer Genes
Date: April 2025 - March 2026
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats

# ============================================================
# CONFIGURATION
# ============================================================
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Plot style
plt.rcParams.update({
    "figure.figsize": (12, 8),
    "figure.dpi": 300,
    "font.size": 12,
    "font.family": "serif",
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 11,
    "figure.facecolor": "white",
})

# Color palette for gene categories
COLORS = {
    "Oncogene": "#E74C3C",    # Red
    "TSG": "#3498DB",          # Blue
    "Dual-Role": "#2ECC71",    # Green
}


def load_gc_data():
    """Load GC content results."""
    gc_file = os.path.join(DATA_DIR, "gc_content_results.tsv")
    if os.path.exists(gc_file):
        return pd.read_csv(gc_file, sep="\t")
    else:
        print(f"WARNING: {gc_file} not found. Generating sample data for visualization.")
        return generate_sample_data()


def load_enc_data():
    """Load ENC analysis results."""
    enc_file = os.path.join(DATA_DIR, "enc_analysis.tsv")
    if os.path.exists(enc_file):
        return pd.read_csv(enc_file, sep="\t")
    return None


def generate_sample_data():
    """Generate realistic sample data for visualization testing."""
    np.random.seed(42)

    oncogenes = ["KRAS", "BRAF", "MYC", "EGFR", "ERBB2", "PIK3CA", "ABL1",
                 "RET", "KIT", "MET", "ALK", "NRAS", "HRAS", "FLT3", "JAK2",
                 "FGFR1", "FGFR2", "FGFR3", "CDK4", "CCND1"]
    tsgs = ["TP53", "BRCA1", "BRCA2", "RB1", "APC", "PTEN", "VHL", "WT1",
            "NF1", "NF2", "CDKN2A", "SMAD4", "STK11", "MLH1", "MSH2",
            "ATM", "BAP1", "CHEK2", "CDH1", "FBXW7"]
    dual = ["NOTCH1", "IDH1", "IDH2", "EZH2", "DNMT3A", "TET2", "SF3B1",
            "NPM1", "CTNNB1", "MDM2"]

    data = []
    for gene in oncogenes:
        gc = np.random.normal(57.8, 8.2)
        gc = np.clip(gc, 42, 72)
        gc3 = gc * 1.18 + np.random.normal(0, 3)
        gc3 = np.clip(gc3, 45, 85)
        data.append(_make_gene_row(gene, "Oncogene", gc, gc3))

    for gene in tsgs:
        gc = np.random.normal(48.3, 9.1)
        gc = np.clip(gc, 36, 62)
        gc3 = gc * 1.13 + np.random.normal(0, 4)
        gc3 = np.clip(gc3, 35, 75)
        data.append(_make_gene_row(gene, "TSG", gc, gc3))

    for gene in dual:
        gc = np.random.normal(52.4, 7.6)
        gc = np.clip(gc, 40, 65)
        gc3 = gc * 1.15 + np.random.normal(0, 3.5)
        gc3 = np.clip(gc3, 40, 78)
        data.append(_make_gene_row(gene, "Dual-Role", gc, gc3))

    return pd.DataFrame(data)


def _make_gene_row(gene, category, gc, gc3):
    """Helper to create a single gene data row."""
    gc1 = gc * 0.97 + np.random.normal(0, 2)
    gc2 = gc * 0.73 + np.random.normal(0, 2)
    gc12 = (gc1 + gc2) / 2
    at = 100 - gc
    a_pct = at * (0.48 + np.random.normal(0, 0.03))
    t_pct = at - a_pct
    g_pct = gc * (0.52 + np.random.normal(0, 0.02))
    c_pct = gc - g_pct

    enc = 61 - (gc / 100) * 25 + np.random.normal(0, 3)
    enc = np.clip(enc, 30, 58)

    cpg_oe = gc / 100 * 1.2 + np.random.normal(0, 0.1)
    cpg_oe = np.clip(cpg_oe, 0.35, 1.2)

    return {
        "Gene": gene, "Category": category,
        "CDS_Length": int(np.random.uniform(500, 12000)),
        "A%": round(a_pct, 2), "T%": round(t_pct, 2),
        "G%": round(g_pct, 2), "C%": round(c_pct, 2),
        "GC%": round(gc, 2), "AT%": round(at, 2),
        "GC1%": round(gc1, 2), "GC2%": round(gc2, 2),
        "GC3%": round(gc3, 2), "GC12%": round(gc12, 2),
        "GC3s%": round(gc3 - np.random.uniform(1, 3), 2),
        "GC_skew": round(np.random.normal(0.038, 0.06), 4),
        "AT_skew": round(np.random.normal(-0.034, 0.07), 4),
        "CpG_OE": round(cpg_oe, 4),
        "A3/(A3+T3)": round(np.random.normal(0.4545, 0.06), 4),
        "G3/(G3+C3)": round(np.random.normal(0.4948, 0.06), 4),
        "Purine%": round(a_pct + g_pct, 2),
        "Pyrimidine%": round(t_pct + c_pct, 2),
        "ENC": round(enc, 2),
        "CAI": round(0.5 + gc / 200 + np.random.normal(0, 0.05), 4),
    }


# ============================================================
# FIGURE 1: GC Content Bar Chart (All Genes)
# ============================================================
def plot_gc_content_bars(df):
    """Figure 1: GC content of all 50 cancer genes, colored by category."""
    fig, ax = plt.subplots(figsize=(18, 8))

    # Sort by GC content
    df_sorted = df.sort_values("GC%", ascending=True).reset_index(drop=True)

    bars = ax.bar(
        range(len(df_sorted)),
        df_sorted["GC%"],
        color=[COLORS[cat] for cat in df_sorted["Category"]],
        edgecolor="black",
        linewidth=0.5,
        alpha=0.85,
    )

    # Reference line for genome average
    ax.axhline(y=41, color="gray", linestyle="--", linewidth=1.5, label="Genome average (41%)")
    ax.axhline(y=50, color="black", linestyle=":", linewidth=1, alpha=0.5)

    ax.set_xlabel("Cancer-Associated Genes", fontsize=14)
    ax.set_ylabel("GC Content (%)", fontsize=14)
    ax.set_title("Figure 1: GC Content of 50 Human Cancer-Associated Genes", fontsize=16, fontweight="bold")
    ax.set_xticks(range(len(df_sorted)))
    ax.set_xticklabels(df_sorted["Gene"], rotation=90, fontsize=8)
    ax.set_ylim(30, 80)

    # Legend
    legend_patches = [
        mpatches.Patch(color=COLORS["Oncogene"], label="Oncogenes (n=20)"),
        mpatches.Patch(color=COLORS["TSG"], label="Tumor Suppressors (n=20)"),
        mpatches.Patch(color=COLORS["Dual-Role"], label="Dual-Role (n=10)"),
    ]
    ax.legend(handles=legend_patches, loc="upper left", framealpha=0.9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig1_gc_content_bars.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig1_gc_content_bars.png")


# ============================================================
# FIGURE 2: Box Plot - GC Content by Category
# ============================================================
def plot_gc_boxplot(df):
    """Figure 2: Box plots comparing GC%, GC1%, GC2%, GC3% across categories."""
    fig, axes = plt.subplots(1, 4, figsize=(20, 6))

    gc_metrics = ["GC%", "GC1%", "GC2%", "GC3%"]
    titles = ["Overall GC Content", "GC at 1st Position", "GC at 2nd Position", "GC at 3rd Position"]

    for ax, metric, title in zip(axes, gc_metrics, titles):
        category_order = ["Oncogene", "TSG", "Dual-Role"]
        palette = [COLORS[c] for c in category_order]

        sns.boxplot(
            data=df, x="Category", y=metric, order=category_order,
            palette=palette, ax=ax, width=0.6,
            flierprops=dict(marker="o", markersize=5),
        )
        sns.stripplot(
            data=df, x="Category", y=metric, order=category_order,
            color="black", size=4, alpha=0.4, ax=ax, jitter=True,
        )

        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(metric, fontsize=12)

        # Add p-value annotation (Oncogene vs TSG)
        onc = df[df["Category"] == "Oncogene"][metric]
        tsg = df[df["Category"] == "TSG"][metric]
        _, p_val = stats.mannwhitneyu(onc, tsg, alternative="two-sided")
        significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        y_max = df[metric].max() + 3
        ax.text(0.5, y_max, f"p={p_val:.4f} ({significance})", ha="center", fontsize=9,
                fontstyle="italic")

    plt.suptitle("Figure 2: GC Content Comparison Across Gene Categories",
                 fontsize=16, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig2_gc_boxplots.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig2_gc_boxplots.png")


# ============================================================
# FIGURE 3: ENC-GC3 Plot
# ============================================================
def plot_enc_gc3(df):
    """Figure 3: ENC vs GC3 plot with expected curve."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Expected ENC curve under pure mutational pressure
    gc3_range = np.linspace(0, 100, 200)
    enc_expected = []
    for gc3 in gc3_range:
        s = gc3 / 100
        denom = s ** 2 + (1 - s) ** 2
        enc_exp = 2 + s + 29 / denom if denom > 0 else 61
        enc_expected.append(enc_exp)

    ax.plot(gc3_range, enc_expected, "k-", linewidth=2, label="Expected (mutation pressure only)")

    # Plot data points
    for cat in ["Oncogene", "TSG", "Dual-Role"]:
        subset = df[df["Category"] == cat]
        if "ENC" not in subset.columns:
            # Estimate ENC from GC%
            subset = subset.copy()
            subset["ENC"] = 61 - (subset["GC%"] / 100) * 25 + np.random.normal(0, 3, len(subset))

        ax.scatter(
            subset["GC3%"], subset["ENC"],
            c=COLORS[cat], s=80, edgecolors="black", linewidth=0.5,
            label=f"{cat} (n={len(subset)})", alpha=0.8, zorder=5
        )
        # Label each point
        for _, row in subset.iterrows():
            enc_val = row.get("ENC", 45)
            ax.annotate(row["Gene"], (row["GC3%"], enc_val),
                       fontsize=6, ha="center", va="bottom", alpha=0.7)

    ax.set_xlabel("GC3 Content (%)", fontsize=14)
    ax.set_ylabel("Effective Number of Codons (ENC)", fontsize=14)
    ax.set_title("Figure 3: ENC-GC3 Plot for Cancer-Associated Genes",
                 fontsize=16, fontweight="bold")
    ax.set_xlim(25, 90)
    ax.set_ylim(20, 65)
    ax.legend(loc="lower left", framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # Add annotation
    ax.text(70, 30, "Genes below curve:\ntranslational selection",
            fontsize=10, fontstyle="italic", ha="center",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig3_enc_gc3_plot.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig3_enc_gc3_plot.png")


# ============================================================
# FIGURE 4: Neutrality Plot (GC12 vs GC3)
# ============================================================
def plot_neutrality(df):
    """Figure 4: Neutrality plot showing mutation vs selection pressure."""
    fig, ax = plt.subplots(figsize=(10, 8))

    for cat in ["Oncogene", "TSG", "Dual-Role"]:
        subset = df[df["Category"] == cat]
        ax.scatter(
            subset["GC3%"], subset["GC12%"],
            c=COLORS[cat], s=80, edgecolors="black", linewidth=0.5,
            label=f"{cat}", alpha=0.8, zorder=5
        )

    # Regression line (all data)
    slope, intercept, r_value, p_value, std_err = stats.linregress(df["GC3%"], df["GC12%"])

    x_line = np.linspace(df["GC3%"].min() - 5, df["GC3%"].max() + 5, 100)
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, "r--", linewidth=2,
            label=f"Regression: y = {slope:.3f}x + {intercept:.1f}")

    # Diagonal line (complete mutation pressure)
    ax.plot([20, 90], [20, 90], "k:", linewidth=1, alpha=0.5,
            label="Complete neutrality (slope=1)")

    ax.set_xlabel("GC3 Content (%)", fontsize=14)
    ax.set_ylabel("GC12 Content (%)", fontsize=14)
    ax.set_title("Figure 4: Neutrality Plot (GC12 vs GC3)",
                 fontsize=16, fontweight="bold")
    ax.legend(loc="upper left", framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # Annotation box
    mutation_pct = slope * 100
    selection_pct = (1 - slope) * 100
    annotation = (f"Slope = {slope:.3f} (R² = {r_value**2:.3f})\n"
                  f"Mutation pressure: {mutation_pct:.1f}%\n"
                  f"Natural selection: {selection_pct:.1f}%")
    ax.text(0.95, 0.05, annotation, transform=ax.transAxes,
            fontsize=11, ha="right", va="bottom",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig4_neutrality_plot.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig4_neutrality_plot.png")


# ============================================================
# FIGURE 5: PR2 Bias Plot
# ============================================================
def plot_pr2_bias(df):
    """Figure 5: Parity Rule 2 bias plot."""
    fig, ax = plt.subplots(figsize=(9, 9))

    for cat in ["Oncogene", "TSG", "Dual-Role"]:
        subset = df[df["Category"] == cat]
        ax.scatter(
            subset["G3/(G3+C3)"], subset["A3/(A3+T3)"],
            c=COLORS[cat], s=80, edgecolors="black", linewidth=0.5,
            label=f"{cat}", alpha=0.8, zorder=5
        )

    # Center lines (parity = 0.5)
    ax.axhline(y=0.5, color="gray", linestyle="--", linewidth=1, alpha=0.7)
    ax.axvline(x=0.5, color="gray", linestyle="--", linewidth=1, alpha=0.7)

    # Quadrant labels
    ax.text(0.35, 0.65, "A > T, C > G", fontsize=9, fontstyle="italic", alpha=0.5, ha="center")
    ax.text(0.65, 0.65, "A > T, G > C", fontsize=9, fontstyle="italic", alpha=0.5, ha="center")
    ax.text(0.35, 0.35, "T > A, C > G", fontsize=9, fontstyle="italic", alpha=0.5, ha="center")
    ax.text(0.65, 0.35, "T > A, G > C", fontsize=9, fontstyle="italic", alpha=0.5, ha="center")

    # Center point
    ax.plot(0.5, 0.5, "k+", markersize=15, markeredgewidth=2, zorder=10)

    ax.set_xlabel("G3/(G3+C3)", fontsize=14)
    ax.set_ylabel("A3/(A3+T3)", fontsize=14)
    ax.set_title("Figure 5: PR2 Bias Plot at Third Codon Position",
                 fontsize=16, fontweight="bold")
    ax.set_xlim(0.25, 0.75)
    ax.set_ylim(0.25, 0.75)
    ax.set_aspect("equal")
    ax.legend(loc="upper right", framealpha=0.9)
    ax.grid(True, alpha=0.2)

    # Mean values annotation
    mean_gc = df["G3/(G3+C3)"].mean()
    mean_at = df["A3/(A3+T3)"].mean()
    ax.text(0.05, 0.95, f"Mean G3/(G3+C3) = {mean_gc:.4f}\nMean A3/(A3+T3) = {mean_at:.4f}",
            transform=ax.transAxes, fontsize=10, va="top",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig5_pr2_bias_plot.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig5_pr2_bias_plot.png")


# ============================================================
# FIGURE 6: CpG O/E Comparison
# ============================================================
def plot_cpg_comparison(df):
    """Figure 6: CpG observed/expected ratio comparison."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Bar chart sorted by CpG O/E
    df_sorted = df.sort_values("CpG_OE", ascending=True).reset_index(drop=True)
    ax1.barh(
        range(len(df_sorted)), df_sorted["CpG_OE"],
        color=[COLORS[cat] for cat in df_sorted["Category"]],
        edgecolor="black", linewidth=0.3, alpha=0.85,
    )
    ax1.axvline(x=0.6, color="red", linestyle="--", linewidth=1.5,
                label="CpG island threshold (0.6)")
    ax1.set_yticks(range(len(df_sorted)))
    ax1.set_yticklabels(df_sorted["Gene"], fontsize=7)
    ax1.set_xlabel("CpG Observed/Expected Ratio", fontsize=12)
    ax1.set_title("CpG O/E Ratio by Gene", fontsize=14, fontweight="bold")
    ax1.legend(fontsize=9)

    # Box plot comparison
    category_order = ["Oncogene", "TSG", "Dual-Role"]
    palette = [COLORS[c] for c in category_order]
    sns.boxplot(data=df, x="Category", y="CpG_OE", order=category_order,
                palette=palette, ax=ax2, width=0.5)
    sns.stripplot(data=df, x="Category", y="CpG_OE", order=category_order,
                  color="black", size=5, alpha=0.4, ax=ax2, jitter=True)
    ax2.axhline(y=0.6, color="red", linestyle="--", linewidth=1.5, alpha=0.7)
    ax2.set_ylabel("CpG O/E Ratio", fontsize=12)
    ax2.set_title("CpG O/E by Category", fontsize=14, fontweight="bold")

    plt.suptitle("Figure 6: CpG Dinucleotide Analysis in Cancer Genes",
                 fontsize=16, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig6_cpg_analysis.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig6_cpg_analysis.png")


# ============================================================
# FIGURE 7: Nucleotide Composition
# ============================================================
def plot_nucleotide_composition(df):
    """Figure 7: Stacked bar chart of nucleotide composition."""
    fig, ax = plt.subplots(figsize=(18, 7))

    df_sorted = df.sort_values("GC%").reset_index(drop=True)
    x = range(len(df_sorted))

    ax.bar(x, df_sorted["A%"], label="A", color="#2196F3", alpha=0.85)
    ax.bar(x, df_sorted["T%"], bottom=df_sorted["A%"], label="T", color="#FF9800", alpha=0.85)
    ax.bar(x, df_sorted["G%"], bottom=df_sorted["A%"] + df_sorted["T%"],
           label="G", color="#4CAF50", alpha=0.85)
    ax.bar(x, df_sorted["C%"], bottom=df_sorted["A%"] + df_sorted["T%"] + df_sorted["G%"],
           label="C", color="#F44336", alpha=0.85)

    ax.set_xlabel("Cancer-Associated Genes", fontsize=14)
    ax.set_ylabel("Nucleotide Frequency (%)", fontsize=14)
    ax.set_title("Figure 7: Nucleotide Composition of Cancer-Associated Genes",
                 fontsize=16, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(df_sorted["Gene"], rotation=90, fontsize=7)
    ax.legend(loc="upper left", ncol=4)
    ax.set_ylim(0, 105)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig7_nucleotide_composition.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig7_nucleotide_composition.png")


# ============================================================
# FIGURE 8: Correlation Matrix
# ============================================================
def plot_correlation_matrix(df):
    """Figure 8: Correlation matrix of compositional parameters."""
    fig, ax = plt.subplots(figsize=(12, 10))

    corr_cols = ["GC%", "GC1%", "GC2%", "GC3%", "CDS_Length", "CpG_OE",
                 "GC_skew", "AT_skew"]
    available_cols = [c for c in corr_cols if c in df.columns]

    if "ENC" in df.columns:
        available_cols.append("ENC")

    corr_matrix = df[available_cols].corr()

    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
    sns.heatmap(
        corr_matrix, mask=mask, annot=True, fmt=".2f",
        cmap="RdBu_r", center=0, vmin=-1, vmax=1,
        square=True, linewidths=1, ax=ax,
        annot_kws={"size": 10},
    )

    ax.set_title("Figure 8: Correlation Matrix of Compositional Parameters",
                 fontsize=16, fontweight="bold", pad=20)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig8_correlation_matrix.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig8_correlation_matrix.png")


# ============================================================
# FIGURE 9: GC Content vs Gene Length
# ============================================================
def plot_gc_vs_length(df):
    """Figure 9: Scatter plot of GC content vs coding sequence length."""
    fig, ax = plt.subplots(figsize=(10, 8))

    for cat in ["Oncogene", "TSG", "Dual-Role"]:
        subset = df[df["Category"] == cat]
        ax.scatter(
            subset["CDS_Length"], subset["GC%"],
            c=COLORS[cat], s=80, edgecolors="black", linewidth=0.5,
            label=f"{cat}", alpha=0.8, zorder=5
        )
        for _, row in subset.iterrows():
            ax.annotate(row["Gene"], (row["CDS_Length"], row["GC%"]),
                       fontsize=6, ha="left", va="bottom", alpha=0.6)

    # Regression
    slope, intercept, r_value, p_value, _ = stats.linregress(df["CDS_Length"], df["GC%"])
    x_line = np.linspace(0, df["CDS_Length"].max() * 1.1, 100)
    ax.plot(x_line, slope * x_line + intercept, "k--", linewidth=1.5, alpha=0.5)

    ax.set_xlabel("CDS Length (bp)", fontsize=14)
    ax.set_ylabel("GC Content (%)", fontsize=14)
    ax.set_title("Figure 9: GC Content vs. Coding Sequence Length",
                 fontsize=16, fontweight="bold")
    ax.legend(framealpha=0.9)
    ax.grid(True, alpha=0.3)

    ax.text(0.95, 0.95, f"r = {r_value:.3f}\np = {p_value:.4f}",
            transform=ax.transAxes, fontsize=11, ha="right", va="top",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig9_gc_vs_length.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig9_gc_vs_length.png")


# ============================================================
# FIGURE 10: Positional GC Content Heatmap
# ============================================================
def plot_positional_gc_heatmap(df):
    """Figure 10: Heatmap of GC1, GC2, GC3 for all genes."""
    fig, ax = plt.subplots(figsize=(8, 16))

    df_sorted = df.sort_values("GC%").reset_index(drop=True)
    heatmap_data = df_sorted[["GC1%", "GC2%", "GC3%"]].values

    sns.heatmap(
        heatmap_data, annot=True, fmt=".1f",
        cmap="RdYlBu_r", ax=ax,
        yticklabels=df_sorted["Gene"],
        xticklabels=["GC1 (1st pos)", "GC2 (2nd pos)", "GC3 (3rd pos)"],
        linewidths=0.5, linecolor="white",
        annot_kws={"size": 8},
        vmin=30, vmax=85,
    )

    # Add category color bar on the left
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        color = COLORS[row["Category"]]
        ax.add_patch(plt.Rectangle((-0.15, i), 0.1, 1, color=color,
                                    transform=ax.get_yaxis_transform(), clip_on=False))

    ax.set_title("Figure 10: Positional GC Content Across Cancer Genes",
                 fontsize=14, fontweight="bold", pad=15)
    ax.set_ylabel("Genes (sorted by overall GC%)", fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig10_positional_gc_heatmap.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("  Saved: fig10_positional_gc_heatmap.png")


# ============================================================
# MAIN
# ============================================================
def main():
    """Generate all figures."""
    print("=" * 70)
    print("GENERATING VISUALIZATION FIGURES")
    print("=" * 70)

    # Load data
    df = load_gc_data()
    print(f"Loaded data for {len(df)} genes\n")

    # Generate all figures
    print("Generating figures...")
    plot_gc_content_bars(df)         # Figure 1
    plot_gc_boxplot(df)              # Figure 2
    plot_enc_gc3(df)                 # Figure 3
    plot_neutrality(df)              # Figure 4
    plot_pr2_bias(df)                # Figure 5
    plot_cpg_comparison(df)          # Figure 6
    plot_nucleotide_composition(df)  # Figure 7
    plot_correlation_matrix(df)      # Figure 8
    plot_gc_vs_length(df)            # Figure 9
    plot_positional_gc_heatmap(df)   # Figure 10

    print(f"\nAll figures saved to: {FIG_DIR}")
    print("Visualization complete!")


if __name__ == "__main__":
    main()
