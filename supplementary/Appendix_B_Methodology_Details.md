# Appendix B: Detailed Methodology and Formulas

## B.1 Nucleotide Composition Formulas

### B.1.1 Individual Base Frequency
```
f(N) = Count(N) / L × 100

where N ∈ {A, T, G, C} and L = total sequence length
```

### B.1.2 Overall GC Content
```
GC% = (nG + nC) / (nA + nT + nG + nC) × 100
```

### B.1.3 Positional GC Content
For a coding sequence divided into codons (triplets), the GC content at each position:
```
GC1% = Σ(1 if codon[0] ∈ {G,C} else 0) / total_codons × 100
GC2% = Σ(1 if codon[1] ∈ {G,C} else 0) / total_codons × 100
GC3% = Σ(1 if codon[2] ∈ {G,C} else 0) / total_codons × 100

GC12% = (GC1% + GC2%) / 2
```

### B.1.4 GC3s (Third Position Synonymous GC)
```
GC3s% = (G3s + C3s) / total_synonymous_third_positions × 100

Excludes: ATG (Met), TGG (Trp), TAA/TAG/TGA (Stop)
```

## B.2 Nucleotide Skew Formulas

### B.2.1 GC Skew
```
GC-skew = (nG - nC) / (nG + nC)

Range: [-1, +1]
Positive: G-enriched on coding strand
Negative: C-enriched on coding strand
```

### B.2.2 AT Skew
```
AT-skew = (nA - nT) / (nA + nT)

Range: [-1, +1]
Positive: A-enriched on coding strand
Negative: T-enriched on coding strand
```

## B.3 Codon Usage Metrics

### B.3.1 Relative Synonymous Codon Usage (RSCU)
Sharp & Li (1987):
```
RSCU(ij) = X(ij) / [(1/ni) × Σ(j=1 to ni) X(ij)]

where:
  X(ij) = observed count of the jth codon for the ith amino acid
  ni = number of synonymous codons for amino acid i

Interpretation:
  RSCU = 1.0: no bias (codon used at expected frequency)
  RSCU > 1.0: overrepresented (preferred)
  RSCU < 1.0: underrepresented (avoided)
  RSCU > 1.6: significantly preferred
  RSCU < 0.6: significantly avoided
```

### B.3.2 Effective Number of Codons (ENC)
Wright (1990):
```
ENC = 2 + 9/F̄₂ + 1/F̄₃ + 5/F̄₄ + 3/F̄₆

where F̄ₖ = average homozygosity for k-fold degenerate amino acids

For each amino acid with k synonymous codons:
  F = Σ(pj²) for j=1 to k
  where pj = frequency of the jth synonymous codon

Corrected (sample-based):
  F̂ = (n × Σpj² - 1) / (n - 1)
  where n = total number of codons for that amino acid

Range: 20 (extreme bias) to 61 (no bias)
Strong bias: ENC < 35
Moderate bias: 35 ≤ ENC < 50
Low bias: ENC ≥ 50
```

### B.3.3 Expected ENC under Mutation Pressure
Novembre (2002):
```
ENC_expected = 2 + s + 29 / (s² + (1-s)²)

where s = GC3s (as fraction, 0 to 1)
```

### B.3.4 Codon Adaptation Index (CAI)
Sharp & Li (1987):
```
CAI = exp[(1/L) × Σ ln(w(codon_i))]

where:
  L = number of codons (excluding Met, Trp, Stop)
  w(codon) = RSCU(codon) / max(RSCU for that amino acid) in reference set

Range: 0 to 1
Higher CAI = better adapted to host translational machinery
```

## B.4 Evolutionary Force Analysis

### B.4.1 ENC-GC3 Plot Interpretation
```
Expected curve: ENC_exp = 2 + s + 29 / (s² + (1-s)²)

If gene falls ON the curve → mutation pressure alone determines codon usage
If gene falls BELOW the curve → additional selection pressure exists
Distance below curve = magnitude of selection effect
```

### B.4.2 Neutrality Plot Analysis
```
Linear regression: GC12 = slope × GC3 + intercept

Slope = 1.0 → complete mutation pressure (neutrality)
Slope = 0.0 → complete selection
Slope = b → mutation contributes b×100%, selection contributes (1-b)×100%

R² indicates how well mutation pressure explains GC12 variation
```

### B.4.3 Parity Rule 2 (PR2) Bias Analysis
Chargaff's second parity rule predicts A=T and G=C within single-stranded DNA:
```
PR2 plot axes:
  X-axis: G3/(G3+C3) → GC bias at 3rd position
  Y-axis: A3/(A3+T3) → AT bias at 3rd position

Center point (0.5, 0.5) = no bias (perfect parity)
Deviations indicate:
  - Strand-specific mutational pressure
  - Transcription-coupled repair bias
  - Translational selection
```

## B.5 CpG Content Analysis

### B.5.1 CpG Observed/Expected Ratio
```
CpG O/E = [f(CpG) × N] / [f(C) × f(G)]

where:
  f(CpG) = frequency of CpG dinucleotides = nCpG / (N-1)
  f(C) = frequency of cytosine = nC / N
  f(G) = frequency of guanine = nG / N
  N = total sequence length

CpG island criteria (Gardiner-Garden & Frommer, 1987):
  - Length ≥ 200 bp
  - GC% ≥ 50%
  - CpG O/E ≥ 0.6
```

## B.6 Sliding Window Analysis
```
For each window position i:
  GC(i) = [G(i) + C(i)] / window_size × 100

Parameters used:
  Window size: 100 bp
  Step size: 25 bp
  Center position: i + window_size/2
```

## B.7 Statistical Tests

### B.7.1 Mann-Whitney U Test
Non-parametric test for comparing two independent groups:
```
H₀: The distributions of the two groups are identical
H₁: The distributions differ

Used for: Oncogenes vs TSGs comparison
Chosen because: Data may not be normally distributed
```

### B.7.2 Kruskal-Wallis H Test
Non-parametric one-way ANOVA:
```
H₀: All groups have the same distribution
H₁: At least one group differs

Used for: Comparing all three categories (Oncogene, TSG, Dual-Role)
```

### B.7.3 Shapiro-Wilk Test
Tests for normality:
```
H₀: Data follows a normal distribution
H₁: Data does not follow a normal distribution

p > 0.05 → fail to reject H₀ (assume normal)
p ≤ 0.05 → reject H₀ (non-normal)
```

### B.7.4 Fisher's Exact Test
Tests association in 2×2 contingency tables:
```
Used for: CpG island classification (Yes/No) × Gene category (Oncogene/TSG)
```

### B.7.5 Pearson Correlation
```
r = Σ[(xi - x̄)(yi - ȳ)] / √[Σ(xi - x̄)² × Σ(yi - ȳ)²]

Interpretation:
  |r| ≥ 0.7: Strong correlation
  0.4 ≤ |r| < 0.7: Moderate correlation
  0.2 ≤ |r| < 0.4: Weak correlation
  |r| < 0.2: Negligible correlation
```

### B.7.6 Bonferroni Correction
```
α_corrected = α / k

where k = number of simultaneous tests
Used to control family-wise error rate in multiple comparisons
```

## B.8 Software and Tools Reference

| Software | Version | Usage | Reference |
|----------|---------|-------|-----------|
| Python | 3.11 | Programming | python.org |
| BioPython | 1.83 | Sequence I/O, GC calculation | Cock et al. (2009) |
| NumPy | 1.26 | Array operations | Harris et al. (2020) |
| Pandas | 2.1 | Data frames | McKinney (2010) |
| Matplotlib | 3.8 | Plotting | Hunter (2007) |
| Seaborn | 0.13 | Statistical plots | Waskom (2021) |
| SciPy | 1.11 | Statistical tests | Virtanen et al. (2020) |
| CodonW | 1.4.2 | RSCU, ENC validation | Peden (1999) |
| MEGA | 11 | Codon statistics | Kumar et al. (2018) |
| EMBOSS | 6.6.0 | Sequence composition | Rice et al. (2000) |
