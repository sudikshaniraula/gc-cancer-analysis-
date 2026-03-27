# Computational Analysis of GC Content Variability and Nucleotide Composition in Human Cancer-Associated Genes

---

**Author:** Sudiksha Niraula
Self Directed undergraduate research
**Date:** March 2026
**Duration of Study:** April 2025 – March 2026

---

## Abstract

The nucleotide composition of protein-coding genes profoundly influences gene expression, mRNA stability, codon usage, and mutational susceptibility. In cancer biology, understanding the compositional architecture of cancer-associated genes provides critical insights into the molecular mechanisms underlying oncogenesis. This study presents a comprehensive computational analysis of GC content variability, nucleotide composition patterns, codon usage bias, and compositional skew in **50 well-characterized human cancer-associated genes** spanning oncogenes, tumor suppressor genes (TSGs), and dual-role genes. Gene sequences were retrieved from the NCBI Gene Database and cross-referenced with the COSMIC Cancer Gene Census (Tier 1). Using BioPython, CodonW, and custom Python scripts, we computed GC content at overall and codon-positional levels (GC1, GC2, GC3), relative synonymous codon usage (RSCU), effective number of codons (ENC), nucleotide skew (GC-skew and AT-skew), and CpG observed/expected ratios. Sliding window analysis revealed positional GC content heterogeneity across gene bodies. Our results demonstrate that **oncogenes exhibit higher mean GC content (51.84 ± 8.80%) compared to tumor suppressor genes (48.81 ± 10.03%)**, although this difference did not reach statistical significance (p = 0.323, Mann-Whitney U test). The trend was more pronounced at the third codon position (GC3: 59.97% vs. 50.61%, p = 0.086). ENC-GC3 plot analysis revealed that 84% of cancer genes fall below the expected curve under pure mutational pressure, indicating that translational selection operates in addition to compositional constraints. Neutrality plot analysis showed a regression coefficient of 0.200, indicating that natural selection is the dominant force shaping codon usage, with mutation pressure contributing approximately 20.0% and natural selection 80.0%. Parity rule 2 (PR2) bias plot analysis revealed that 92% of genes deviate from equal nucleotide usage at the third codon position, suggesting asymmetric mutational pressures and selection constraints acting on cancer-associated genes. This work provides a foundational computational framework for understanding the compositional biology of cancer genes and its implications for gene regulation, epigenetic modification, and mutational vulnerability.

**Keywords:** GC content, nucleotide composition, cancer genes, codon usage bias, RSCU, ENC, CpG islands, oncogenes, tumor suppressor genes, computational biology

---

## 1. Introduction

### 1.1 Background and Rationale

Cancer remains one of the leading causes of mortality worldwide, with an estimated 19.3 million new cases and nearly 10 million deaths annually (Sung et al., 2021). At the molecular level, cancer is driven by somatic and germline mutations in genes that regulate cell growth, differentiation, apoptosis, and genomic stability. These cancer-associated genes broadly fall into two functional categories: **oncogenes**, which promote cell proliferation when activated, and **tumor suppressor genes (TSGs)**, which normally restrain cell growth and whose loss of function contributes to malignancy (Vogelstein et al., 2013).

The nucleotide composition of a gene — specifically its GC content (the proportion of guanine and cytosine bases) — is far from a trivial structural characteristic. GC content profoundly influences multiple biological processes:

- **Gene expression levels:** GC-rich genes tend to be more highly expressed, partly due to their association with open chromatin structures and active transcription (Kudla et al., 2006).
- **mRNA stability and secondary structure:** Higher GC content increases mRNA thermodynamic stability, which can affect translational efficiency (Mauger et al., 2019).
- **Codon usage bias:** The third codon position GC content (GC3) is a major determinant of synonymous codon preferences, which in turn affects translational accuracy and speed (Hershberg & Petrov, 2008).
- **CpG dinucleotide frequency:** CpG sites are critical substrates for DNA methylation, an epigenetic modification frequently disrupted in cancer through promoter hypermethylation of tumor suppressors or global hypomethylation (Baylin & Jones, 2016).
- **Mutational vulnerability:** Methylated CpG dinucleotides are mutation hotspots due to spontaneous deamination of 5-methylcytosine to thymine, contributing to C→T transition mutations that are among the most common somatic mutations in cancer (Alexandrov et al., 2013).

### 1.2 GC Content in Genomic Context

The human genome has an average GC content of approximately 41% (Piovesan et al., 2019), but this value varies dramatically across chromosomal regions, forming compositionally distinct domains called **isochores** — megabase-scale regions of relatively homogeneous base composition (Bernardi, 2000). GC-rich isochores (H isochores) are gene-dense, early-replicating, and associated with active chromatin, while GC-poor isochores (L isochores) are gene-poor and late-replicating.

Within coding sequences, GC content shows characteristic positional variation across the three codon positions:
- **GC1 (first position):** Moderately constrained by amino acid identity
- **GC2 (second position):** Most constrained, as changes here almost always alter the encoded amino acid
- **GC3 (third position):** Least constrained due to codon degeneracy, making it the most sensitive indicator of mutational and selectional forces

### 1.3 Codon Usage Bias and Cancer

Codon usage bias — the non-random preferential use of synonymous codons — is a universal phenomenon observed across all domains of life (Sharp et al., 2010). In human genes, codon usage is influenced by:

1. **Mutational pressure** (compositional bias toward AT or GC)
2. **Translational selection** (preference for codons recognized by abundant tRNAs)
3. **GC-biased gene conversion** (a recombination-associated process favoring GC alleles)

Several metrics quantify codon usage bias:
- **Relative Synonymous Codon Usage (RSCU):** The ratio of observed codon frequency to expected frequency under uniform usage. RSCU > 1.0 indicates overrepresentation; RSCU < 1.0 indicates underrepresentation (Sharp & Li, 1987).
- **Effective Number of Codons (ENC):** Ranges from 20 (extreme bias, one codon per amino acid) to 61 (no bias, all synonymous codons used equally). ENC < 35 indicates strong bias (Wright, 1990).
- **Codon Adaptation Index (CAI):** Measures the extent to which codon usage matches that of highly expressed genes (Sharp & Li, 1987).

Recent studies have revealed that codon usage patterns differ between oncogenes and tumor suppressors, and between cancer cells and normal tissue, suggesting that translational reprogramming may contribute to oncogenic transformation (Supek et al., 2014).

### 1.4 CpG Islands and Cancer Epigenetics

CpG islands (CGIs) are regions of the genome with CpG dinucleotide frequency at or above the statistically expected level (~4-6%), in contrast to the bulk genome where CpG is ~5-fold depleted (~1%) (Gardiner-Garden & Frommer, 1987). Approximately 70% of human gene promoters are associated with CGIs. In cancer:

- **Promoter hypermethylation** of CGIs silences tumor suppressor genes (e.g., RB1, CDKN2A, MLH1, BRCA1)
- **Global hypomethylation** activates oncogenes and transposable elements
- **Aberrant CpG methylation patterns** serve as diagnostic and prognostic biomarkers

Understanding the CpG content and distribution within cancer gene coding sequences provides insights into their epigenetic regulation and mutational susceptibility.

### 1.5 Objectives of the Study

This study aims to:

1. Retrieve and curate coding sequences (CDS) of 50 well-characterized human cancer-associated genes from the NCBI database
2. Compute overall GC content and positional GC content (GC1, GC2, GC3) for each gene
3. Analyze nucleotide composition patterns (individual base frequencies, purine/pyrimidine ratios)
4. Calculate codon usage metrics (RSCU, ENC, CAI) and identify preferred/avoided codons
5. Perform ENC-GC3 plot analysis to assess the relative contributions of mutation pressure and translational selection
6. Conduct neutrality plot analysis (GC12 vs. GC3) to quantify evolutionary forces
7. Perform Parity Rule 2 (PR2) bias plot analysis to detect asymmetric mutational pressures
8. Analyze GC-skew and AT-skew patterns across gene bodies using sliding window approaches
9. Calculate CpG observed/expected ratios to assess CpG island characteristics within coding regions
10. Compare compositional features between oncogenes and tumor suppressor genes using statistical tests

### 1.6 Significance of the Study

This research contributes to the growing field of cancer genomics by providing a systematic computational characterization of the compositional landscape of cancer genes. Understanding GC content variability and nucleotide composition patterns has practical implications for:

- Identifying mutational hotspots susceptible to C→T transitions
- Understanding differential gene expression regulation in cancer
- Improving cancer gene panel design for diagnostic sequencing
- Informing codon optimization strategies for cancer gene therapy constructs
- Contributing to the broader understanding of genome evolution and compositional heterogeneity

---

## 2. Literature Review

### 2.1 Genome Composition and Isochore Theory

The concept of compositional heterogeneity in mammalian genomes was first proposed by Bernardi and colleagues, who identified large-scale chromosomal domains (isochores) with distinct GC content (Bernardi, 2000). The human genome contains five families of isochores: two GC-poor families (L1, L2) and three GC-rich families (H1, H2, H3). GC-rich isochores harbor a disproportionate share of genes, particularly housekeeping genes, and are associated with early replication timing, open chromatin, and high recombination rates (Costantini et al., 2006).

The isochore structure has direct implications for cancer genomics, as somatic mutation rates vary across isochores, with GC-rich regions showing distinct mutational signatures compared to GC-poor regions (Lawrence et al., 2013).

### 2.2 GC Content and Gene Expression

Multiple studies have established a positive correlation between GC content and gene expression levels. Kudla et al. (2006) demonstrated that synonymous codon changes that increase GC content systematically increase protein expression in human cells, independent of mRNA levels. This effect is partly attributable to enhanced mRNA export and translation initiation efficiency.

In cancer, expression levels of oncogenes and tumor suppressors are critical determinants of cellular phenotype. The observation that oncogenes tend to have higher GC content than tumor suppressors (D'Onofrio et al., 1999) raises the possibility that compositional differences contribute to their differential expression patterns.

### 2.3 Codon Usage in Human Disease Genes

Several investigators have examined codon usage patterns in disease-associated genes. Plotkin and Kudla (2011) provided a comprehensive review of synonymous codon usage and its biological consequences. More recently, Buhr et al. (2016) demonstrated that synonymous codons can direct differential protein folding, establishing a direct mechanistic link between codon usage and protein function.

In the context of cancer specifically, Supek et al. (2014) showed that codon usage changes in tumors are associated with tRNA expression changes, suggesting a co-adaptive mechanism. Goodarzi et al. (2016) identified a regulatory program where tRNA modifications promote metastasis by enhancing translation of specific codon-enriched transcripts.

### 2.4 Nucleotide Composition Analysis in Disease Genes

Prabha et al. (2017) analyzed codon usage patterns in genes associated with various human cancers and found significant compositional biases. Their work revealed that cancer-related genes show distinct GC3 profiles compared to non-cancer genes. Similarly, comprehensive studies on autophagy genes (Malik et al., 2022), apoptosis genes, and immune response genes have employed the analytical framework that this study extends to cancer-associated genes.

The study by Chakraborty et al. (2019) on oral cancer genes demonstrated that nucleotide composition analysis can reveal evolutionary forces acting on cancer genes, with their ENC values ranging from 34.6 to 55.9 (mean 49.03 ± 4.22), indicating generally low codon usage bias.

### 2.5 CpG Depletion and Mutational Patterns in Cancer

The depletion of CpG dinucleotides in the mammalian genome is a well-documented consequence of cytosine methylation and subsequent deamination. In cancer, this process generates C→T transition mutations at CpG sites, which represent approximately 28% of all somatic point mutations in cancer (Alexandrov et al., 2013). The CpG observed/expected ratio (CpG O/E) serves as a quantitative measure of CpG depletion, with values < 1.0 indicating depletion relative to random expectation.

### 2.6 Computational Tools for Nucleotide Composition Analysis

Modern bioinformatics provides robust tools for sequence composition analysis:
- **BioPython** (Cock et al., 2009): Python library for biological computation, including GC content calculation via `Bio.SeqUtils.gc_fraction()`
- **CodonW** (Peden, 1999): Standalone program for codon usage analysis including RSCU, ENC, and correspondence analysis
- **EMBOSS** (Rice et al., 2000): Suite of bioinformatics tools including `cusp` for codon usage and `compseq` for composition analysis
- **MEGA** (Kumar et al., 2018): Molecular Evolutionary Genetics Analysis software for codon usage statistics
- **CAIcal** (Puigbò et al., 2008): Web server for codon adaptation index calculation

---

## 3. Materials and Methods

### 3.1 Gene Selection and Classification

A panel of 50 human cancer-associated genes was selected based on the following criteria:
- Listed in the **COSMIC Cancer Gene Census Tier 1** (Sondka et al., 2018)
- Well-characterized functional role as oncogene, tumor suppressor, or dual-role gene
- Validated across multiple cancer types with strong experimental evidence
- Representation across major cancer signaling pathways

The selected genes were classified into three functional categories:

**Table 1: Selected Cancer-Associated Genes**

| Category | Genes (n) | Gene Symbols |
|----------|-----------|--------------|
| **Oncogenes** (n=20) | 20 | KRAS, BRAF, MYC, EGFR, ERBB2, PIK3CA, ABL1, RET, KIT, MET, ALK, NRAS, HRAS, FLT3, JAK2, FGFR1, FGFR2, FGFR3, CDK4, CCND1 |
| **Tumor Suppressor Genes** (n=20) | 20 | TP53, BRCA1, BRCA2, RB1, APC, PTEN, VHL, WT1, NF1, NF2, CDKN2A, SMAD4, STK11, MLH1, MSH2, ATM, BAP1, CHEK2, CDH1, FBXW7 |
| **Dual-Role / Other** (n=10) | 10 | NOTCH1, IDH1, IDH2, EZH2, DNMT3A, TET2, SF3B1, NPM1, CTNNB1, MDM2 |

### 3.2 Sequence Retrieval

Coding DNA sequences (CDS) for all 50 genes were retrieved from the **NCBI Gene Database** (https://www.ncbi.nlm.nih.gov/gene/) and **NCBI Nucleotide Database** (https://www.ncbi.nlm.nih.gov/nuccore/) using the following protocol:

1. Each gene was searched by its official HGNC symbol in the NCBI Gene database
2. The canonical RefSeq mRNA transcript (NM_ accession) was identified
3. The complete coding sequence (CDS) was extracted from the corresponding RefSeq entry
4. Sequences were downloaded in FASTA format
5. Sequence integrity was verified by confirming:
   - Start codon (ATG) presence
   - Stop codon (TAA, TAG, or TGA) presence
   - Coding sequence length divisible by 3 (complete codons)
   - No internal stop codons

All sequences correspond to the human reference genome assembly GRCh38 (hg38).

### 3.3 Nucleotide Composition Analysis

#### 3.3.1 Overall Base Frequencies
For each gene, the frequency of each nucleotide (A, T, G, C) was calculated as:

```
f(N) = Count(N) / Total_Length × 100
```

where N ∈ {A, T, G, C}

#### 3.3.2 GC Content Calculation
Overall GC content was computed as:

```
GC% = (G + C) / (A + T + G + C) × 100
```

#### 3.3.3 Positional GC Content
GC content at each codon position was calculated separately:

- **GC1:** GC content at the first codon position
- **GC2:** GC content at the second codon position
- **GC3:** GC content at the third codon position

```
GCn% = (Gn + Cn) / Total_codons × 100
```

where n ∈ {1, 2, 3} represents the codon position.

#### 3.3.4 GC Content at Third Synonymous Position (GC3s)
GC3s was calculated considering only synonymous third-position sites (excluding Met and Trp codons, which have no synonymous alternatives).

### 3.4 Nucleotide Skew Analysis

#### 3.4.1 GC Skew
GC skew measures the asymmetry between G and C on the coding strand:

```
GC-skew = (G - C) / (G + C)
```

Positive values indicate G enrichment; negative values indicate C enrichment.

#### 3.4.2 AT Skew
AT skew measures the asymmetry between A and T:

```
AT-skew = (A - T) / (A + T)
```

#### 3.4.3 Sliding Window Analysis
GC content was computed across each gene using a sliding window approach:
- **Window size:** 100 bp
- **Step size:** 25 bp
- GC content was calculated for each window and plotted against genomic position to visualize positional heterogeneity

### 3.5 Codon Usage Analysis

#### 3.5.1 Relative Synonymous Codon Usage (RSCU)
RSCU was calculated using the formula (Sharp & Li, 1987):

```
RSCU(ij) = X(ij) / [(1/n_i) × Σ X(ij)]
```

where:
- X(ij) = observed frequency of the jth codon for the ith amino acid
- n_i = number of synonymous codons for the ith amino acid

RSCU = 1.0 indicates no bias. RSCU > 1.6 indicates an overrepresented (preferred) codon. RSCU < 0.6 indicates an underrepresented (avoided) codon.

#### 3.5.2 Effective Number of Codons (ENC)
ENC was calculated using Wright's formula (Wright, 1990):

```
ENC = 2 + 9/F₂ + 1/F₃ + 5/F₄ + 3/F₆
```

where F_k represents the average homozygosity for amino acids with k-fold degeneracy. ENC ranges from 20 (maximum bias) to 61 (no bias).

#### 3.5.3 Codon Adaptation Index (CAI)
CAI was calculated using the reference codon usage table for highly expressed human genes, obtained from the Codon Usage Database (Kazusa, http://www.kazusa.or.jp/codon/).

### 3.6 Evolutionary Force Analysis

#### 3.6.1 ENC-GC3 Plot Analysis
The ENC values were plotted against GC3 values to assess whether codon usage bias is driven by compositional constraints (mutation pressure) alone. The expected ENC under pure mutational pressure follows the curve (Novembre, 2002):

```
ENC_expected = 2 + s + 29 / [s² + (1-s)²]
```

where s = GC3 content. Genes falling on or near the expected curve are predominantly influenced by mutation pressure, while genes falling significantly below the curve are subject to additional selection pressure.

#### 3.6.2 Neutrality Plot Analysis (GC12 vs. GC3)
The average GC content at the first and second codon positions (GC12) was plotted against GC3. Under complete mutational pressure, GC12 and GC3 should show a perfect correlation with a regression slope of 1.0. Deviations from slope = 1.0 indicate the influence of natural selection:

```
Mutation pressure contribution = Regression slope × 100%
Natural selection contribution = (1 - Regression slope) × 100%
```

#### 3.6.3 Parity Rule 2 (PR2) Bias Plot
PR2 analysis examines deviations from Chargaff's second parity rule (A=T, G=C within single-stranded DNA). A PR2 bias plot displays:
- **X-axis:** G/(G+C) — GC bias at the third codon position
- **Y-axis:** A/(A+T) — AT bias at the third codon position

The center point (0.5, 0.5) represents no bias. Deviations from this center indicate asymmetric mutational pressures and/or selection.

### 3.7 CpG Content Analysis

The CpG observed/expected ratio was calculated as:

```
CpG O/E = (f(CpG) × N) / (f(C) × f(G))
```

where:
- f(CpG) = frequency of CpG dinucleotides
- f(C) = frequency of cytosine
- f(G) = frequency of guanine
- N = total sequence length

CpG O/E ≥ 0.6 with GC content ≥ 50% is indicative of CpG island characteristics (Gardiner-Garden & Frommer, 1987).

### 3.8 Statistical Analysis

All statistical analyses were performed using Python (SciPy library) and SPSS 23.0:

- **Descriptive statistics:** Mean, standard deviation, range for all compositional parameters
- **Correlation analysis:** Karl Pearson's correlation coefficient for relationships between compositional parameters (GC content vs. ENC, GC3 vs. GC12, gene length vs. GC content, etc.)
- **Comparative statistics:** Mann-Whitney U test for comparing oncogene vs. TSG compositional features (non-parametric, as distributions may not be normal)
- **Kruskal-Wallis test:** For comparing across all three gene categories (oncogenes, TSGs, dual-role)
- **Multiple comparison correction:** Bonferroni correction for multiple testing
- **Significance level:** p < 0.05 was considered statistically significant

### 3.9 Computational Tools and Software

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.11 | Primary programming language |
| BioPython | 1.83 | Sequence parsing, GC content calculation |
| NumPy | 1.26 | Numerical computations |
| Pandas | 2.1 | Data manipulation and analysis |
| Matplotlib | 3.8 | Data visualization |
| Seaborn | 0.13 | Statistical visualization |
| SciPy | 1.11 | Statistical testing |
| CodonW | 1.4.2 | RSCU, ENC, CAI calculation |
| MEGA | 11 | Codon usage statistics validation |
| EMBOSS | 6.6.0 | Sequence composition analysis |
| CAIcal | 2.0 | CAI server for validation |
| SPSS | 23.0 | Statistical analysis |
| Origin Pro | 2024 | Plot generation |

### 3.10 Workflow Summary

```
Gene Selection (COSMIC CGC Tier 1)
        ↓
Sequence Retrieval (NCBI RefSeq CDS, FASTA format)
        ↓
Quality Control (CDS integrity verification)
        ↓
    ┌───────────────────┬──────────────────────┐
    ↓                   ↓                      ↓
Nucleotide          Codon Usage           CpG Content
Composition         Analysis              Analysis
(GC%, GC1/2/3,     (RSCU, ENC, CAI)      (CpG O/E ratio)
 skew, sliding
 window)
    ↓                   ↓                      ↓
    └───────────────────┼──────────────────────┘
                        ↓
            Evolutionary Force Analysis
            (ENC-GC3, Neutrality plot, PR2)
                        ↓
            Statistical Comparison
            (Oncogenes vs. TSGs)
                        ↓
            Visualization & Interpretation
```

---

## 4. Results

### 4.1 Gene Characteristics Overview

The 50 cancer-associated genes analyzed varied substantially in coding sequence length, ranging from 471 bp (CDKN2A) to 10,257 bp (BRCA2). The mean CDS length was 2,940 ± 2,289 bp. Oncogenes had a mean CDS length of 2,511 ± 1,301 bp, while tumor suppressor genes had a longer mean CDS length of 3,418 ± 3,127 bp (Table 2).

**Table 2: Summary of Gene Characteristics**

| Parameter | Oncogenes (n=20) | TSGs (n=20) | Dual-Role (n=10) | All Genes (n=50) |
|-----------|-------------------|-------------|-------------------|-------------------|
| Mean CDS length (bp) | 2,511 ± 1,301 | 3,418 ± 3,127 | 2,992 ± 2,242 | 2,940 ± 2,289 |
| Min CDS length (bp) | 567 (KRAS) | 471 (CDKN2A) | 885 (NPM1) | 471 |
| Max CDS length (bp) | 4,863 (ALK) | 10,257 (BRCA2) | 7,668 (NOTCH1) | 10,257 |

### 4.2 Overall GC Content

The overall GC content of the 50 cancer genes ranged from 35.73% (BRCA2) to 71.97% (CDKN2A), with a mean of 49.93 ± 9.24%. This is notably higher than the genome-wide average of ~41%, consistent with the known gene-dense, GC-rich isochore localization of many cancer genes.

**Table 3: GC Content by Gene Category**

| Parameter | Oncogenes | TSGs | Dual-Role | p-value* |
|-----------|-----------|------|-----------|----------|
| Mean GC% | 51.84 ± 8.80 | 48.81 ± 10.03 | 49.13 ± 8.91 | 0.3235 |
| Mean GC1% | 54.90 ± 5.88 | 54.89 ± 8.48 | 55.21 ± 4.50 | 0.6750 |
| Mean GC2% | 40.66 ± 4.90 | 40.92 ± 6.65 | 41.20 ± 4.83 | 0.8817 |
| Mean GC3% | 59.97 ± 18.28 | 50.61 ± 18.43 | 50.97 ± 20.68 | 0.0859 |
| Mean GC3s% | 58.53 ± 19.12 | 49.13 ± 19.02 | 49.14 ± 21.38 | 0.0962 |

*Mann-Whitney U test (Oncogenes vs. TSGs)

**Key Finding:** Oncogenes showed a trend toward higher GC content than tumor suppressor genes, particularly at the third codon position (GC3: 59.97% vs. 50.61%, p = 0.086), approaching but not reaching conventional statistical significance (α = 0.05). The first and second codon positions showed nearly identical means across categories, indicating that the compositional divergence is primarily driven by synonymous (wobble) position variation.

#### 4.2.1 GC Content of Individual Genes

**Highest GC content genes (>60%):**
- CDKN2A (72.0%), FGFR3 (65.8%), NOTCH1 (65.8%), WT1 (61.7%), VHL (61.5%), STK11 (61.3%), CCND1 (61.1%), ERBB2 (60.8%), RET (60.2%), DNMT3A (60.2%)
- These include both oncogenes (FGFR3, CCND1, ERBB2, RET) and TSGs (CDKN2A, WT1, VHL, STK11), indicating that high GC content is not exclusive to oncogenes

**Lowest GC content genes (<40%):**
- BRCA2 (35.7%), KRAS (37.6%), PTEN (38.4%), ATM (38.8%), JAK2 (38.7%), RB1 (39.0%), PIK3CA (39.1%)
- Predominantly large tumor suppressor genes (BRCA2, ATM, RB1) involved in DNA repair, but also includes oncogenes (KRAS, PIK3CA, JAK2)

### 4.3 Nucleotide Composition

#### 4.3.1 Individual Base Frequencies

The mean nucleotide frequencies across all 50 genes were:
- **A:** 27.11%
- **T:** 22.57%
- **G:** 25.87%
- **C:** 24.45%

The purine/pyrimidine ratio (AG/TC) was approximately 1.13, showing a purine excess on the coding strand consistent with transcription-associated mutational asymmetry.

#### 4.3.2 Compositional Differences Between Gene Categories

| Nucleotide | Oncogenes | TSGs | p-value |
|------------|-----------|------|---------|
| A% | 26.22 | 28.53 | 0.1478 |
| T% | 21.94 | 22.66 | 0.6168 |
| G% | 26.59 | 25.27 | 0.2085 |
| C% | 25.25 | 23.54 | 0.3369 |

Oncogenes showed a trend toward higher G and C frequencies and lower A frequencies compared to TSGs, although these differences were not statistically significant at α = 0.05.

### 4.4 Nucleotide Skew Analysis

#### 4.4.1 GC Skew
The mean GC skew across all genes was positive, indicating slight G-enrichment on the coding strand, consistent with known transcription-associated mutational asymmetry. The GC-skew and AT-skew showed a weak positive correlation (r = 0.261, p = 0.067), suggesting partially coordinated strand-specific biases.

#### 4.4.2 AT Skew
The mean AT skew was positive (A > T on coding strand), consistent with the general pattern of purine loading on the coding strand of human genes.

#### 4.4.3 Sliding Window Analysis
Sliding window analysis (100 bp windows, 25 bp steps) was performed on representative genes from each category, revealing substantial intra-genic GC content heterogeneity. Several notable patterns emerged:

- **FGFR3** (GC = 65.8%): Consistently high GC throughout the coding region with peaks exceeding 75% in the kinase domain
- **BRCA2** (GC = 35.7%): Low GC throughout with fluctuations between 28-45%, reflecting its location in a GC-poor isochore
- **NOTCH1** (GC = 65.8%): High GC with characteristic peaks in the EGF-like repeat domain-encoding region
- **KRAS** (GC = 37.6%): Surprisingly low GC for an oncogene, with modest peaks near the codon 12/13 mutation hotspot region

### 4.5 Codon Usage Analysis

#### 4.5.1 RSCU Analysis

RSCU values were computed for all 59 sense codons (excluding ATG for Met and TGG for Trp). Across all 50 cancer genes:

**Overrepresented codons (RSCU > 1.6):**
| Codon | Amino Acid | Mean RSCU | Category |
|-------|------------|-----------|----------|
| CTG | Leucine | 2.182 | GC-ending |
| GTG | Valine | 1.763 | GC-ending |

**Underrepresented codons (RSCU < 0.6):**
| Codon | Amino Acid | Mean RSCU | Category |
|-------|------------|-----------|----------|
| GCG | Alanine | 0.355 | CG-ending |
| TCG | Serine | 0.386 | CG-ending |
| CCG | Proline | 0.430 | CG-ending |
| CTA | Leucine | 0.457 | AT-ending |
| ACG | Threonine | 0.471 | CG-ending |
| CGT | Arginine | 0.531 | AT-ending |
| ATA | Isoleucine | 0.551 | AT-ending |
| CAA | Glutamine | 0.554 | AT-ending |

**Key observation:** Only two codons were strongly overrepresented (CTG and GTG, both GC-ending), reflecting the preference for specific GC-ending codons rather than a broad trend. The most avoided codons include both CG-ending codons (GCG, TCG, CCG, ACG — reflecting CpG dinucleotide avoidance due to methylation-mediated depletion) and AT-ending codons (CTA, ATA, CAA — reflecting compositional bias).

#### 4.5.2 RSCU Comparison: Oncogenes vs. TSGs

Oncogenes showed a trend toward higher RSCU values for GC-ending codons compared to TSGs, consistent with their higher mean GC3 content, although the difference was not individually significant for most codons due to high within-group variance.

#### 4.5.3 ENC Analysis

The ENC values ranged from 33.4 (NOTCH1) to 58.1 (IDH1), with a mean of 48.0 ± 6.2.

| Category | Mean ENC | Range |
|----------|----------|-------|
| Oncogenes | 46.73 ± 6.42 | 33.64 - 55.29 |
| TSGs | 49.54 ± 4.74 | 36.84 - 56.11 |
| Dual-Role | 48.55 ± 7.48 | 33.44 - 58.15 |

No gene showed extreme codon usage bias (ENC < 30), but several genes exhibited moderate bias (ENC 33-40), including NOTCH1 (33.4), CCND1 (33.6), FGFR3 (34.8), HRAS (35.4), and STK11 (36.8). Oncogenes showed a trend toward lower ENC values (stronger codon usage bias) than TSGs, though this difference did not reach statistical significance (p = 0.126).

### 4.6 Evolutionary Force Analysis

#### 4.6.1 ENC-GC3 Plot

The ENC-GC3s plot revealed that the majority of cancer genes (42/50, 84%) fell below the expected curve under pure mutational pressure. This indicates that while compositional bias (mutation pressure) is a significant factor, additional selective forces also influence codon usage in cancer genes.

- **On or near the curve (within 10% of expected):** 8 genes (16%) — predominantly genes with intermediate GC3s values
- **Below the curve:** 42 genes (84%) — distributed across all three functional categories

**Interpretation:** The high proportion of genes falling below the Wright (1990) expected curve indicates that translational selection operates in addition to mutational pressure across all cancer gene categories, consistent with the functional importance and high expression requirements of cancer-associated genes.

#### 4.6.2 Neutrality Plot Analysis

The neutrality plot (GC12 vs. GC3) yielded a regression equation:

```
GC12 = 0.200 × GC3 + 37.02    (R² = 0.474, p = 3.29 × 10⁻⁸)
```

The regression slope of 0.200 indicates:
- **Mutation pressure contribution:** 20.0%
- **Natural selection contribution:** 80.0%

This demonstrates that natural selection is the dominant force shaping codon usage in human cancer genes, contributing approximately 80% of the codon usage variation. This is consistent with the functional importance and high expression levels of many cancer-associated genes.

When analyzed separately:
- Oncogenes: slope = 0.201 (R² = 0.576), mutation pressure = 20.1%, selection = 79.9%
- TSGs: slope = 0.277 (R² = 0.548), mutation pressure = 27.7%, selection = 72.3%
- Dual-Role: slope = 0.130 (R² = 0.447), mutation pressure = 13.0%, selection = 87.0%

#### 4.6.3 PR2 Bias Plot

The PR2 bias plot (A3/(A3+T3) vs. G3/(G3+C3)) showed:

- **Mean A3/(A3+T3):** 0.4603 ± 0.064 (deviation from 0.5 indicates slight T3 > A3)
- **Mean G3/(G3+C3):** 0.5113 ± 0.055 (slight G3 > C3)

The majority of genes (46/50, 92%) deviated from the center point (0.5, 0.5), confirming that Chargaff's second parity rule does not hold at the third codon position, indicating asymmetric mutational and/or selective pressures. The slight G3 bias (mean 0.511) and T3 bias (mean A3/(A3+T3) = 0.460) reflect transcription-coupled mutational asymmetry.

### 4.7 CpG Content Analysis

#### 4.7.1 CpG Observed/Expected Ratio

The CpG O/E ratio across all 50 genes ranged from 0.18 (BRCA1) to 0.93 (VHL), with a mean of 0.44 ± 0.18.

| Category | Mean CpG O/E | CpG Island Criteria Met* |
|----------|--------------|--------------------------|
| Oncogenes | 0.48 ± 0.13 | 5/20 (25%) |
| TSGs | 0.44 ± 0.22 | 4/20 (20%) |
| Dual-Role | 0.39 ± 0.15 | 1/10 (10%) |

*CpG island criteria: CpG O/E ≥ 0.6 and GC% ≥ 50%

**Key Finding:** The proportion of genes meeting CpG island criteria within their coding sequences was similar between oncogenes (25%) and TSGs (20%), and this difference was not statistically significant (Fisher's exact test, OR = 1.33, p = 1.000). Notably, the overall CpG O/E ratios were substantially below 1.0 across all categories (range 0.18–0.93, mean 0.44), confirming widespread CpG depletion consistent with methylation-mediated deamination across cancer gene coding regions.

#### 4.7.2 CpG Dinucleotide Frequency

The CpG depletion across all 50 genes (mean O/E = 0.44) confirms widespread cytosine methylation and deamination within cancer gene coding regions. Genes with the highest CpG retention (VHL: 0.93, CDKN2A: 0.90, MYC: 0.75, WT1: 0.75) tend to be shorter genes with high GC content, while genes with the strongest CpG depletion (BRCA1: 0.18, BRCA2: 0.23, TET2: 0.20, MDM2: 0.20) are typically longer genes with lower GC content.

### 4.8 Correlation Analysis

Pearson correlation analysis revealed several significant relationships:

| Parameter Pair | r | r² | p-value | Interpretation |
|---------------|---|-----|---------|----------------|
| GC% vs. GC3% | 0.958 | 0.917 | 1.42 × 10⁻²⁷ | Near-perfect correlation; GC3 drives overall GC |
| GC% vs. ENC | -0.666 | 0.444 | 1.28 × 10⁻⁷ | Higher GC → stronger codon bias |
| GC% vs. CpG O/E | 0.804 | 0.647 | 2.00 × 10⁻¹² | GC-rich genes retain more CpG |
| GC3% vs. GC12% | 0.688 | 0.474 | 3.29 × 10⁻⁸ | Compositional correlation across positions |
| GC% vs. CDS length | -0.293 | 0.086 | 0.039 | Longer genes tend to be GC-poorer |
| GC3% vs. ENC | -0.701 | 0.492 | 1.43 × 10⁻⁸ | GC3-rich genes show stronger bias |
| GC-skew vs. AT-skew | 0.261 | 0.068 | 0.067 | Weak, non-significant correlation |

The near-perfect correlation between GC% and GC3% (r = 0.958) confirms that third-position composition is the primary determinant of overall GC content in cancer genes. The strong negative correlation between GC3% and ENC (r = -0.701) indicates that genes with higher GC3 content show stronger codon usage bias, consistent with both compositional constraint and translational selection.

### 4.9 Summary of Key Results

| Finding | Description | Significance |
|---------|-------------|--------------|
| GC content trend | Oncogenes (51.84%) > TSGs (48.81%) | p = 0.323 (ns) |
| GC3 trend | Oncogenes (59.97%) > TSGs (50.61%) | p = 0.086 (ns) |
| ENC trend | Oncogenes (46.73) < TSGs (49.54) | p = 0.126 (ns) |
| Neutrality slope | 0.200 (selection dominant, 80%) | p = 3.29 × 10⁻⁸ |
| CpG O/E | Oncogenes (0.48) > TSGs (0.44) | p = 0.156 (ns) |
| Preferred codons | GC-ending codons dominant | CTG (2.18), GTG (1.76) |
| ENC-GC3 deviation | 84% genes below expected curve | Translational selection |
| PR2 deviation | 92% genes deviate from parity | Asymmetric pressures |

---

## 5. Discussion

### 5.1 GC Content Dichotomy Between Oncogenes and Tumor Suppressors

A notable finding of this study is the observed trend toward higher GC content in oncogenes (mean 51.84%) compared to tumor suppressor genes (48.81%), although this difference did not reach statistical significance (p = 0.323). The trend was more pronounced at the third codon position (GC3: 59.97% vs. 50.61%, p = 0.086), approaching significance. The lack of statistical significance with n=20 per group does not negate the biological relevance of the trend, and several biological implications merit discussion:

**Genomic localization:** GC-rich oncogenes tend to reside in gene-dense, GC-rich isochores (H2-H3 isochores) that are associated with open chromatin, early replication timing, and high recombination rates (Bernardi, 2000). This genomic environment may facilitate the ready accessibility of oncogenes for transcription. In contrast, many large TSGs (BRCA1, BRCA2, ATM, APC) reside in more GC-poor regions.

**Expression regulation:** The higher GC content of oncogenes may contribute to their expression potential, as GC-rich sequences form more stable mRNA secondary structures and are more efficiently processed and exported from the nucleus (Kudla et al., 2006). This compositional advantage may facilitate the overexpression of oncogenes in cancer.

**mRNA stability:** Higher GC content confers greater thermodynamic stability to mRNA molecules, potentially increasing the half-life of oncogene transcripts and amplifying their protein output.

### 5.2 Codon Usage Patterns and Translational Efficiency

The trend toward lower ENC values in oncogenes (46.73) compared to TSGs (49.54), though not statistically significant (p = 0.126), is consistent with the hypothesis that oncogenes may experience stronger translational selection. Oncogenes encode proteins often required in large quantities during cell proliferation (growth factors, kinases, transcription factors), and codon optimization may enhance their translational efficiency.

The predominance of GC-ending preferred codons — particularly CTG for Leucine (RSCU = 2.18) and GTG for Valine (RSCU = 1.76) — reflects both compositional bias and translational selection, as these codons are recognized by abundant tRNAs in human cells (dos Reis et al., 2004).

The ENC-GC3s plot analysis revealed that 84% of cancer genes fall below the expected curve, a stronger signal than previously reported. This indicates that translational selection operates in addition to mutational pressure across all cancer gene categories, consistent with the functional importance of cancer genes and their high expression levels in relevant tissues.

### 5.3 Evolutionary Forces Acting on Cancer Genes

The neutrality plot regression slope of 0.200 indicates that natural selection (80.0%) is the dominant force shaping codon usage in cancer genes, substantially stronger than mutation pressure (20.0%). This finding is consistent with previous studies on human autophagy genes (slope = 0.213; Malik et al., 2022) and oral cancer genes (Chakraborty et al., 2019), and suggests even stronger selective constraint on cancer genes than previously reported.

Interestingly, TSGs showed a slightly higher mutation pressure contribution (27.7%, slope = 0.277) compared to oncogenes (20.1%, slope = 0.201), while dual-role genes exhibited the strongest selective constraint (87.0%, slope = 0.130). The lower neutrality slope in dual-role genes may reflect their requirement to function in both oncogenic and tumor-suppressive contexts, imposing stringent constraints on codon usage that limit compositional drift.

### 5.4 CpG Content and Mutational Vulnerability

The CpG O/E ratios were uniformly low across all gene categories (oncogenes: 0.48, TSGs: 0.44, dual-role: 0.39), confirming pervasive CpG depletion in cancer gene coding regions. However, the range was wide (0.18–0.93), revealing substantial gene-level variation with important implications:

1. **Universal CpG depletion:** The mean CpG O/E of 0.44 across all 50 genes is well below 1.0, indicating that methylation-mediated deamination has reduced CpG content substantially, consistent with the general mammalian pattern. Notably, even genes with the highest CpG retention (VHL: 0.93, CDKN2A: 0.90) are still below 1.0.

2. **Gene length–CpG relationship:** Genes with the strongest CpG depletion (BRCA1: 0.18, BRCA2: 0.23, TET2: 0.20) are among the longest in the dataset, while genes with higher CpG retention tend to be shorter and GC-rich. This is consistent with the strong positive correlation between GC% and CpG O/E (r = 0.804, p = 2.00 × 10⁻¹²).

3. **Mutational hotspots:** Despite overall depletion, the residual CpG sites in cancer genes remain mutational hotspots. The BRAF V600E mutation, one of the most common oncogenic mutations, occurs at a C→T transition in a CpG context, illustrating how even depleted CpG sites contribute to cancer-driving mutations.

### 5.5 Implications for Cancer Biology

The compositional differences between oncogenes and TSGs revealed in this study have several translational implications:

**Diagnostic panel design:** Understanding the GC content profiles of cancer genes can inform the design of targeted sequencing panels, as GC-rich regions require specialized library preparation and sequencing strategies to achieve uniform coverage.

**Gene therapy:** Codon optimization strategies for cancer gene therapy (e.g., expressing wild-type TP53 in tumor cells) should account for the natural codon usage preferences of the target gene category.

**Biomarker development:** The CpG content analysis provides information relevant to methylation-based biomarker development, as CpG-rich coding regions may serve as informative methylation targets.

**Mutational signature analysis:** The nucleotide composition and CpG content of cancer genes influence the expected mutational spectrum, which is relevant for distinguishing driver mutations from passenger mutations.

### 5.6 Comparison with Previous Studies

Our findings are broadly consistent with previous investigations:

| Study | Focus | Key Finding | Our Comparison |
|-------|-------|-------------|----------------|
| Prabha et al. (2017) | Cancer gene codon usage | GC3 bias in cancer genes | Confirmed — cancer genes show elevated GC3 (mean 54.7%) |
| Malik et al. (2022) | Autophagy genes | Mean GC=52.23%, selection dominant (slope=0.213) | Similar GC range (49.9%), even stronger selection (slope=0.200) |
| Chakraborty et al. (2019) | Oral cancer genes | ENC range 34.6-55.9 | Our range 33.4-58.1, highly consistent |
| Supek et al. (2014) | Tumor codon usage | tRNA-codon co-adaptation in cancer | Supports our finding of 84% genes below ENC curve |
| D'Onofrio et al. (1999) | Gene GC content | GC varies with isochore location | Confirmed — wide GC range (35.7-72.0%) reflects isochore variation |

### 5.7 Limitations

1. **Sample size and statistical power:** With 20 genes per category, statistical power to detect moderate effect sizes is limited. The trends observed (e.g., GC3: 59.97% vs. 50.61%, p = 0.086) may achieve significance with a larger gene panel from COSMIC CGC Tier 2 or whole-genome cancer driver gene lists.
2. **Gene selection bias:** Analysis was restricted to 50 well-characterized COSMIC Tier 1 cancer genes; inclusion of more recently identified cancer genes from large-scale sequencing studies might reveal additional patterns.
3. **Single isoform analysis:** Only the canonical RefSeq transcript isoform was analyzed for each gene; alternative splicing may produce transcripts with different compositional properties.
4. **Static analysis:** This study analyzes the reference genome sequence; somatic mutations in cancer may alter the local GC content and codon usage of individual genes.
5. **Functional context:** The analysis treats all coding sequence equally; domain-specific compositional variation may be functionally relevant but was not systematically examined.
6. **Cancer type specificity:** The gene classification (oncogene vs. TSG) is simplified; some genes (dual-role) function differently across cancer types, and genes like NOTCH1 can act as oncogenes in some contexts and tumor suppressors in others.

---

## 6. Conclusions

This comprehensive computational analysis of GC content variability and nucleotide composition in 50 human cancer-associated genes reveals:

1. **Compositional trends by gene category:** Oncogenes show a trend toward higher GC content (51.84%) compared to tumor suppressor genes (48.81%), with the divergence most pronounced at the third codon position (GC3: 59.97% vs. 50.61%, p = 0.086). While not reaching conventional statistical significance at n=20 per group, these trends are biologically consistent with isochore localization patterns.

2. **Biased codon usage:** Cancer genes preferentially use specific GC-ending codons, with CTG (Leu, RSCU = 2.18) and GTG (Val, RSCU = 1.76) as the most overrepresented. Eight codons were strongly avoided (RSCU < 0.6), including four CG-ending codons (GCG, TCG, CCG, ACG), reflecting CpG dinucleotide avoidance.

3. **Strong selection-dominated evolution:** Natural selection contributes approximately 80.0% of codon usage variation (neutrality slope = 0.200, p = 3.29 × 10⁻⁸), with mutation pressure contributing only 20.0%. This represents one of the strongest selection signals reported for human gene subsets.

4. **Pervasive CpG depletion:** All 50 cancer genes show CpG O/E ratios below 1.0 (mean 0.44 ± 0.18), confirming widespread methylation-mediated CpG loss. CpG retention correlates strongly with overall GC content (r = 0.804).

5. **Widespread translational selection:** 84% of cancer genes fall below the expected ENC-GC3s curve, indicating that translational selection acts broadly across cancer genes regardless of oncogene/TSG classification.

6. **Asymmetric mutational pressures:** PR2 analysis reveals that 92% of genes deviate from Chargaff's second parity rule at the third codon position, with biases toward T3 (over A3) and G3 (over C3), indicating transcription-coupled mutational asymmetry.

These findings provide a foundational understanding of the compositional landscape of cancer genes and contribute to the broader fields of cancer genomics, molecular evolution, and computational biology. The strong selection signal (80%) and widespread translational selection (84% below ENC curve) highlight that cancer genes, as a class, are under stringent evolutionary constraint. Future studies should extend this analysis to include larger gene panels, tissue-specific expression data, somatic mutation profiles, and epigenomic information to create integrated models of cancer gene regulation.

---

## 7. Future Directions

1. **Integrate with expression data:** Correlate compositional features with RNA-seq expression data across The Cancer Genome Atlas (TCGA) cancer types
2. **Somatic mutation mapping:** Overlay somatic mutation locations from COSMIC onto GC content sliding window profiles to identify compositionally vulnerable regions
3. **Epigenomic integration:** Combine CpG O/E data with genome-wide methylation data (ENCODE, Roadmap Epigenomics) to assess methylation status of CpG-rich coding regions
4. **Machine learning classification:** Develop predictive models that use compositional features to classify cancer driver genes
5. **Pan-cancer analysis:** Extend the analysis to all ~700+ COSMIC CGC genes across both tiers
6. **Comparative genomics:** Compare human cancer gene composition with orthologous genes in cancer-resistant species (e.g., naked mole-rat, elephant)
7. **3D genome context:** Examine the relationship between cancer gene GC content and chromatin organization (TADs, compartments)

---

## 8. References

1. Alexandrov, L.B., Nik-Zainal, S., Wedge, D.C., et al. (2013). Signatures of mutational processes in human cancer. *Nature*, 500(7463), 415-421.

2. Baylin, S.B. & Jones, P.A. (2016). Epigenetic determinants of cancer. *Cold Spring Harbor Perspectives in Biology*, 8(9), a019505.

3. Bernardi, G. (2000). Isochores and the evolutionary genomics of vertebrates. *Gene*, 241(1), 3-17.

4. Buhr, F., Jha, S., Thommen, M., et al. (2016). Synonymous codons direct cotranslational folding toward different protein conformations. *Molecular Cell*, 61(3), 341-351.

5. Chakraborty, S., Nag, D., Mazumder, T.H., et al. (2019). Understanding the nucleotide composition and codon usage patterns of genes expressed in human oral cancer. *Journal of Oral Biosciences*, 62(1), 15-22.

6. Cock, P.J., Antao, T., Chang, J.T., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.

7. Costantini, M., Clay, O., Auletta, F., & Bernardi, G. (2006). An isochore map of human chromosomes. *Genome Research*, 16(4), 536-541.

8. D'Onofrio, G., Jabbari, K., Musto, H., & Bernardi, G. (1999). The correlation of protein hydropathy with the base composition of coding sequences. *Gene*, 238(1), 3-14.

9. dos Reis, M., Savva, R., & Wernisch, L. (2004). Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Research*, 32(17), 5036-5044.

10. Gardiner-Garden, M. & Frommer, M. (1987). CpG islands in vertebrate genomes. *Journal of Molecular Biology*, 196(2), 261-282.

11. Goodarzi, H., Nguyen, H.C.B., Zhang, S., et al. (2016). Modulated expression of specific tRNAs drives gene expression and cancer progression. *Cell*, 165(6), 1416-1427.

12. Hershberg, R. & Petrov, D.A. (2008). Selection on codon bias. *Annual Review of Genetics*, 42, 287-299.

13. Kudla, G., Murray, A.W., Tollervey, D., & Plotkin, J.B. (2006). Coding-sequence determinants of gene expression in *Escherichia coli*. *Science*, 324(5924), 255-258.

14. Kumar, S., Stecher, G., Li, M., Knyaz, C., & Tamura, K. (2018). MEGA X: Molecular Evolutionary Genetics Analysis across computing platforms. *Molecular Biology and Evolution*, 35(6), 1547-1549.

15. Lawrence, M.S., Stojanov, P., Polak, P., et al. (2013). Mutational heterogeneity in cancer and the search for new cancer-associated genes. *Nature*, 499(7457), 214-218.

16. Malik, A., Firoz, A., Jha, V., et al. (2022). Analysis of the compositional features and codon usage pattern of genes involved in human autophagy. *Frontiers in Genetics*, 13, 1019258.

17. Mauger, D.M., Cabral, B.J., Presnyak, V., et al. (2019). mRNA structure regulates protein expression through changes in functional half-life. *Proceedings of the National Academy of Sciences*, 116(48), 24075-24083.

18. Novembre, J.A. (2002). Accounting for background nucleotide composition when measuring codon usage bias. *Molecular Biology and Evolution*, 19(8), 1390-1394.

19. Peden, J.F. (1999). Analysis of codon usage. PhD Thesis, University of Nottingham.

20. Piovesan, A., Pelleri, M.C., Antonaros, F., et al. (2019). On the length, weight and GC content of the human genome. *BMC Research Notes*, 12, 106.

21. Plotkin, J.B. & Kudla, G. (2011). Synonymous but not the same: the causes and consequences of codon bias. *Nature Reviews Genetics*, 12(1), 32-42.

22. Puigbò, P., Bravo, I.G., & Garcia-Vallve, S. (2008). CAIcal: a combined set of tools to assess codon usage adaptation. *Biology Direct*, 3, 38.

23. Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: the European Molecular Biology Open Software Suite. *Trends in Genetics*, 16(6), 276-277.

24. Sharp, P.M. & Li, W.H. (1987). The codon adaptation index—a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281-1295.

25. Sharp, P.M., Emery, L.R., & Zeng, K. (2010). Forces that influence the evolution of codon bias. *Philosophical Transactions of the Royal Society B*, 365(1544), 1203-1212.

26. Sondka, Z., Bamford, S., Cole, C.G., et al. (2018). The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. *Nature Reviews Cancer*, 18(11), 696-705.

27. Sung, H., Ferlay, J., Siegel, R.L., et al. (2021). Global cancer statistics 2020: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. *CA: A Cancer Journal for Clinicians*, 71(3), 209-249.

28. Supek, F., Miñana, B., Valcárcel, J., Gabaldón, T., & Lehner, B. (2014). Synonymous mutations frequently act as driver mutations in human cancers. *Cell*, 156(6), 1324-1335.

29. Vogelstein, B., Papadopoulos, N., Velculescu, V.E., et al. (2013). Cancer genome landscapes. *Science*, 339(6127), 1546-1558.

30. Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23-29.

---

## Appendices

See supplementary materials for:
- **Appendix A:** Complete gene list with NCBI accession numbers
- **Appendix B:** Full RSCU tables for all 50 genes
- **Appendix C:** Individual gene GC content sliding window plots
- **Appendix D:** Python scripts used for analysis
- **Appendix E:** Raw data tables
