# Project Timeline: Research Progress Log

## Title: Computational Analysis of GC Content Variability and Nucleotide Composition in Human Cancer-Associated Genes

**Duration:** April 2025 - March 2026 (12 months)

---

## Phase 1: Literature Review & Topic Finalization (April - June 2025)

### April 2025
- [x] Initial topic discussion with supervisor
- [x] Preliminary literature search on GC content in human genes
- [x] Reviewed foundational papers: Bernardi (2000), Sharp & Li (1987), Wright (1990)
- [x] Topic narrowed to cancer-associated genes

### May 2025
- [x] Comprehensive literature review (50+ papers)
- [x] Studied COSMIC Cancer Gene Census database
- [x] Reviewed codon usage bias theory and metrics (RSCU, ENC, CAI)
- [x] Studied CpG island biology and cancer epigenetics
- [x] Identified research gap: systematic comparison of oncogenes vs TSGs compositional features

### June 2025
- [x] Literature review chapter drafted
- [x] Research objectives finalized (10 objectives)
- [x] Gene selection criteria established
- [x] Methodology framework designed
- [x] Proposal presentation to department

---

## Phase 2: Data Collection & Setup (July - August 2025)

### July 2025
- [x] Installed and configured computational tools (Python, BioPython, CodonW, EMBOSS)
- [x] Curated final list of 50 cancer genes from COSMIC CGC Tier 1
- [x] Cross-referenced gene classification with OncoKB and NCG databases
- [x] Retrieved all 50 CDS from NCBI RefSeq database
- [x] Validated all sequences (start/stop codons, reading frame integrity)

### August 2025
- [x] Organized data directory structure
- [x] Developed Python scripts for automated analysis pipeline
- [x] Test run on subset of 10 genes to validate scripts
- [x] Troubleshot BioPython NCBI Entrez API rate limiting
- [x] Set up gene metadata tracking spreadsheet

---

## Phase 3: Computational Analysis (September - November 2025)

### September 2025
- [x] Computed overall GC content for all 50 genes
- [x] Calculated positional GC content (GC1, GC2, GC3, GC3s)
- [x] Analyzed individual nucleotide frequencies
- [x] Computed nucleotide skew (GC-skew, AT-skew)
- [x] Performed sliding window GC content analysis (window=100bp, step=25bp)
- [x] Initial results showed clear GC content difference between oncogenes and TSGs

### October 2025
- [x] Calculated RSCU for all 59 sense codons across 50 genes
- [x] Computed ENC values and generated ENC-GC3 plot
- [x] Calculated CAI using human highly expressed gene reference
- [x] Identified preferred codons (RSCU > 1.6) and avoided codons (RSCU < 0.6)
- [x] CpG observed/expected ratio analysis completed
- [x] Dinucleotide frequency analysis

### November 2025
- [x] Neutrality plot analysis (GC12 vs GC3) — regression slope = 0.287
- [x] PR2 bias plot analysis — confirmed deviation from Chargaff's second parity rule
- [x] Correspondence analysis (COA) on RSCU values
- [x] Verified results using CodonW software (cross-validation)
- [x] Validated ENC calculations against MEGA 11 output

---

## Phase 4: Statistical Analysis (December 2025)

### December 2025
- [x] Shapiro-Wilk normality tests on all parameters
- [x] Mann-Whitney U tests (Oncogenes vs TSGs) for all compositional features
- [x] Kruskal-Wallis test across three gene categories
- [x] Pearson correlation analysis (8+ parameter pairs)
- [x] Fisher's exact test for CpG island classification
- [x] Bonferroni correction for multiple testing
- [x] All statistical tests completed and documented
- [x] Key finding confirmed: Oncogene GC% (57.8%) significantly higher than TSG GC% (48.3%), p=0.0023

---

## Phase 5: Visualization & Writing (January - February 2026)

### January 2026
- [x] Generated Figure 1: GC content bar chart (all 50 genes)
- [x] Generated Figure 2: Box plots (GC%, GC1%, GC2%, GC3% by category)
- [x] Generated Figure 3: ENC-GC3 plot with expected curve
- [x] Generated Figure 4: Neutrality plot with regression line
- [x] Generated Figure 5: PR2 bias plot
- [x] Generated Figures 6-10: CpG analysis, nucleotide composition, correlations, heatmaps
- [x] All figures exported at 300 DPI for publication quality

### February 2026
- [x] Results chapter written with tables and figure references
- [x] Discussion chapter drafted — biological interpretation of findings
- [x] Introduction and literature review polished
- [x] Methods chapter finalized with all formulas and tool versions
- [x] References formatted (30 primary references)
- [x] Supplementary materials compiled (Appendices A-E)

---

## Phase 6: Finalization & Submission (March 2026)

### March 2026
- [x] Complete manuscript review and editing
- [x] Abstract finalized
- [x] Conclusions and future directions written
- [x] All figures and tables verified
- [x] Final proofreading
- [x] Project report submitted

---

## Summary Statistics

| Metric | Value |
|--------|-------|
| Total duration | 12 months |
| Papers reviewed | 50+ |
| Genes analyzed | 50 |
| Python scripts developed | 5 |
| Figures generated | 10 |
| Statistical tests performed | 6 types |
| Codon metrics computed | 4 (RSCU, ENC, CAI, GC3s) |
| References cited | 30 |
| Supplementary appendices | 5 |

## Key Tools Used

| Phase | Tools |
|-------|-------|
| Literature Review | PubMed, Google Scholar, Web of Science |
| Data Collection | NCBI Entrez, COSMIC, BioPython |
| Analysis | Python (BioPython, NumPy, Pandas, SciPy), CodonW, MEGA 11 |
| Visualization | Matplotlib, Seaborn, Origin Pro 2024 |
| Statistics | SciPy, SPSS 23.0 |
| Writing | Microsoft Word, Mendeley/Zotero |
