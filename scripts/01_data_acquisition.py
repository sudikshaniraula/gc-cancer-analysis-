"""
Script 01: Data Acquisition
============================
Downloads coding sequences (CDS) of 50 cancer-associated genes from NCBI.
Uses BioPython's Entrez module to fetch RefSeq sequences in FASTA format.

Author: Sudiksha
Project: GC Content Variability in Cancer Genes
Date: April 2025 - March 2026
"""

import os
import time
from Bio import Entrez, SeqIO

# ============================================================
# CONFIGURATION
# ============================================================
Entrez.email = "your_email@university.edu"  # REQUIRED: Set your email
Entrez.api_key = None  # Optional: NCBI API key for faster requests

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# CANCER GENE LIST
# Based on COSMIC Cancer Gene Census Tier 1
# ============================================================
CANCER_GENES = {
    # === ONCOGENES ===
    "KRAS":   {"refseq": "NM_004985.5", "category": "oncogene", "cancer_types": "Pancreatic, Colorectal, Lung"},
    "BRAF":   {"refseq": "NM_004333.6", "category": "oncogene", "cancer_types": "Melanoma, Thyroid, Colorectal"},
    "MYC":    {"refseq": "NM_002467.6", "category": "oncogene", "cancer_types": "Burkitt Lymphoma, Breast, Lung"},
    "EGFR":   {"refseq": "NM_005228.5", "category": "oncogene", "cancer_types": "Lung, Glioblastoma, Colorectal"},
    "ERBB2":  {"refseq": "NM_004448.4", "category": "oncogene", "cancer_types": "Breast, Gastric, Ovarian"},
    "PIK3CA": {"refseq": "NM_006218.4", "category": "oncogene", "cancer_types": "Breast, Colorectal, Endometrial"},
    "ABL1":   {"refseq": "NM_005157.6", "category": "oncogene", "cancer_types": "CML, ALL"},
    "RET":    {"refseq": "NM_020975.6", "category": "oncogene", "cancer_types": "Thyroid, MEN2, Lung"},
    "KIT":    {"refseq": "NM_000222.3", "category": "oncogene", "cancer_types": "GIST, Melanoma, AML"},
    "MET":    {"refseq": "NM_000245.4", "category": "oncogene", "cancer_types": "Renal, Lung, Hepatocellular"},
    "ALK":    {"refseq": "NM_004304.5", "category": "oncogene", "cancer_types": "Lung, Neuroblastoma, ALCL"},
    "NRAS":   {"refseq": "NM_002524.5", "category": "oncogene", "cancer_types": "Melanoma, AML, Thyroid"},
    "HRAS":   {"refseq": "NM_005343.4", "category": "oncogene", "cancer_types": "Bladder, Thyroid, Salivary"},
    "FLT3":   {"refseq": "NM_004119.3", "category": "oncogene", "cancer_types": "AML"},
    "JAK2":   {"refseq": "NM_004972.4", "category": "oncogene", "cancer_types": "MPN, AML, ALL"},
    "FGFR1":  {"refseq": "NM_023110.3", "category": "oncogene", "cancer_types": "MPN, Breast, Lung"},
    "FGFR2":  {"refseq": "NM_000141.5", "category": "oncogene", "cancer_types": "Endometrial, Breast, Gastric"},
    "FGFR3":  {"refseq": "NM_000142.5", "category": "oncogene", "cancer_types": "Bladder, Myeloma, Cervical"},
    "CDK4":   {"refseq": "NM_000075.4", "category": "oncogene", "cancer_types": "Melanoma, Sarcoma, Glioblastoma"},
    "CCND1":  {"refseq": "NM_053056.3", "category": "oncogene", "cancer_types": "Breast, Mantle Cell Lymphoma"},

    # === TUMOR SUPPRESSOR GENES ===
    "TP53":   {"refseq": "NM_000546.6", "category": "tsg", "cancer_types": "Pan-cancer (>50% of all cancers)"},
    "BRCA1":  {"refseq": "NM_007294.4", "category": "tsg", "cancer_types": "Breast, Ovarian"},
    "BRCA2":  {"refseq": "NM_000059.4", "category": "tsg", "cancer_types": "Breast, Ovarian, Pancreatic"},
    "RB1":    {"refseq": "NM_000321.3", "category": "tsg", "cancer_types": "Retinoblastoma, Osteosarcoma, SCLC"},
    "APC":    {"refseq": "NM_000038.6", "category": "tsg", "cancer_types": "Colorectal (FAP)"},
    "PTEN":   {"refseq": "NM_000314.8", "category": "tsg", "cancer_types": "Endometrial, Glioblastoma, Prostate"},
    "VHL":    {"refseq": "NM_000551.4", "category": "tsg", "cancer_types": "Renal Cell Carcinoma"},
    "WT1":    {"refseq": "NM_024426.6", "category": "tsg", "cancer_types": "Wilms Tumor, AML"},
    "NF1":    {"refseq": "NM_001042492.3", "category": "tsg", "cancer_types": "Neurofibromatosis, MPNST, Glioma"},
    "NF2":    {"refseq": "NM_000268.4", "category": "tsg", "cancer_types": "Meningioma, Schwannoma"},
    "CDKN2A": {"refseq": "NM_000077.5", "category": "tsg", "cancer_types": "Melanoma, Pancreatic, Glioblastoma"},
    "SMAD4":  {"refseq": "NM_005359.6", "category": "tsg", "cancer_types": "Pancreatic, Colorectal"},
    "STK11":  {"refseq": "NM_000455.5", "category": "tsg", "cancer_types": "Lung, Peutz-Jeghers Syndrome"},
    "MLH1":   {"refseq": "NM_000249.4", "category": "tsg", "cancer_types": "Colorectal (Lynch Syndrome)"},
    "MSH2":   {"refseq": "NM_000251.3", "category": "tsg", "cancer_types": "Colorectal (Lynch Syndrome)"},
    "ATM":    {"refseq": "NM_000051.4", "category": "tsg", "cancer_types": "Breast, Lymphoma, Leukemia"},
    "BAP1":   {"refseq": "NM_004656.4", "category": "tsg", "cancer_types": "Mesothelioma, Uveal Melanoma, Renal"},
    "CHEK2":  {"refseq": "NM_007194.4", "category": "tsg", "cancer_types": "Breast, Prostate, Colorectal"},
    "CDH1":   {"refseq": "NM_004360.5", "category": "tsg", "cancer_types": "Gastric, Lobular Breast"},
    "FBXW7":  {"refseq": "NM_033632.3", "category": "tsg", "cancer_types": "T-ALL, Colorectal, Endometrial"},

    # === DUAL-ROLE / OTHER CANCER GENES ===
    "NOTCH1": {"refseq": "NM_017617.5", "category": "dual", "cancer_types": "T-ALL, HNSCC, CLL"},
    "IDH1":   {"refseq": "NM_005896.4", "category": "dual", "cancer_types": "Glioma, AML, Cholangiocarcinoma"},
    "IDH2":   {"refseq": "NM_002168.4", "category": "dual", "cancer_types": "AML, Glioma, Angioimmunoblastic"},
    "EZH2":   {"refseq": "NM_004456.5", "category": "dual", "cancer_types": "DLBCL, Follicular Lymphoma"},
    "DNMT3A": {"refseq": "NM_175629.2", "category": "dual", "cancer_types": "AML, MDS, T-cell Lymphoma"},
    "TET2":   {"refseq": "NM_001127208.3", "category": "dual", "cancer_types": "AML, MDS, MPN"},
    "SF3B1":  {"refseq": "NM_012433.4", "category": "dual", "cancer_types": "MDS, CLL, Uveal Melanoma"},
    "NPM1":   {"refseq": "NM_002520.7", "category": "dual", "cancer_types": "AML, NHL"},
    "CTNNB1": {"refseq": "NM_001904.4", "category": "dual", "cancer_types": "Colorectal, Hepatoblastoma, Desmoid"},
    "MDM2":   {"refseq": "NM_002392.6", "category": "dual", "cancer_types": "Sarcoma, Glioblastoma"},
}


def fetch_gene_sequence(gene_name, refseq_id):
    """Fetch CDS from NCBI Nucleotide database using RefSeq accession."""
    try:
        # Search for the RefSeq record
        handle = Entrez.efetch(
            db="nucleotide",
            id=refseq_id,
            rettype="fasta_cds_na",
            retmode="text"
        )
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        if records:
            # Take the first (longest) CDS
            record = max(records, key=lambda r: len(r.seq))
            return record
        else:
            print(f"  WARNING: No CDS found for {gene_name} ({refseq_id})")
            return None

    except Exception as e:
        print(f"  ERROR fetching {gene_name}: {e}")
        return None


def validate_cds(sequence):
    """Validate that the sequence is a proper CDS."""
    seq_str = str(sequence).upper()
    issues = []

    # Check start codon
    if not seq_str.startswith("ATG"):
        issues.append("Missing start codon (ATG)")

    # Check stop codon
    if seq_str[-3:] not in ("TAA", "TAG", "TGA"):
        issues.append("Missing stop codon")

    # Check length divisible by 3
    if len(seq_str) % 3 != 0:
        issues.append(f"Length ({len(seq_str)}) not divisible by 3")

    # Check for internal stop codons
    for i in range(3, len(seq_str) - 3, 3):
        codon = seq_str[i:i+3]
        if codon in ("TAA", "TAG", "TGA"):
            issues.append(f"Internal stop codon ({codon}) at position {i}")
            break

    return issues


def main():
    """Main function to download and validate all gene sequences."""
    print("=" * 70)
    print("CANCER GENE SEQUENCE ACQUISITION")
    print("=" * 70)
    print(f"Total genes to fetch: {len(CANCER_GENES)}")
    print(f"Output directory: {OUTPUT_DIR}")
    print()

    # Track results
    successful = []
    failed = []
    validation_warnings = []

    # Combined output file
    combined_file = os.path.join(OUTPUT_DIR, "all_cancer_genes_cds.fasta")

    all_records = []

    for i, (gene_name, info) in enumerate(CANCER_GENES.items(), 1):
        print(f"[{i:02d}/{len(CANCER_GENES)}] Fetching {gene_name} ({info['refseq']})...", end=" ")

        record = fetch_gene_sequence(gene_name, info["refseq"])

        if record:
            # Validate CDS
            issues = validate_cds(record.seq)

            # Update record description with gene info
            record.id = gene_name
            record.name = gene_name
            record.description = (
                f"{gene_name} | {info['category'].upper()} | "
                f"RefSeq:{info['refseq']} | "
                f"Length:{len(record.seq)}bp | "
                f"Cancer:{info['cancer_types']}"
            )

            # Save individual FASTA
            individual_file = os.path.join(OUTPUT_DIR, f"{gene_name}_cds.fasta")
            SeqIO.write(record, individual_file, "fasta")

            all_records.append(record)
            successful.append(gene_name)

            if issues:
                validation_warnings.append((gene_name, issues))
                print(f"OK (with warnings: {'; '.join(issues)})")
            else:
                print(f"OK ({len(record.seq)} bp)")
        else:
            failed.append(gene_name)
            print("FAILED")

        # Rate limiting: NCBI allows 3 requests/second without API key
        time.sleep(0.4)

    # Write combined FASTA file
    if all_records:
        SeqIO.write(all_records, combined_file, "fasta")

    # Write gene metadata file
    metadata_file = os.path.join(OUTPUT_DIR, "gene_metadata.tsv")
    with open(metadata_file, "w") as f:
        f.write("Gene\tCategory\tRefSeq\tCDS_Length\tCancer_Types\n")
        for record in all_records:
            gene = record.id
            info = CANCER_GENES[gene]
            f.write(f"{gene}\t{info['category']}\t{info['refseq']}\t"
                    f"{len(record.seq)}\t{info['cancer_types']}\n")

    # Print summary
    print()
    print("=" * 70)
    print("ACQUISITION SUMMARY")
    print("=" * 70)
    print(f"Successfully fetched: {len(successful)}/{len(CANCER_GENES)}")
    print(f"Failed: {len(failed)}")
    if failed:
        print(f"  Failed genes: {', '.join(failed)}")
    if validation_warnings:
        print(f"Validation warnings: {len(validation_warnings)}")
        for gene, issues in validation_warnings:
            print(f"  {gene}: {'; '.join(issues)}")
    print(f"\nOutput files:")
    print(f"  Combined FASTA: {combined_file}")
    print(f"  Metadata TSV:   {metadata_file}")
    print(f"  Individual CDS: {OUTPUT_DIR}/<GENE>_cds.fasta")


if __name__ == "__main__":
    main()
