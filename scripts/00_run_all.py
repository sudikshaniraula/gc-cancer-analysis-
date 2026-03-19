"""
Master Script: Run Complete Analysis Pipeline
===============================================
Executes all analysis scripts in sequence.

Usage:
    python 00_run_all.py

Prerequisites:
    pip install -r requirements.txt
    Set your NCBI email in 01_data_acquisition.py

Author: Sudiksha
Project: GC Content Variability in Cancer Genes
"""

import subprocess
import sys
import os
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

SCRIPTS = [
    ("01_data_acquisition.py", "Downloading gene sequences from NCBI"),
    ("02_gc_content_analysis.py", "Computing GC content and nucleotide composition"),
    ("03_codon_usage_analysis.py", "Analyzing codon usage (RSCU, ENC, CAI)"),
    ("04_visualization.py", "Generating all figures"),
    ("05_statistical_tests.py", "Running statistical analyses"),
]


def run_script(script_name, description):
    """Run a single analysis script."""
    script_path = os.path.join(SCRIPT_DIR, script_name)

    if not os.path.exists(script_path):
        print(f"  SKIPPED: {script_name} not found")
        return False

    print(f"\n{'='*60}")
    print(f"Running: {script_name}")
    print(f"Purpose: {description}")
    print(f"{'='*60}")

    start_time = time.time()

    try:
        result = subprocess.run(
            [sys.executable, script_path],
            capture_output=True, text=True, timeout=600
        )

        elapsed = time.time() - start_time

        if result.returncode == 0:
            print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)
            print(f"\nCompleted in {elapsed:.1f} seconds")
            return True
        else:
            print(f"ERROR in {script_name}:")
            print(result.stderr[-500:] if len(result.stderr) > 500 else result.stderr)
            return False

    except subprocess.TimeoutExpired:
        print(f"TIMEOUT: {script_name} exceeded 10 minutes")
        return False
    except Exception as e:
        print(f"EXCEPTION: {e}")
        return False


def main():
    print("=" * 60)
    print("CANCER GENE COMPOSITIONAL ANALYSIS - FULL PIPELINE")
    print("=" * 60)
    print(f"Scripts to run: {len(SCRIPTS)}")
    print(f"Script directory: {SCRIPT_DIR}")

    results = {}
    total_start = time.time()

    for script_name, description in SCRIPTS:
        success = run_script(script_name, description)
        results[script_name] = "SUCCESS" if success else "FAILED"

    # Summary
    total_time = time.time() - total_start
    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    for script, status in results.items():
        icon = "[OK]" if status == "SUCCESS" else "[FAIL]"
        print(f"  {icon} {script}: {status}")
    print(f"\nTotal time: {total_time:.1f} seconds")


if __name__ == "__main__":
    main()
