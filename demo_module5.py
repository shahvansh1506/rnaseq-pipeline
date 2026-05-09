"""
demo_module5.py — Demo for Module 5: DESeq2 Differential Expression.

Run from your project root:
    python demo_module5.py

This demo:
1. Simulates a realistic count matrix (like GSE60450)
2. Builds sample metadata
3. Runs DESeq2 via the R wrapper
4. Prints a summary of significant genes
5. Shows where plots are saved

NOTE: This demo actually runs DESeq2 in R!
Make sure you have run install_r_packages.R first:
    Rscript install_r_packages.R
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from src.parsers.count_parser import CountParser
from src.analysis.run_deseq2 import DESeq2Runner


# ══════════════════════════════════════════════════════════════
# SIMULATE REALISTIC COUNT DATA
# ══════════════════════════════════════════════════════════════

def simulate_rna_seq_counts() -> tuple:
    """
    Simulate a realistic RNA-seq count matrix similar to GSE60450.

    Creates:
    - 6 samples: basal_virgin x2, basal_pregnant x2, luminal_virgin x2
    - 1000 genes with realistic count distributions
    - ~100 genuinely DE genes between conditions

    Returns:
        (counts DataFrame, metadata DataFrame)
    """
    rng = np.random.default_rng(42)

    samples = [
        "basal_virgin_1", "basal_virgin_2",
        "basal_pregnant_1", "basal_pregnant_2",
        "luminal_virgin_1", "luminal_virgin_2",
    ]

    n_genes = 1000
    genes   = [f"Gene_{i:04d}" for i in range(n_genes)]

    # Replace first few with real gene names from the tutorial
    real_genes = ["Csn1s2a", "Csn1s1", "Csn2", "Glycam1",
                  "Wap", "Trf", "Eef1a1", "Actb", "Gapdh", "Brca1"]
    genes[:len(real_genes)] = real_genes

    # Base expression — negative binomial (realistic for RNA-seq)
    counts = rng.negative_binomial(5, 0.3, size=(n_genes, len(samples)))

    # Add housekeeping genes (high, stable expression)
    # Actb, Gapdh, Eef1a1 should NOT be DE
    counts[6:9] = rng.integers(5000, 15000, size=(3, len(samples)))

    # Add DE genes — upregulated in basal_pregnant vs basal_virgin
    # Csn genes are milk proteins, highly expressed in pregnancy
    de_genes_up = list(range(0, 5))    # Csn1s2a, Csn1s1, Csn2, Glycam1, Wap
    for i in de_genes_up:
        counts[i, 0:2] = rng.integers(100, 500,   size=2)   # low in virgin
        counts[i, 2:4] = rng.integers(5000, 20000, size=2)   # high in pregnant
        counts[i, 4:6] = rng.integers(50, 200,    size=2)   # low in luminal

    # Add some randomly DE genes
    random_de = range(20, 120)
    for i in random_de:
        if rng.random() > 0.5:
            counts[i, :3]  *= rng.integers(3, 8)   # up in first 3 samples
        else:
            counts[i, 3:]  *= rng.integers(3, 8)   # up in last 3 samples

    # Ensure all counts are non-negative integers
    counts = np.maximum(counts, 0).astype(int)

    count_df = pd.DataFrame(counts, index=genes, columns=samples)

    # Build metadata
    metadata = pd.DataFrame({
        "condition": [
            "basal_virgin", "basal_virgin",
            "basal_pregnant", "basal_pregnant",
            "luminal_virgin", "luminal_virgin",
        ]
    }, index=samples)

    return count_df, metadata


# ══════════════════════════════════════════════════════════════
# MAIN DEMO
# ══════════════════════════════════════════════════════════════

def main():
    print("\n" + "="*60)
    print("  MODULE 5 DEMO — DESeq2 Differential Expression")
    print("="*60)

    # ── 1. Simulate count data ────────────────────────────────
    print("\n📍 Step 1: Simulating count matrix (1000 genes, 6 samples)")
    counts, metadata = simulate_rna_seq_counts()
    print(f"   Count matrix: {counts.shape[0]} genes x {counts.shape[1]} samples")
    print(f"   Conditions  : {metadata['condition'].unique().tolist()}")

    # ── 2. Filter low count genes ─────────────────────────────
    print("\n📍 Step 2: Filtering low-count genes")
    parser = CountParser(min_total_count=10)
    counts = parser._filter_low_counts(counts)
    print(f"   Kept {len(counts):,} genes after filtering")

    # ── 3. Run DESeq2 ─────────────────────────────────────────
    print("\n📍 Step 3: Running DESeq2 analysis")
    print("   (This requires R and DESeq2 to be installed)")
    print("   If not installed yet, run: Rscript install_r_packages.R\n")

    output_dir = Path("results/deseq2")
    runner     = DESeq2Runner(output_dir=output_dir)

    try:
        results = runner.run(
            counts        = counts,
            metadata      = metadata,
            condition_col = "condition",
        )

        # ── 4. Show results summary ───────────────────────────
        print("\n📍 Step 4: Results Summary")
        runner.print_summary(results)

        # ── 5. Show top significant genes ─────────────────────
        print("📍 Step 5: Top 10 Most Significant Genes")
        sig = runner.get_significant_genes(results)

        if len(sig) > 0:
            top10 = sig.nsmallest(10, "padj")[
                ["gene", "log2FoldChange", "padj",
                 "direction", "comparison"]
            ]
            print(top10.to_string(index=False))
        else:
            print("   No significant genes found with current thresholds.")

        # ── 6. Show output files ──────────────────────────────
        print(f"\n📍 Step 6: Output Files Saved to {output_dir}/")
        for f in ["de_results.csv", "volcano_plot.png",
                  "pca_plot.png", "heatmap_top50.png"]:
            fpath  = output_dir / f
            status = "✅" if fpath.exists() else "❌ not found"
            print(f"   {status} {f}")

        print("\n✅ Module 5 complete!\n")

    except RuntimeError as e:
        print(f"\n⚠️  DESeq2 could not run: {e}")
        print("\nTo fix this, run in Anaconda Prompt:")
        print("    Rscript install_r_packages.R")
        print("\nThen run this demo again.")
        print("\nThe Python tests still work without R — run:")
        print("    pytest tests/test_deseq2.py -v\n")


if __name__ == "__main__":
    main()
