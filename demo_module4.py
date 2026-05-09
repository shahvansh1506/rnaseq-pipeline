"""
demo_module4.py — Demo for Module 4: Alignment & Quantification.

Run from your project root:
    python demo_module4.py

Since HISAT2 and featureCounts are Linux/Mac only tools,
this demo simulates their outputs and focuses on:
1. Parsing simulated alignment statistics
2. Building and filtering a count matrix
3. Building sample metadata
4. Normalizing counts for visualization
5. Validating counts + metadata
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from src.alignment.aligner import Aligner, AlignResult
from src.parsers.count_parser import CountParser


# ══════════════════════════════════════════════════════════════
# SIMULATE ALIGNMENT STATS
# ══════════════════════════════════════════════════════════════

def simulate_alignment_stats() -> list:
    """Simulate HISAT2 alignment statistics for 6 samples."""

    # Simulated HISAT2 stderr output for each sample
    fake_stderr_template = """
{total} reads; of these:
  {total} (100.00%) were unpaired; of these:
    {unaligned} ({unalign_pct:.2f}%) aligned 0 times
    {once} ({once_pct:.2f}%) aligned exactly 1 time
    {multi} ({multi_pct:.2f}%) aligned >1 times
{rate:.2f}% overall alignment rate
"""
    samples = [
        ("basal_virgin_1",    1_250_000, 95.2),
        ("basal_virgin_2",    1_180_000, 94.8),
        ("basal_pregnant_1",  1_310_000, 96.1),
        ("basal_pregnant_2",  1_290_000, 95.7),
        ("luminal_virgin_1",  980_000,   93.4),
        ("luminal_virgin_2",  1_050_000, 94.0),
    ]

    aligner = Aligner.__new__(Aligner)
    results = []

    for name, total, rate in samples:
        unaligned   = int(total * (100 - rate) / 100)
        aligned     = total - unaligned
        once        = int(aligned * 0.90)
        multi       = aligned - once

        stderr = fake_stderr_template.format(
            total       = total,
            unaligned   = unaligned,
            unalign_pct = (100 - rate),
            once        = once,
            once_pct    = (once / total * 100),
            multi       = multi,
            multi_pct   = (multi / total * 100),
            rate        = rate,
        )

        result = aligner._parse_hisat2_stats(stderr, name, f"{name}_sorted.bam")
        results.append(result)

    return results


# ══════════════════════════════════════════════════════════════
# SIMULATE COUNT MATRIX
# ══════════════════════════════════════════════════════════════

def simulate_count_matrix(tmp_path: Path) -> Path:
    """
    Simulate a realistic RNA-seq count matrix.
    Uses the same GSE60450 study structure you've been working with.
    """
    rng = np.random.default_rng(42)

    samples = [
        "basal_virgin_1", "basal_virgin_2",
        "basal_pregnant_1", "basal_pregnant_2",
        "luminal_virgin_1", "luminal_virgin_2",
    ]

    n_genes = 500

    # Most genes have low/zero counts (realistic)
    # A few genes have very high counts (housekeeping genes)
    base_counts = rng.negative_binomial(5, 0.3, size=(n_genes, len(samples)))

    # Add some highly expressed genes (top 20)
    base_counts[:20] = rng.integers(5000, 50000, size=(20, len(samples)))

    # Add some zero-count genes (bottom 50) — will be filtered
    base_counts[-50:] = 0

    genes = [f"Gene_{i:04d}" for i in range(n_genes)]

    # Add some real gene names for realism
    real_genes = ["Csn1s2a", "Csn1s1", "Csn2", "Glycam1",
                  "Wap", "Trf", "Eef1a1", "Actb"]
    genes[:len(real_genes)] = real_genes

    df = pd.DataFrame(base_counts, index=genes, columns=samples)

    csv_path = tmp_path / "count_matrix.csv"
    df.to_csv(csv_path)
    return csv_path


# ══════════════════════════════════════════════════════════════
# MAIN DEMO
# ══════════════════════════════════════════════════════════════

def main():
    import tempfile
    tmp = Path(tempfile.mkdtemp())

    print("\n" + "="*60)
    print("  MODULE 4 DEMO — Alignment & Quantification")
    print("="*60)

    # ── 1. Alignment statistics ───────────────────────────────
    print("\n📍 Step 1: Simulated Alignment Statistics")
    aligner = Aligner.__new__(Aligner)
    results = simulate_alignment_stats()
    aligner.print_summary(results)

    # ── 2. Load count matrix ──────────────────────────────────
    print("📍 Step 2: Loading & Filtering Count Matrix")
    csv_path = simulate_count_matrix(tmp)
    parser   = CountParser(min_total_count=10)
    counts   = parser.load_and_filter(csv_path)
    print(f"   Loaded count matrix: {counts.shape[0]} genes x "
          f"{counts.shape[1]} samples")

    # ── 3. Build metadata ─────────────────────────────────────
    print("\n📍 Step 3: Building Sample Metadata")
    condition_map = {
        "basal_virgin_1"  : "basal_virgin",
        "basal_virgin_2"  : "basal_virgin",
        "basal_pregnant_1": "basal_pregnant",
        "basal_pregnant_2": "basal_pregnant",
        "luminal_virgin_1": "luminal_virgin",
        "luminal_virgin_2": "luminal_virgin",
    }
    metadata = parser.build_metadata(counts, condition_map)
    print(f"   Built metadata for {len(metadata)} samples")
    print(f"   Conditions: {metadata['condition'].unique().tolist()}")

    # ── 4. Validate ───────────────────────────────────────────
    print("\n📍 Step 4: Validating Count Matrix + Metadata")
    parser.validate(counts, metadata)
    print("   Validation passed!")

    # ── 5. Normalize ──────────────────────────────────────────
    print("\n📍 Step 5: Normalizing Counts")
    log_cpm = parser.normalize_log_cpm(counts)
    print(f"   log2(CPM+1) range: "
          f"{log_cpm.values.min():.2f} - {log_cpm.values.max():.2f}")

    # ── 6. Full summary ───────────────────────────────────────
    print("\n📍 Step 6: Count Matrix Summary")
    parser.print_summary(counts, metadata)

    # ── 7. Save count matrix ──────────────────────────────────
    out_path = tmp / "filtered_counts.csv"
    counts.to_csv(out_path)
    print(f"   Count matrix saved to: {out_path}")

    print("✅ Module 4 demo complete!\n")


if __name__ == "__main__":
    main()
