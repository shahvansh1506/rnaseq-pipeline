"""
demo_module3.py — Demo script for Module 3: Quality Control.

Run from your project root:
    python demo_module3.py

This script:
1. Creates demo FASTQ files with mixed quality reads
2. Runs the Python fallback trimmer on them
3. Compares stats before and after trimming
4. Simulates parsing a FastQC report
5. Shows QC summary table
"""

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from src.parsers.fastq_parser import get_fastq_stats, print_stats
from src.qc.trimmer import Trimmer
from src.qc.quality_control import QualityController, QCReport


# ══════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════

def create_mixed_quality_fastq(path: Path, num_reads: int = 20) -> Path:
    """
    Create a FASTQ file with a mix of good and bad reads.
    - Some reads have great quality throughout
    - Some reads have quality that drops at the 3' end
    - Some reads are very short
    - Some reads contain adapter sequences
    """
    import random
    random.seed(42)

    bases   = ["A", "T", "G", "C"]
    adapter = "AGATCGGAAGAGCACACGTCTGAAC"   # TruSeq adapter

    with gzip.open(path, "wt") as f:
        for i in range(num_reads):
            read_type = i % 4

            if read_type == 0:
                # Good quality read — 75bp, all high quality
                seq  = "".join(random.choices(bases, k=75))
                qual = "I" * 75

            elif read_type == 1:
                # Quality degrades at 3' end — common in Illumina
                seq  = "".join(random.choices(bases, k=75))
                qual = "I" * 50 + "+" * 15 + "!" * 10

            elif read_type == 2:
                # Contains adapter contamination
                seq  = "".join(random.choices(bases, k=40)) + adapter[:20]
                qual = "I" * len(seq)

            else:
                # Short read — only 30bp
                seq  = "".join(random.choices(bases, k=30))
                qual = "I" * 30

            f.write(f"@sample_read_{i+1} type={read_type}\n")
            f.write(f"{seq}\n")
            f.write(f"+\n")
            f.write(f"{qual}\n")

    return path


def create_fake_qc_reports() -> list:
    """Create fake QCReport objects to demonstrate the summary table."""
    qc = QualityController.__new__(QualityController)
    qc.reports = []

    fake_data = [
        ("Sample_basal_virgin_1",   1_250_000, 52, True,  []),
        ("Sample_basal_pregnant_1", 1_180_000, 49, True,  ["WARN: Per sequence GC content"]),
        ("Sample_luminal_virgin_1", 980_000,   55, True,  []),
        ("Sample_low_quality_1",    750_000,   71, False, ["FAIL: Per base sequence quality",
                                                            "WARN: Adapter Content"]),
    ]

    for name, reads, gc, passes, warnings in fake_data:
        report = QCReport(
            sample_name      = name,
            total_sequences  = reads,
            gc_percent       = gc,
            passes_qc        = passes,
            warnings         = warnings,
            per_base_quality = "PASS" if passes else "FAIL",
            adapter_content  = "PASS",
        )
        qc.reports.append(report)

    return qc


# ══════════════════════════════════════════════════════════════
# MAIN DEMO
# ══════════════════════════════════════════════════════════════

def main():
    print("\n" + "="*60)
    print("  MODULE 3 DEMO — Quality Control + Trimming")
    print("="*60)

    # ── 1. Create demo files ──────────────────────────────────
    raw_file     = Path("demo_raw.fastq.gz")
    trimmed_file = Path("demo_trimmed.fastq.gz")
    create_mixed_quality_fastq(raw_file, num_reads=20)
    print(f"\n✅ Created demo FASTQ with 20 mixed-quality reads")

    # ── 2. Stats BEFORE trimming ──────────────────────────────
    print("\n📊 Stats BEFORE trimming:")
    stats_before = get_fastq_stats(raw_file, min_quality=20, min_length=50)
    print_stats(stats_before)

    # ── 3. Run trimmer ────────────────────────────────────────
    print("✂️  Running Python trimmer...")
    trimmer = Trimmer(
        output_dir  = Path("."),
        min_quality = 20,
        min_length  = 50,
    )
    result = trimmer.trim_python_fallback(raw_file, trimmed_file)
    print(f"\n✅ Trimming complete!")
    print(f"   {result.summary()}")

    # ── 4. Stats AFTER trimming ───────────────────────────────
    print("\n📊 Stats AFTER trimming:")
    stats_after = get_fastq_stats(trimmed_file, min_quality=20, min_length=50)
    print_stats(stats_after)

    # ── 5. Before vs After comparison ────────────────────────
    print("📈 Before vs After comparison:")
    print(f"   Reads kept  : {stats_before['total_reads']} → {stats_after['total_reads']}")
    print(f"   Mean quality: {stats_before['mean_quality']} → {stats_after['mean_quality']}")
    print(f"   Pass filter : {stats_before['pass_filter_pct']}% → {stats_after['pass_filter_pct']}%")

    # ── 6. QC summary table ───────────────────────────────────
    print("\n📋 Simulated QC Summary Table (4 samples):")
    qc = create_fake_qc_reports()
    qc.output_dir = Path(".")
    qc.print_summary()

    # ── 7. Clean up ───────────────────────────────────────────
    raw_file.unlink(missing_ok=True)
    trimmed_file.unlink(missing_ok=True)
    print("✅ Module 3 demo complete!\n")


if __name__ == "__main__":
    main()
