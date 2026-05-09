"""
demo_module2.py — Quick demo to test Module 2 parsers.

Run from your project root:
    python demo_module2.py

This script:
1. Creates a tiny fake FASTQ file
2. Parses it using our FastqRecord generator
3. Prints per-read stats
4. Prints overall QC summary
5. (Optional) Fetches real GEO metadata — uncomment to try!
"""

import gzip
import tempfile
from pathlib import Path

# Add project root to path so imports work
import sys
sys.path.insert(0, str(Path(__file__).parent))

from src.parsers.fastq_parser import parse_fastq, get_fastq_stats, print_stats


def create_demo_fastq(path: Path) -> Path:
    """Create a small demo FASTQ file with 5 reads of mixed quality."""
    reads = [
    ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", "I" * 52),
    ("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG", "I" * 52),
    ("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "!" * 52),
    ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I" * 53),
    ("ATGC", "IIII"),
    ]
    with gzip.open(path, "wt") as f:
        for i, (seq, qual) in enumerate(reads, 1):
            f.write(f"@demo_read_{i} length={len(seq)}\n")
            f.write(f"{seq}\n")
            f.write(f"+\n")
            f.write(f"{qual}\n")

    return path


def main():
    print("\n" + "="*55)
    print("  MODULE 2 DEMO — FASTQ Parser")
    print("="*55)

    # ── 1. Create demo FASTQ ──────────────────────────────────
    demo_file = Path("demo_reads.fastq.gz")
    create_demo_fastq(demo_file)
    print(f"\n✅ Created demo FASTQ: {demo_file}")

    # ── 2. Parse reads one by one ────────────────────────────
    print("\n📖 Parsing reads individually:\n")
    for record in parse_fastq(demo_file):
        status = "✅ PASS" if record.is_good_quality() else "❌ FAIL"
        print(f"  {status} | Read: {record.header:<20} | "
              f"Length: {record.length:>3}bp | "
              f"GC: {record.gc_content:>5}% | "
              f"MeanQ: {record.mean_quality:>5}")

    # ── 3. Get overall statistics ─────────────────────────────
    print()
    stats = get_fastq_stats(demo_file, min_quality=20, min_length=50)
    print_stats(stats)

    # ── 4. Clean up demo file ─────────────────────────────────
    demo_file.unlink()
    print("✅ Demo complete! Module 2 parsers are working correctly.\n")

    # ── 5. OPTIONAL: Fetch real GEO data ─────────────────────
    # Uncomment the lines below to test with real GEO data.
    # This downloads metadata from the mammary gland study
    # you've been working with in your R tutorial!
    #
    # from src.parsers.geo_fetcher import GEOFetcher
    # fetcher = GEOFetcher("GSE60450")
    # metadata = fetcher.download_metadata()
    # fetcher.print_summary()


if __name__ == "__main__":
    main()
