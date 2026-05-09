"""
fastq_parser.py — Memory-efficient FASTQ file parser.

What is a FASTQ file?
- Standard format for storing raw DNA/RNA sequencing reads
- Each read = 4 lines:
    Line 1: @ + read ID      e.g. @SRR1552444.1
    Line 2: DNA sequence     e.g. ATCGATCGATCG...
    Line 3: + (separator)    e.g. +
    Line 4: Quality scores   e.g. IIIIIIIIIIII...

What this file does:
- Parses FASTQ files one read at a time (generator pattern)
- Handles both plain (.fastq) and compressed (.fastq.gz) files
- Computes QC statistics (GC content, quality scores, read lengths)
- Never loads the whole file into memory — safe for files of any size

Usage:
    # Iterate through reads
    for record in parse_fastq("sample.fastq.gz"):
        print(record.header, record.gc_content)

    # Get summary statistics
    stats = get_fastq_stats("sample.fastq.gz")
    print(stats)
"""

import gzip
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Generator
from collections import defaultdict

from src.config import LOG_FORMAT, LOG_LEVEL

# ── Setup logger ──────────────────────────────────────────────
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════
# DATA CLASS — represents one sequencing read
# ══════════════════════════════════════════════════════════════

@dataclass
class FastqRecord:
    """
    Holds all information for a single sequencing read.

    A dataclass automatically generates __init__, __repr__ etc.
    Think of it as a neat container for one FASTQ read's 4 lines.
    """
    header  : str    # Read ID (without the leading @)
    sequence: str    # DNA/RNA sequence
    plus    : str    # Always just "+" — separator line
    quality : str    # Phred quality score string

    # ── Computed properties (calculated on demand) ─────────

    @property
    def length(self) -> int:
        """Number of bases in this read."""
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        """
        Percentage of bases that are G or C.
        Typical RNA-seq reads: 45-60% GC content.
        Very high (>70%) or low (<30%) may indicate problems.
        """
        if self.length == 0:
            return 0.0
        gc = self.sequence.count("G") + self.sequence.count("C")
        return round(gc / self.length * 100, 2)

    @property
    def mean_quality(self) -> float:
        """
        Average Phred quality score across all bases.

        Phred score conversion:
        - ASCII character → subtract 33 → Phred score
        - Phred 20 = 99%   base call accuracy
        - Phred 30 = 99.9% base call accuracy  ← good quality
        - Phred 40 = 99.99% base call accuracy ← excellent

        We want mean quality > 20 for good RNA-seq data.
        """
        if not self.quality:
            return 0.0
        scores = [ord(char) - 33 for char in self.quality]
        return round(sum(scores) / len(scores), 2)

    @property
    def base_composition(self) -> dict:
        """Count of each nucleotide (A, T, G, C, N) in the read."""
        counts = defaultdict(int)
        for base in self.sequence:
            counts[base] += 1
        return dict(counts)

    def is_good_quality(self, min_quality: float = 20.0,
                        min_length: int = 50) -> bool:
        """
        Quick check: is this read good enough to use?

        Args:
            min_quality: Minimum mean Phred score (default 20)
            min_length : Minimum read length in bp (default 50)

        Returns:
            True if read passes both filters
        """
        return self.mean_quality >= min_quality and self.length >= min_length


# ══════════════════════════════════════════════════════════════
# GENERATOR FUNCTION — the core parser
# ══════════════════════════════════════════════════════════════

def parse_fastq(filepath: Path) -> Generator[FastqRecord, None, None]:
    """
    Generator that yields one FastqRecord at a time.

    WHY a generator?
    - A normal function returns ALL results at once → loads everything into RAM
    - A generator yields ONE result at a time → uses constant, tiny RAM
    - Perfect for FASTQ files that can be 50GB+

    Handles both:
    - Plain FASTQ files (.fastq)
    - Gzip compressed FASTQ files (.fastq.gz) ← most common in practice

    Args:
        filepath: Path to FASTQ or FASTQ.gz file

    Yields:
        FastqRecord for each read in the file

    Example:
        for read in parse_fastq("data/raw/sample.fastq.gz"):
            if read.mean_quality > 30:
                print(f"High quality read: {read.header}")
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"FASTQ file not found: {filepath}")

    # Choose opener based on file extension
    is_gzipped = str(filepath).endswith(".gz")
    opener     = gzip.open if is_gzipped else open

    logger.info(f"Parsing {'gzipped ' if is_gzipped else ''}FASTQ: {filepath.name}")

    with opener(filepath, "rt") as f:
        while True:
            # Read exactly 4 lines = 1 FASTQ record
            header   = f.readline().strip()
            sequence = f.readline().strip()
            plus     = f.readline().strip()
            quality  = f.readline().strip()

            # Empty header = end of file
            if not header:
                break

            # Validate FASTQ format
            if not header.startswith("@"):
                logger.warning(f"Skipping malformed record — header: {header[:50]}")
                continue

            if len(sequence) != len(quality):
                logger.warning(f"Skipping read {header[:30]}: "
                               f"sequence/quality length mismatch")
                continue

            # yield pauses here, returns the record, resumes on next iteration
            yield FastqRecord(
                header   = header[1:],   # Strip leading "@"
                sequence = sequence,
                plus     = plus,
                quality  = quality,
            )


# ══════════════════════════════════════════════════════════════
# STATISTICS FUNCTION
# ══════════════════════════════════════════════════════════════

def get_fastq_stats(filepath: Path,
                    min_quality: float = 20.0,
                    min_length : int   = 50) -> dict:
    """
    Compute QC statistics for an entire FASTQ file.

    Streams the file using parse_fastq() — RAM usage stays
    constant no matter how large the file is.

    Args:
        filepath   : Path to FASTQ or FASTQ.gz file
        min_quality: Phred score cutoff for "good quality"
        min_length : Minimum read length cutoff

    Returns:
        dict with these keys:
        - total_reads      : Total number of reads
        - total_bases      : Total bases sequenced
        - mean_read_length : Average read length
        - gc_content_pct   : GC% across all reads
        - mean_quality     : Average Phred quality score
        - low_quality_reads: Reads below min_quality
        - low_quality_pct  : Percentage of low quality reads
        - short_reads      : Reads below min_length
        - pass_filter_pct  : % of reads passing both filters
    """
    # Counters — all start at zero
    total_reads       = 0
    total_bases       = 0
    total_gc          = 0
    quality_sum       = 0.0
    low_quality_reads = 0
    short_reads       = 0
    pass_filter       = 0

    logger.info(f"Computing statistics for: {Path(filepath).name}")

    # Stream through file — one read at a time
    for record in parse_fastq(filepath):
        total_reads += 1
        total_bases += record.length
        total_gc    += record.sequence.count("G") + record.sequence.count("C")
        quality_sum += record.mean_quality

        if record.mean_quality < min_quality:
            low_quality_reads += 1
        if record.length < min_length:
            short_reads += 1
        if record.is_good_quality(min_quality, min_length):
            pass_filter += 1

    # Avoid division by zero if file is empty
    if total_reads == 0:
        logger.warning("No reads found in file!")
        return {}

    stats = {
        "file"             : Path(filepath).name,
        "total_reads"      : total_reads,
        "total_bases"      : total_bases,
        "mean_read_length" : round(total_bases / total_reads, 1),
        "gc_content_pct"   : round(total_gc / total_bases * 100, 2),
        "mean_quality"     : round(quality_sum / total_reads, 2),
        "low_quality_reads": low_quality_reads,
        "low_quality_pct"  : round(low_quality_reads / total_reads * 100, 2),
        "short_reads"      : short_reads,
        "pass_filter"      : pass_filter,
        "pass_filter_pct"  : round(pass_filter / total_reads * 100, 2),
    }

    logger.info(f"Stats complete: {total_reads:,} reads, "
                f"GC={stats['gc_content_pct']}%, "
                f"MeanQ={stats['mean_quality']}")

    return stats


def print_stats(stats: dict) -> None:
    """Pretty-print FASTQ statistics to the console."""
    if not stats:
        print("No statistics to display.")
        return

    print(f"\n{'='*45}")
    print(f"  FASTQ QC Report: {stats.get('file', 'unknown')}")
    print(f"{'='*45}")
    print(f"  Total reads      : {stats['total_reads']:>12,}")
    print(f"  Total bases      : {stats['total_bases']:>12,}")
    print(f"  Mean read length : {stats['mean_read_length']:>11} bp")
    print(f"  GC content       : {stats['gc_content_pct']:>11} %")
    print(f"  Mean quality     : {stats['mean_quality']:>11} (Phred)")
    print(f"  Low quality reads: {stats['low_quality_reads']:>12,} "
          f"({stats['low_quality_pct']}%)")
    print(f"  Short reads      : {stats['short_reads']:>12,}")
    print(f"  Pass filter      : {stats['pass_filter']:>12,} "
          f"({stats['pass_filter_pct']}%)")
    print(f"{'='*45}\n")
