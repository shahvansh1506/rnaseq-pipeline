"""
aligner.py — Aligns trimmed RNA-seq reads to a reference genome.

What is alignment?
- After trimming, we have millions of short reads (75-150bp each)
- We need to find WHERE in the genome each read came from
- HISAT2 is a fast, splice-aware aligner for RNA-seq data
- "Splice-aware" means it understands that RNA reads can span
  across exon-exon junctions (introns are spliced out in mRNA)

What this file does:
- Downloads reference genome index (or builds one from FASTA)
- Runs HISAT2 to align reads -> produces SAM/BAM files
- Converts and sorts BAM files using samtools
- Calculates alignment statistics

Output:
- BAM file: Binary Alignment Map — stores where each read aligned
- BAI file: BAM index — allows fast random access to BAM file

Usage:
    aligner = Aligner()
    result  = aligner.align_single_end(
        fastq_file = Path("data/processed/sample_trimmed.fastq.gz"),
        index_path = Path("data/reference/genome"),
    )
"""

import subprocess
import logging
import re
from pathlib import Path
from dataclasses import dataclass

from src.config import (
    HISAT2_PATH, HISAT2_BUILD_PATH, SAMTOOLS_PATH,
    ALIGNED_DIR, REF_DIR, THREADS, LOG_FORMAT, LOG_LEVEL
)

# ── Setup logger ──────────────────────────────────────────────
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════
# DATA CLASS — alignment results
# ══════════════════════════════════════════════════════════════

@dataclass
class AlignResult:
    """
    Stores alignment statistics for one sample.

    Good RNA-seq alignment rates:
    - Overall alignment > 80% is acceptable
    - Overall alignment > 90% is good
    - Overall alignment > 95% is excellent
    """
    sample_name         : str
    total_reads         : int   = 0
    aligned_once        : int   = 0
    aligned_multiple    : int   = 0
    unaligned           : int   = 0
    overall_align_rate  : float = 0.0
    bam_file            : str   = ""

    def summary(self) -> str:
        status = "GOOD" if self.overall_align_rate >= 80 else "LOW"
        return (f"[{status}] {self.sample_name:<30} | "
                f"Total: {self.total_reads:>10,} | "
                f"Aligned: {self.overall_align_rate:.1f}%")


# ══════════════════════════════════════════════════════════════
# MAIN CLASS
# ══════════════════════════════════════════════════════════════

class Aligner:
    """
    Aligns RNA-seq reads to a reference genome using HISAT2.

    HISAT2 workflow:
    1. Build genome index (one time only)
    2. Align reads -> SAM file
    3. Convert SAM -> BAM (binary, smaller)
    4. Sort BAM by coordinate
    5. Index BAM for fast access

    Usage:
        aligner = Aligner()
        result  = aligner.align_single_end(fastq, index_path)
    """

    def __init__(self, output_dir: Path = None,
                 threads: int = THREADS):
        """
        Args:
            output_dir: Where to save BAM files
            threads   : CPU threads for alignment
        """
        self.output_dir = output_dir or ALIGNED_DIR
        self.threads    = threads
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ──────────────────────────────────────────────────────────
    # BUILD GENOME INDEX
    # ──────────────────────────────────────────────────────────

    def build_index(self, fasta_file: Path,
                    index_name: str = "genome") -> Path:
        """
        Build a HISAT2 genome index from a FASTA file.

        This only needs to be done ONCE per genome.
        The index allows HISAT2 to rapidly search for
        where each read comes from.

        For mouse (mm10) or human (hg38), you can also
        download pre-built indexes from:
        https://daehwankimlab.github.io/hisat2/download/

        Args:
            fasta_file: Path to genome FASTA file
            index_name: Base name for index files

        Returns:
            Path to index directory
        """
        fasta_file = Path(fasta_file)
        index_dir  = REF_DIR / "hisat2_index"
        index_dir.mkdir(parents=True, exist_ok=True)
        index_path = index_dir / index_name

        if not fasta_file.exists():
            raise FileNotFoundError(f"FASTA not found: {fasta_file}")

        logger.info(f"Building HISAT2 index from: {fasta_file.name}")
        logger.info("This may take 30-60 minutes for a full genome...")

        cmd = [
            HISAT2_BUILD_PATH,
            "-p", str(self.threads),
            str(fasta_file),
            str(index_path),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"Index build failed:\n{result.stderr}")
            raise RuntimeError("HISAT2 index build failed")

        logger.info(f"Index built successfully: {index_path}")
        return index_path

    # ──────────────────────────────────────────────────────────
    # ALIGNMENT
    # ──────────────────────────────────────────────────────────

    def align_single_end(self, fastq_file : Path,
                         index_path : Path) -> AlignResult:
        """
        Align single-end reads to the genome.

        Steps:
        1. Run HISAT2 -> SAM file
        2. Convert SAM to BAM (samtools view)
        3. Sort BAM by coordinate (samtools sort)
        4. Index BAM (samtools index)
        5. Delete SAM file (save disk space)

        Args:
            fastq_file: Trimmed FASTQ file
            index_path: Path to HISAT2 genome index

        Returns:
            AlignResult with statistics
        """
        fastq_file  = Path(fastq_file)
        sample_name = fastq_file.stem.replace(".fastq", "").replace(".gz", "")

        # Output file paths
        sam_file = self.output_dir / f"{sample_name}.sam"
        bam_file = self.output_dir / f"{sample_name}.bam"
        sorted_bam = self.output_dir / f"{sample_name}_sorted.bam"

        logger.info(f"Aligning: {fastq_file.name}")

        # ── Step 1: Run HISAT2 ────────────────────────────────
        cmd_hisat2 = [
            HISAT2_PATH,
            "-x", str(index_path),     # Genome index
            "-U", str(fastq_file),     # Input reads (U = unpaired)
            "-S", str(sam_file),       # Output SAM file
            "-p", str(self.threads),   # Threads
            "--dta",                   # Downstream transcriptome assembly
            "--no-unal",               # Don't output unaligned reads
        ]

        result = subprocess.run(
            cmd_hisat2, capture_output=True, text=True
        )

        # Parse alignment stats from HISAT2 stderr
        align_result = self._parse_hisat2_stats(
            result.stderr, sample_name, str(sorted_bam)
        )

        if result.returncode != 0:
            logger.error(f"HISAT2 failed:\n{result.stderr}")
            raise RuntimeError(f"Alignment failed for {sample_name}")

        # ── Step 2: SAM -> BAM ────────────────────────────────
        self._sam_to_bam(sam_file, bam_file)

        # ── Step 3: Sort BAM ──────────────────────────────────
        self._sort_bam(bam_file, sorted_bam)

        # ── Step 4: Index BAM ─────────────────────────────────
        self._index_bam(sorted_bam)

        # ── Step 5: Clean up SAM and unsorted BAM ────────────
        sam_file.unlink(missing_ok=True)
        bam_file.unlink(missing_ok=True)

        logger.info(f"Alignment complete: {align_result.summary()}")
        return align_result

    def align_paired_end(self, fastq_r1  : Path,
                         fastq_r2  : Path,
                         index_path: Path) -> AlignResult:
        """
        Align paired-end reads to the genome.

        Paired-end = two FASTQ files (R1 and R2) per sample.
        HISAT2 uses both files together for better accuracy.

        Args:
            fastq_r1  : Read 1 trimmed FASTQ file
            fastq_r2  : Read 2 trimmed FASTQ file
            index_path: Path to HISAT2 genome index

        Returns:
            AlignResult with statistics
        """
        fastq_r1    = Path(fastq_r1)
        fastq_r2    = Path(fastq_r2)
        sample_name = fastq_r1.stem.replace("_R1", "").replace(".fastq", "")

        sam_file   = self.output_dir / f"{sample_name}.sam"
        bam_file   = self.output_dir / f"{sample_name}.bam"
        sorted_bam = self.output_dir / f"{sample_name}_sorted.bam"

        logger.info(f"Aligning paired-end: {fastq_r1.name} + {fastq_r2.name}")

        cmd_hisat2 = [
            HISAT2_PATH,
            "-x", str(index_path),
            "-1", str(fastq_r1),       # Read 1
            "-2", str(fastq_r2),       # Read 2
            "-S", str(sam_file),
            "-p", str(self.threads),
            "--dta",
            "--no-unal",
        ]

        result = subprocess.run(
            cmd_hisat2, capture_output=True, text=True
        )

        align_result = self._parse_hisat2_stats(
            result.stderr, sample_name, str(sorted_bam)
        )

        if result.returncode != 0:
            logger.error(f"HISAT2 failed:\n{result.stderr}")
            raise RuntimeError(f"Paired alignment failed for {sample_name}")

        self._sam_to_bam(sam_file, bam_file)
        self._sort_bam(bam_file, sorted_bam)
        self._index_bam(sorted_bam)

        sam_file.unlink(missing_ok=True)
        bam_file.unlink(missing_ok=True)

        logger.info(f"Paired alignment complete: {align_result.summary()}")
        return align_result

    def align_batch(self, fastq_dir  : Path,
                    index_path : Path,
                    paired_end : bool = False) -> list:
        """
        Align all FASTQ files in a directory.

        Args:
            fastq_dir : Directory with trimmed FASTQ files
            index_path: HISAT2 genome index path
            paired_end: True if reads are paired-end

        Returns:
            List of AlignResult objects
        """
        fastq_dir = Path(fastq_dir)
        results   = []

        if paired_end:
            # Find R1 files and pair with R2
            r1_files = sorted(fastq_dir.glob("*_R1*.fastq.gz"))
            for r1 in r1_files:
                r2 = Path(str(r1).replace("_R1", "_R2"))
                if r2.exists():
                    try:
                        result = self.align_paired_end(r1, r2, index_path)
                        results.append(result)
                    except Exception as e:
                        logger.error(f"Failed: {r1.name}: {e}")
        else:
            # Single-end: align each file independently
            fastq_files = list(fastq_dir.glob("*.fastq.gz"))
            for fastq in fastq_files:
                try:
                    result = self.align_single_end(fastq, index_path)
                    results.append(result)
                except Exception as e:
                    logger.error(f"Failed: {fastq.name}: {e}")

        return results

    # ──────────────────────────────────────────────────────────
    # SAMTOOLS HELPERS
    # ──────────────────────────────────────────────────────────

    def _sam_to_bam(self, sam_file: Path, bam_file: Path) -> None:
        """
        Convert SAM to BAM format.

        SAM = text format, human readable but large
        BAM = binary format, compressed, much smaller
        """
        logger.info(f"Converting SAM to BAM: {bam_file.name}")
        cmd = [
            SAMTOOLS_PATH, "view",
            "-bS",              # b=output BAM, S=input SAM
            "-o", str(bam_file),
            str(sam_file),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"SAM to BAM failed: {result.stderr}")

    def _sort_bam(self, bam_file: Path, sorted_bam: Path) -> None:
        """
        Sort BAM by genomic coordinate.
        Required for indexing and for featureCounts.
        """
        logger.info(f"Sorting BAM: {sorted_bam.name}")
        cmd = [
            SAMTOOLS_PATH, "sort",
            "-o", str(sorted_bam),
            "-@", str(self.threads),
            str(bam_file),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"BAM sort failed: {result.stderr}")

    def _index_bam(self, bam_file: Path) -> None:
        """
        Create BAM index (.bai file).
        Allows fast random access to any region of the genome.
        """
        logger.info(f"Indexing BAM: {bam_file.name}")
        cmd = [SAMTOOLS_PATH, "index", str(bam_file)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"BAM index failed: {result.stderr}")

    # ──────────────────────────────────────────────────────────
    # PARSING
    # ──────────────────────────────────────────────────────────

    def _parse_hisat2_stats(self, stderr     : str,
                             sample_name : str,
                             bam_file    : str) -> AlignResult:
        """
        Parse HISAT2 alignment statistics from stderr output.

        HISAT2 prints something like:
            1000000 reads; of these:
              1000000 (100.00%) were unpaired; of these:
                50000 (5.00%) aligned 0 times
                900000 (90.00%) aligned exactly 1 time
                50000 (5.00%) aligned >1 times
            95.00% overall alignment rate
        """
        result = AlignResult(
            sample_name = sample_name,
            bam_file    = bam_file,
        )

        # Total reads
        match = re.search(r"(\d+) reads; of these", stderr)
        if match:
            result.total_reads = int(match.group(1))

        # Aligned exactly once
        match = re.search(r"(\d+).*aligned exactly 1 time", stderr)
        if match:
            result.aligned_once = int(match.group(1))

        # Aligned multiple times
        match = re.search(r"(\d+).*aligned >1 times", stderr)
        if match:
            result.aligned_multiple = int(match.group(1))

        # Unaligned
        match = re.search(r"(\d+).*aligned 0 times", stderr)
        if match:
            result.unaligned = int(match.group(1))

        # Overall alignment rate
        match = re.search(r"([0-9.]+)% overall alignment rate", stderr)
        if match:
            result.overall_align_rate = float(match.group(1))

        return result

    def print_summary(self, results: list) -> None:
        """Print alignment summary for all samples."""
        if not results:
            print("No alignment results.")
            return

        good = sum(1 for r in results if r.overall_align_rate >= 80)

        print(f"\n{'='*65}")
        print(f"  ALIGNMENT SUMMARY  |  {len(results)} samples  |  "
              f"Good (>80%): {good}")
        print(f"{'='*65}")
        for r in results:
            print(f"  {r.summary()}")
        print(f"{'='*65}\n")
