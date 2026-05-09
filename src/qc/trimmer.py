"""
trimmer.py — Trims low-quality bases and adapter sequences from reads.

What is trimming?
- Raw reads often have low quality bases at the 3' end
- Adapter sequences (lab reagents) can get sequenced by mistake
- Reads that are too short after trimming are discarded
- We remove all of these before alignment to improve accuracy

What this file does:
- Wraps Trimmomatic (industry standard trimming tool)
- Supports both single-end and paired-end reads
- Falls back to a pure Python trimmer if Trimmomatic is not available

Usage:
    trimmer = Trimmer()
    result  = trimmer.trim_single_end(
        input_file  = Path("data/raw/sample.fastq.gz"),
        output_file = Path("data/processed/sample_trimmed.fastq.gz"),
    )
    print(result.summary())
"""

import subprocess
import logging
import re
import gzip
from pathlib import Path
from dataclasses import dataclass

from src.config import (
    TRIMMOMATIC_PATH, PROCESSED_DIR,
    MIN_QUALITY, MIN_LENGTH, LOG_FORMAT, LOG_LEVEL
)
from src.parsers.fastq_parser import parse_fastq, FastqRecord

# ── Setup logger ──────────────────────────────────────────────
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════
# DATA CLASS — trimming results
# ══════════════════════════════════════════════════════════════

@dataclass
class TrimResult:
    """Stores statistics from a trimming run."""
    sample_name     : str
    input_reads     : int   = 0
    surviving_reads : int   = 0
    dropped_reads   : int   = 0
    survival_pct    : float = 0.0
    output_file     : str   = ""

    def summary(self) -> str:
        return (f"{self.sample_name}: "
                f"{self.surviving_reads:,} / {self.input_reads:,} "
                f"reads kept ({self.survival_pct:.1f}%)")


# ══════════════════════════════════════════════════════════════
# MAIN CLASS
# ══════════════════════════════════════════════════════════════

class Trimmer:
    """
    Trims reads using Trimmomatic or a built-in Python fallback.

    Trimmomatic steps applied (in order):
    1. LEADING       - remove low quality bases from start
    2. TRAILING      - remove low quality bases from end
    3. SLIDINGWINDOW - cut when average quality drops in a window
    4. MINLEN        - discard reads shorter than minimum length
    """

    # Common Illumina adapter sequences
    ADAPTER_SEQUENCES = [
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",   # TruSeq Read 1
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",   # TruSeq Read 2
        "CTGTCTCTTATACACATCT",                  # Nextera
    ]

    def __init__(self,
                 output_dir  : Path = None,
                 min_quality : int  = MIN_QUALITY,
                 min_length  : int  = MIN_LENGTH):
        """
        Args:
            output_dir : Where to save trimmed files
            min_quality: Minimum Phred quality score to keep
            min_length : Minimum read length to keep
        """
        self.output_dir  = output_dir or PROCESSED_DIR
        self.min_quality = min_quality
        self.min_length  = min_length
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ──────────────────────────────────────────────────────────
    # TRIMMOMATIC (external tool)
    # ──────────────────────────────────────────────────────────

    def trim_single_end(self, input_file : Path,
                        output_file: Path = None) -> TrimResult:
        """
        Trim a single-end FASTQ file using Trimmomatic.

        Single-end = one FASTQ file per sample.

        Args:
            input_file : Path to raw FASTQ.gz file
            output_file: Path for trimmed output. Auto-named if None.

        Returns:
            TrimResult with statistics
        """
        input_file  = Path(input_file)

        if output_file is None:
            stem        = input_file.stem.replace(".fastq", "")
            output_file = self.output_dir / f"{stem}_trimmed.fastq.gz"

        sample_name = input_file.stem.replace(".fastq", "").replace(".gz", "")
        logger.info(f"Trimming: {input_file.name}")

        cmd = [
            TRIMMOMATIC_PATH, "SE",
            "-phred33",
            str(input_file),
            str(output_file),
            f"LEADING:{self.min_quality}",
            f"TRAILING:{self.min_quality}",
            f"SLIDINGWINDOW:4:{self.min_quality}",
            f"MINLEN:{self.min_length}",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        trim_result = self._parse_trimmomatic_output(
            result.stderr, sample_name, str(output_file)
        )

        if result.returncode != 0 and trim_result.input_reads == 0:
            logger.warning("Trimmomatic not available. "
                           "Using Python fallback trimmer...")
            return self.trim_python_fallback(input_file, output_file)

        logger.info(f"Trimming complete: {trim_result.summary()}")
        return trim_result

    def trim_paired_end(self, input_r1  : Path,
                        input_r2  : Path,
                        output_r1 : Path = None,
                        output_r2 : Path = None) -> TrimResult:
        """
        Trim paired-end reads using Trimmomatic.

        Paired-end = two FASTQ files per sample (R1 and R2).
        Both files must be trimmed together to keep pairs in sync.

        Args:
            input_r1 : Read 1 FASTQ file
            input_r2 : Read 2 FASTQ file
            output_r1: Trimmed Read 1 output
            output_r2: Trimmed Read 2 output

        Returns:
            TrimResult with statistics
        """
        input_r1 = Path(input_r1)
        input_r2 = Path(input_r2)

        if output_r1 is None:
            stem      = input_r1.stem.replace(".fastq", "")
            output_r1 = self.output_dir / f"{stem}_trimmed_R1.fastq.gz"
        if output_r2 is None:
            stem      = input_r2.stem.replace(".fastq", "")
            output_r2 = self.output_dir / f"{stem}_trimmed_R2.fastq.gz"

        unpaired_r1 = self.output_dir / "unpaired_R1.fastq.gz"
        unpaired_r2 = self.output_dir / "unpaired_R2.fastq.gz"
        sample_name = input_r1.stem.replace(".fastq", "").replace(".gz", "")

        logger.info(f"Trimming paired-end: {input_r1.name} + {input_r2.name}")

        cmd = [
            TRIMMOMATIC_PATH, "PE", "-phred33",
            str(input_r1), str(input_r2),
            str(output_r1), str(unpaired_r1),
            str(output_r2), str(unpaired_r2),
            f"LEADING:{self.min_quality}",
            f"TRAILING:{self.min_quality}",
            f"SLIDINGWINDOW:4:{self.min_quality}",
            f"MINLEN:{self.min_length}",
        ]

        result      = subprocess.run(cmd, capture_output=True, text=True)
        trim_result = self._parse_trimmomatic_output(
            result.stderr, sample_name, str(output_r1)
        )

        logger.info(f"Paired trimming complete: {trim_result.summary()}")
        return trim_result

    # ──────────────────────────────────────────────────────────
    # PYTHON FALLBACK TRIMMER
    # ──────────────────────────────────────────────────────────

    def trim_python_fallback(self, input_file : Path,
                             output_file: Path) -> TrimResult:
        """
        Pure Python trimmer used when Trimmomatic is not installed.

        Applies:
        1. Quality trimming from 3' end
        2. Adapter trimming
        3. Minimum length filter

        Args:
            input_file : Raw FASTQ file
            output_file: Trimmed output file

        Returns:
            TrimResult with statistics
        """
        input_file  = Path(input_file)
        output_file = Path(output_file)
        sample_name = input_file.stem.replace(".fastq", "").replace(".gz", "")

        logger.info(f"Python fallback trimmer: {input_file.name}")

        input_reads     = 0
        surviving_reads = 0

        with gzip.open(output_file, "wt") as out_f:
            for record in parse_fastq(input_file):
                input_reads += 1

                # Step 1: Trim low quality bases from 3' end
                trimmed = self._quality_trim_3prime(record)

                # Step 2: Remove adapter sequences
                trimmed = self._trim_adapters(trimmed)

                # Step 3: Keep only reads long enough
                if trimmed.length >= self.min_length:
                    surviving_reads += 1
                    out_f.write(f"@{trimmed.header}\n")
                    out_f.write(f"{trimmed.sequence}\n")
                    out_f.write(f"+\n")
                    out_f.write(f"{trimmed.quality}\n")

        dropped  = input_reads - surviving_reads
        survival = (surviving_reads / input_reads * 100
                    if input_reads > 0 else 0.0)

        result = TrimResult(
            sample_name     = sample_name,
            input_reads     = input_reads,
            surviving_reads = surviving_reads,
            dropped_reads   = dropped,
            survival_pct    = round(survival, 2),
            output_file     = str(output_file),
        )

        logger.info(f"Python trimmer complete: {result.summary()}")
        return result

    # ──────────────────────────────────────────────────────────
    # PRIVATE HELPERS
    # ──────────────────────────────────────────────────────────

    def _quality_trim_3prime(self, record: FastqRecord) -> FastqRecord:
        """
        Trim low quality bases from the 3' (right) end of a read.

        Scans from right and removes bases below min_quality.
        """
        scores  = [ord(c) - 33 for c in record.quality]
        cut_pos = len(scores)

        for i in range(len(scores) - 1, -1, -1):
            if scores[i] >= self.min_quality:
                cut_pos = i + 1
                break
            cut_pos = i

        return FastqRecord(
            header   = record.header,
            sequence = record.sequence[:cut_pos],
            plus     = record.plus,
            quality  = record.quality[:cut_pos],
        )

    def _trim_adapters(self, record: FastqRecord) -> FastqRecord:
        """
        Remove adapter sequences from a read.

        Checks for any known adapter and trims from that point.
        """
        seq = record.sequence

        for adapter in self.ADAPTER_SEQUENCES:
            for start in range(len(seq) - 10):
                if adapter.startswith(seq[start:]):
                    seq = seq[:start]
                    return FastqRecord(
                        header   = record.header,
                        sequence = seq,
                        plus     = record.plus,
                        quality  = record.quality[:len(seq)],
                    )
        return record

    def _parse_trimmomatic_output(self, stderr     : str,
                                   sample_name : str,
                                   output_file : str) -> TrimResult:
        """
        Parse Trimmomatic statistics from its stderr output.

        Trimmomatic prints:
        'Input Reads: 1000000 Surviving: 950000 (95.00%) Dropped: 50000'
        """
        result = TrimResult(
            sample_name = sample_name,
            output_file = output_file,
        )

        match = re.search(
            r"Input Reads:\s+(\d+)\s+Surviving:\s+(\d+)\s+"
            r"\(([0-9.]+)%\)\s+Dropped:\s+(\d+)",
            stderr
        )

        if match:
            result.input_reads     = int(match.group(1))
            result.surviving_reads = int(match.group(2))
            result.survival_pct    = float(match.group(3))
            result.dropped_reads   = int(match.group(4))

        return result