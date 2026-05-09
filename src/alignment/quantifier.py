"""
quantifier.py — Counts reads per gene to build a count matrix.

What is quantification?
- After alignment, we know WHERE each read mapped in the genome
- Now we need to count HOW MANY reads mapped to each gene
- This gives us a count matrix:
    rows    = genes
    columns = samples
    values  = number of reads mapping to that gene
- This count matrix is the INPUT to DESeq2

What this file does:
- Runs featureCounts on BAM files to count reads per gene
- Parses the output into a clean pandas DataFrame
- Saves the count matrix as a CSV file

Usage:
    quantifier = Quantifier()
    count_matrix = quantifier.count_reads(
        bam_files = [Path("data/aligned/sample1_sorted.bam")],
        gtf_file  = Path("data/reference/genome.gtf"),
    )
"""

import subprocess
import logging
import pandas as pd
from pathlib import Path

from src.config import (
    FEATURECOUNTS_PATH, COUNTS_DIR,
    REF_DIR, THREADS, LOG_FORMAT, LOG_LEVEL
)

# ── Setup logger ──────────────────────────────────────────────
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


class Quantifier:
    """
    Counts reads per gene using featureCounts.

    featureCounts (from the Subread package) is the fastest
    and most widely used read counting tool for RNA-seq.

    It uses a GTF annotation file to know where each gene is,
    then counts how many reads overlap each gene's exons.

    Usage:
        q = Quantifier()
        matrix = q.count_reads(bam_files, gtf_file)
        q.save_count_matrix(matrix, "data/counts/counts.csv")
    """

    def __init__(self, output_dir: Path = None,
                 threads: int = THREADS):
        """
        Args:
            output_dir: Where to save count files
            threads   : CPU threads for counting
        """
        self.output_dir = output_dir or COUNTS_DIR
        self.threads    = threads
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ──────────────────────────────────────────────────────────
    # MAIN COUNTING METHOD
    # ──────────────────────────────────────────────────────────

    def count_reads(self, bam_files  : list,
                    gtf_file   : Path,
                    paired_end : bool = False,
                    strand     : int  = 0) -> pd.DataFrame:
        """
        Count reads per gene across all samples.

        Args:
            bam_files : List of sorted BAM file paths
            gtf_file  : GTF annotation file (gene coordinates)
            paired_end: True if reads are paired-end
            strand    : Strandedness
                        0 = unstranded (most common)
                        1 = stranded
                        2 = reverse stranded

        Returns:
            DataFrame with genes as rows, samples as columns
        """
        gtf_file  = Path(gtf_file)
        bam_files = [Path(b) for b in bam_files]

        # Validate inputs
        if not gtf_file.exists():
            raise FileNotFoundError(f"GTF file not found: {gtf_file}")

        missing = [b for b in bam_files if not b.exists()]
        if missing:
            raise FileNotFoundError(
                f"BAM files not found: {[b.name for b in missing]}"
            )

        output_file = self.output_dir / "raw_counts.txt"

        logger.info(f"Counting reads in {len(bam_files)} BAM files...")

        # Build featureCounts command
        cmd = [
            FEATURECOUNTS_PATH,
            "-T", str(self.threads),   # Threads
            "-a", str(gtf_file),       # Annotation GTF
            "-o", str(output_file),    # Output file
            "-s", str(strand),         # Strandedness
            "-t", "exon",              # Feature type to count
            "-g", "gene_id",           # Attribute to group by
        ]

        # Add paired-end flag if needed
        if paired_end:
            cmd.append("-p")           # Paired end mode
            cmd.append("-B")           # Both reads must map

        # Add all BAM files at the end
        cmd += [str(b) for b in bam_files]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"featureCounts failed:\n{result.stderr}")
            raise RuntimeError("featureCounts failed")

        logger.info("featureCounts complete. Parsing output...")

        # Parse the output file into a clean DataFrame
        count_matrix = self._parse_featurecounts(output_file, bam_files)

        logger.info(f"Count matrix: {count_matrix.shape[0]:,} genes x "
                    f"{count_matrix.shape[1]} samples")

        return count_matrix

    # ──────────────────────────────────────────────────────────
    # PARSING
    # ──────────────────────────────────────────────────────────

    def _parse_featurecounts(self, output_file: Path,
                              bam_files: list) -> pd.DataFrame:
        """
        Parse featureCounts output into a clean count matrix.

        featureCounts output format:
        - First line: comment starting with #
        - Second line: headers
        - Remaining lines: gene_id, chr, start, end, strand, length, counts...

        We keep only the gene_id and count columns.

        Args:
            output_file: Path to featureCounts output file
            bam_files  : List of BAM files (used to name columns)

        Returns:
            Clean DataFrame: genes x samples
        """
        # Read the raw output (skip the first comment line)
        df = pd.read_csv(output_file, sep="\t", comment="#")

        # The first column is gene ID
        # Columns 2-6 are annotation info (chr, start, end, strand, length)
        # Columns 7+ are the actual counts for each sample

        # Extract gene IDs
        gene_ids = df.iloc[:, 0]

        # Extract count columns (everything after column 6)
        count_cols = df.iloc[:, 6:]

        # Rename columns to clean sample names
        sample_names = [Path(b).stem.replace("_sorted", "")
                        for b in bam_files]
        count_cols.columns = sample_names

        # Build clean DataFrame
        count_matrix = count_cols.copy()
        count_matrix.index = gene_ids
        count_matrix.index.name = "gene_id"

        return count_matrix

    # ──────────────────────────────────────────────────────────
    # SAVING & LOADING
    # ──────────────────────────────────────────────────────────

    def save_count_matrix(self, count_matrix : pd.DataFrame,
                          output_path  : Path = None) -> Path:
        """
        Save count matrix to CSV file.

        Args:
            count_matrix: DataFrame with gene counts
            output_path : Where to save. Auto-named if None.

        Returns:
            Path to saved file
        """
        output_path = output_path or self.output_dir / "count_matrix.csv"
        count_matrix.to_csv(output_path)
        logger.info(f"Count matrix saved: {output_path}")
        return output_path

    def load_count_matrix(self, filepath: Path) -> pd.DataFrame:
        """
        Load a saved count matrix from CSV.

        Args:
            filepath: Path to count matrix CSV

        Returns:
            DataFrame with gene counts
        """
        df = pd.read_csv(filepath, index_col=0)
        logger.info(f"Loaded count matrix: {df.shape[0]:,} genes x "
                    f"{df.shape[1]} samples")
        return df

    def print_summary(self, count_matrix: pd.DataFrame) -> None:
        """Print a summary of the count matrix."""
        print(f"\n{'='*55}")
        print(f"  COUNT MATRIX SUMMARY")
        print(f"{'='*55}")
        print(f"  Genes (rows)   : {count_matrix.shape[0]:>10,}")
        print(f"  Samples (cols) : {count_matrix.shape[1]:>10,}")
        print(f"  Total counts   : {count_matrix.values.sum():>10,.0f}")
        print(f"  Zero-count genes: {(count_matrix.sum(axis=1) == 0).sum():>9,}")
        print(f"\n  Samples:")
        for col in count_matrix.columns:
            total = count_matrix[col].sum()
            print(f"    {col:<35} {total:>12,.0f} reads")
        print(f"{'='*55}\n")
