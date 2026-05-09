"""
run_deseq2.py — Python wrapper that calls the DESeq2 R script.

What this file does:
- Takes count matrix and metadata as inputs
- Saves them as CSV files DESeq2 can read
- Calls de_analysis.R via subprocess
- Reads back and returns the results
- Handles errors clearly

Why a Python wrapper?
- Keeps our pipeline in Python
- Allows us to pass data programmatically
- Makes it easy to integrate with Streamlit dashboard
- Handles file I/O and error checking cleanly

Usage:
    runner = DESeq2Runner()
    results = runner.run(
        counts       = count_matrix_df,
        metadata     = metadata_df,
        condition_col= "condition",
    )
"""

import subprocess
import logging
import pandas as pd
from pathlib import Path

from src.config import (
    RSCRIPT_PATH, RESULTS_DIR, COUNTS_DIR,
    LOG_FORMAT, LOG_LEVEL
)

# ── Setup logger ──────────────────────────────────────────────
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


class DESeq2Runner:
    """
    Runs the DESeq2 R script from Python.

    Takes a count matrix DataFrame and metadata DataFrame,
    saves them to CSV, runs the R script, and returns results.

    Usage:
        runner  = DESeq2Runner()
        results = runner.run(counts, metadata, "condition")
        print(results.head())
    """

    def __init__(self, output_dir: Path = None):
        """
        Args:
            output_dir: Where to save results and plots.
                        Defaults to results/deseq2/
        """
        self.output_dir = output_dir or RESULTS_DIR / "deseq2"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Path to the R script
        self.r_script = Path(__file__).parent / "de_analysis.R"

    # ──────────────────────────────────────────────────────────
    # MAIN RUN METHOD
    # ──────────────────────────────────────────────────────────

    def run(self, counts       : pd.DataFrame,
            metadata     : pd.DataFrame,
            condition_col: str = "condition") -> pd.DataFrame:
        """
        Run full DESeq2 differential expression analysis.

        Steps:
        1. Save counts and metadata to CSV
        2. Run de_analysis.R via Rscript
        3. Load and return results

        Args:
            counts       : Count matrix (genes x samples)
            metadata     : Sample metadata (samples x conditions)
            condition_col: Column in metadata defining groups

        Returns:
            DataFrame with DE results for all comparisons
        """
        # ── Step 1: Save inputs to CSV ────────────────────────
        counts_path   = self.output_dir / "input_counts.csv"
        metadata_path = self.output_dir / "input_metadata.csv"

        logger.info("Saving inputs for DESeq2...")
        counts.to_csv(counts_path)
        metadata.to_csv(metadata_path)
        logger.info(f"  Counts  : {counts_path}")
        logger.info(f"  Metadata: {metadata_path}")

        # ── Step 2: Run R script ──────────────────────────────
        logger.info("Running DESeq2 R script...")
        self._run_r_script(
            counts_path   = counts_path,
            metadata_path = metadata_path,
            condition_col = condition_col,
        )

        # ── Step 3: Load results ──────────────────────────────
        results_path = self.output_dir / "de_results.csv"

        if not results_path.exists():
            raise RuntimeError(
                f"DESeq2 results not found at {results_path}. "
                f"Check the R script output for errors."
            )

        results = pd.read_csv(results_path)
        logger.info(f"Loaded {len(results):,} results from DESeq2")

        return results

    def run_from_files(self, counts_path   : Path,
                       metadata_path : Path,
                       condition_col : str = "condition") -> pd.DataFrame:
        """
        Run DESeq2 directly from CSV files on disk.

        Useful when you already have the CSV files saved.

        Args:
            counts_path  : Path to count matrix CSV
            metadata_path: Path to metadata CSV
            condition_col: Condition column name

        Returns:
            DataFrame with DE results
        """
        counts_path   = Path(counts_path)
        metadata_path = Path(metadata_path)

        if not counts_path.exists():
            raise FileNotFoundError(f"Counts file not found: {counts_path}")
        if not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

        logger.info("Running DESeq2 from existing files...")
        self._run_r_script(counts_path, metadata_path, condition_col)

        results_path = self.output_dir / "de_results.csv"
        return pd.read_csv(results_path)

    # ──────────────────────────────────────────────────────────
    # PRIVATE HELPERS
    # ──────────────────────────────────────────────────────────

    def _run_r_script(self, counts_path   : Path,
                      metadata_path : Path,
                      condition_col : str) -> None:
        """
        Call de_analysis.R via subprocess.

        Args:
            counts_path  : Path to count matrix CSV
            metadata_path: Path to metadata CSV
            condition_col: Condition column name
        """
        if not self.r_script.exists():
            raise FileNotFoundError(
                f"R script not found: {self.r_script}\n"
                f"Make sure de_analysis.R is in src/analysis/"
            )

        cmd = [
            RSCRIPT_PATH,
            str(self.r_script),
            str(counts_path),
            str(metadata_path),
            str(self.output_dir),
            condition_col,
        ]

        logger.info(f"Running: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output = True,
            text           = True,
        )

        # Log R output so we can see what DESeq2 is doing
        if result.stdout:
            for line in result.stdout.strip().split("\n"):
                logger.info(f"R: {line}")

        # R warnings go to stderr — log them but don't always fail
        if result.stderr:
            for line in result.stderr.strip().split("\n"):
                if line.strip():
                    logger.debug(f"R stderr: {line}")

        if result.returncode != 0:
            logger.error(f"R script failed with code {result.returncode}")
            logger.error(f"stderr:\n{result.stderr}")
            raise RuntimeError(
                f"DESeq2 R script failed. "
                f"Return code: {result.returncode}\n"
                f"Last error: {result.stderr[-500:]}"
            )

        logger.info("R script completed successfully")

    # ──────────────────────────────────────────────────────────
    # RESULTS HELPERS
    # ──────────────────────────────────────────────────────────

    def get_significant_genes(self, results     : pd.DataFrame,
                               pval_cutoff  : float = 0.05,
                               lfc_cutoff   : float = 1.5) -> pd.DataFrame:
        """
        Filter results to only significant DE genes.

        Args:
            results    : Full DESeq2 results DataFrame
            pval_cutoff: Adjusted p-value threshold
            lfc_cutoff : Log2 fold change threshold

        Returns:
            Filtered DataFrame with only significant genes
        """
        sig = results[
            (results["padj"] < pval_cutoff) &
            (results["log2FoldChange"].abs() > lfc_cutoff)
        ].copy()

        sig["direction"] = sig["log2FoldChange"].apply(
            lambda x: "Up" if x > 0 else "Down"
        )

        logger.info(f"Significant genes: {len(sig):,} "
                    f"(padj<{pval_cutoff}, |LFC|>{lfc_cutoff})")
        return sig

    def print_summary(self, results: pd.DataFrame,
                      pval_cutoff: float = 0.05,
                      lfc_cutoff : float = 1.5) -> None:
        """Print a summary of DE results."""
        comparisons = results["comparison"].unique() \
            if "comparison" in results.columns else ["all"]

        print(f"\n{'='*60}")
        print(f"  DESEQ2 RESULTS SUMMARY")
        print(f"{'='*60}")

        for comp in comparisons:
            if "comparison" in results.columns:
                subset = results[results["comparison"] == comp]
            else:
                subset = results

            sig  = subset[
                (subset["padj"] < pval_cutoff) &
                (subset["log2FoldChange"].abs() > lfc_cutoff)
            ]
            up   = sig[sig["log2FoldChange"] > 0]
            down = sig[sig["log2FoldChange"] < 0]

            print(f"\n  Comparison: {comp}")
            print(f"  Total genes tested  : {len(subset):>8,}")
            print(f"  Significant (total) : {len(sig):>8,}")
            print(f"  Upregulated         : {len(up):>8,}")
            print(f"  Downregulated       : {len(down):>8,}")

        print(f"\n  Output files in: {self.output_dir}")
        print(f"  - de_results.csv")
        print(f"  - volcano_plot.png")
        print(f"  - pca_plot.png")
        print(f"  - heatmap_top50.png")
        print(f"{'='*60}\n")
