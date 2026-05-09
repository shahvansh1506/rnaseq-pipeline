"""
count_parser.py — Loads, cleans and prepares count matrices for DESeq2.

What this file does:
- Loads count matrix CSV files
- Filters out low-count genes (noise)
- Normalises counts for visualization (NOT for DESeq2 — DESeq2 needs raw counts)
- Builds sample metadata table required by DESeq2
- Validates that count matrix and metadata match

Why filter low-count genes?
- Genes with very few counts across all samples are noise
- They slow down DESeq2 without adding useful information
- Standard practice: remove genes with fewer than 10 total counts

Usage:
    parser   = CountParser()
    counts   = parser.load_and_filter("data/counts/count_matrix.csv")
    metadata = parser.build_metadata(counts, condition_map)
    parser.validate(counts, metadata)
"""

import logging
import pandas as pd
import numpy as np
from pathlib import Path

from src.config import COUNTS_DIR, LOG_FORMAT, LOG_LEVEL

# ── Setup logger ──────────────────────────────────────────────
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


class CountParser:
    """
    Loads and prepares count matrices for differential expression.

    The count matrix format required by DESeq2:
    - Rows   = genes (Ensembl IDs or gene symbols)
    - Columns = samples
    - Values = raw integer counts (NOT normalized)

    Usage:
        parser   = CountParser()
        counts   = parser.load_and_filter("count_matrix.csv")
        metadata = parser.build_metadata(counts, condition_map)
    """

    def __init__(self, min_total_count: int = 10):
        """
        Args:
            min_total_count: Minimum total counts across all samples
                             to keep a gene. Default is 10.
        """
        self.min_total_count = min_total_count

    # ──────────────────────────────────────────────────────────
    # LOADING
    # ──────────────────────────────────────────────────────────

    def load_and_filter(self, filepath: Path) -> pd.DataFrame:
        """
        Load count matrix and filter out low-count genes.

        Args:
            filepath: Path to count matrix CSV

        Returns:
            Filtered DataFrame: genes x samples
        """
        filepath = Path(filepath)

        if not filepath.exists():
            raise FileNotFoundError(f"Count matrix not found: {filepath}")

        logger.info(f"Loading count matrix: {filepath.name}")
        counts = pd.read_csv(filepath, index_col=0)

        # Make sure all values are integers
        counts = counts.astype(int)

        before = len(counts)
        logger.info(f"Loaded {before:,} genes x {counts.shape[1]} samples")

        # Filter low count genes
        counts = self._filter_low_counts(counts)
        after  = len(counts)

        logger.info(f"After filtering: {after:,} genes kept "
                    f"({before - after:,} removed)")

        return counts

    def _filter_low_counts(self, counts: pd.DataFrame) -> pd.DataFrame:
        """
        Remove genes with very low total counts.

        A gene needs at least min_total_count reads summed
        across ALL samples to be kept.

        Args:
            counts: Raw count matrix

        Returns:
            Filtered count matrix
        """
        # Sum counts across all samples for each gene
        gene_totals = counts.sum(axis=1)

        # Keep only genes above threshold
        keep_mask = gene_totals >= self.min_total_count
        filtered  = counts[keep_mask]

        logger.info(f"Filtering at min total count = {self.min_total_count}")
        logger.info(f"Kept {keep_mask.sum():,} / {len(counts):,} genes")

        return filtered

    # ──────────────────────────────────────────────────────────
    # METADATA
    # ──────────────────────────────────────────────────────────

    def build_metadata(self, counts       : pd.DataFrame,
                       condition_map : dict) -> pd.DataFrame:
        """
        Build sample metadata table for DESeq2.

        DESeq2 needs a metadata table that maps each sample
        to its experimental condition (e.g. control vs treated).

        Args:
            counts       : Count matrix (columns = sample names)
            condition_map: Dict mapping sample name -> condition
                          e.g. {"sample1": "control",
                                "sample2": "treated"}

        Returns:
            Metadata DataFrame with columns: sample, condition

        Example:
            condition_map = {
                "basal_virgin_1"   : "basal_virgin",
                "basal_pregnant_1" : "basal_pregnant",
                "luminal_virgin_1" : "luminal_virgin",
            }
        """
        samples = counts.columns.tolist()

        # Check all samples have a condition assigned
        missing = [s for s in samples if s not in condition_map]
        if missing:
            raise ValueError(
                f"These samples have no condition assigned: {missing}\n"
                f"Please add them to condition_map."
            )

        metadata = pd.DataFrame({
            "sample"   : samples,
            "condition": [condition_map[s] for s in samples],
        })
        metadata = metadata.set_index("sample")

        logger.info(f"Built metadata for {len(metadata)} samples")
        logger.info(f"Conditions: {metadata['condition'].unique().tolist()}")

        return metadata

    def build_metadata_from_pattern(self, counts   : pd.DataFrame,
                                    separator : str = "_",
                                    parts     : list = None) -> pd.DataFrame:
        """
        Auto-build metadata by parsing sample names.

        Useful when sample names follow a pattern like:
        "basal_virgin_rep1" -> condition = "basal_virgin"

        Args:
            counts   : Count matrix
            separator: Character splitting name parts (default "_")
            parts    : Which parts of the name to use as condition
                       e.g. [0, 1] uses first two parts

        Returns:
            Metadata DataFrame
        """
        samples = counts.columns.tolist()

        if parts is None:
            # Use everything except the last part (usually replicate number)
            conditions = []
            for s in samples:
                name_parts = s.split(separator)
                condition  = separator.join(name_parts[:-1])
                conditions.append(condition)
        else:
            conditions = []
            for s in samples:
                name_parts = s.split(separator)
                condition  = separator.join([name_parts[i] for i in parts
                                             if i < len(name_parts)])
                conditions.append(condition)

        metadata = pd.DataFrame({
            "sample"   : samples,
            "condition": conditions,
        })
        metadata = metadata.set_index("sample")

        logger.info(f"Auto-built metadata: "
                    f"{metadata['condition'].unique().tolist()}")
        return metadata

    # ──────────────────────────────────────────────────────────
    # VALIDATION
    # ──────────────────────────────────────────────────────────

    def validate(self, counts  : pd.DataFrame,
                 metadata: pd.DataFrame) -> bool:
        """
        Validate that count matrix and metadata are compatible.

        Checks:
        1. Sample names match between counts and metadata
        2. No negative counts
        3. No NaN values
        4. At least 2 conditions for comparison

        Args:
            counts  : Count matrix
            metadata: Sample metadata

        Returns:
            True if valid, raises ValueError if not
        """
        logger.info("Validating count matrix and metadata...")
        errors = []

        # Check 1: Sample names match
        count_samples    = set(counts.columns)
        metadata_samples = set(metadata.index)

        if count_samples != metadata_samples:
            in_counts_not_meta = count_samples - metadata_samples
            in_meta_not_counts = metadata_samples - count_samples

            if in_counts_not_meta:
                errors.append(f"In counts but not metadata: "
                               f"{in_counts_not_meta}")
            if in_meta_not_counts:
                errors.append(f"In metadata but not counts: "
                               f"{in_meta_not_counts}")

        # Check 2: No negative counts
        if (counts < 0).any().any():
            errors.append("Count matrix contains negative values")

        # Check 3: No NaN values
        if counts.isnull().any().any():
            errors.append("Count matrix contains NaN values")

        # Check 4: At least 2 conditions
        if "condition" in metadata.columns:
            n_conditions = metadata["condition"].nunique()
            if n_conditions < 2:
                errors.append(f"Need at least 2 conditions for DE analysis, "
                               f"found {n_conditions}")

        if errors:
            for err in errors:
                logger.error(f"Validation error: {err}")
            raise ValueError(f"Validation failed:\n" +
                              "\n".join(errors))

        logger.info("Validation passed!")
        return True

    # ──────────────────────────────────────────────────────────
    # NORMALIZATION (for visualization only)
    # ──────────────────────────────────────────────────────────

    def normalize_cpm(self, counts: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize counts to CPM (Counts Per Million).

        CPM is used ONLY for visualization and comparison.
        Always use RAW counts for DESeq2!

        CPM formula:
            CPM = (count / total_counts_in_sample) * 1,000,000

        This corrects for differences in library size between samples.

        Args:
            counts: Raw count matrix

        Returns:
            CPM normalized DataFrame
        """
        # Total counts per sample (library size)
        library_sizes = counts.sum(axis=0)

        # Divide each column by its library size, multiply by 1M
        cpm = counts.divide(library_sizes, axis=1) * 1_000_000

        logger.info("Normalized counts to CPM")
        return cpm

    def normalize_log_cpm(self, counts  : pd.DataFrame,
                          pseudocount: int = 1) -> pd.DataFrame:
        """
        Log2(CPM + pseudocount) normalization.

        Adding a pseudocount (1) prevents log(0) errors.
        Log transformation makes the data more normally distributed.
        Used for PCA and heatmap visualization.

        Args:
            counts    : Raw count matrix
            pseudocount: Value added before log transform

        Returns:
            Log2 CPM normalized DataFrame
        """
        cpm     = self.normalize_cpm(counts)
        log_cpm = np.log2(cpm + pseudocount)

        logger.info("Normalized counts to log2(CPM + 1)")
        return log_cpm

    # ──────────────────────────────────────────────────────────
    # SUMMARY
    # ──────────────────────────────────────────────────────────

    def print_summary(self, counts  : pd.DataFrame,
                      metadata: pd.DataFrame = None) -> None:
        """Print a summary of the count matrix and metadata."""
        print(f"\n{'='*55}")
        print(f"  COUNT MATRIX SUMMARY")
        print(f"{'='*55}")
        print(f"  Genes          : {counts.shape[0]:>10,}")
        print(f"  Samples        : {counts.shape[1]:>10,}")
        print(f"  Total counts   : {counts.values.sum():>10,.0f}")
        print(f"  Zero-count genes: {(counts.sum(axis=1) == 0).sum():>9,}")
        print(f"  Min count      : {counts.values.min():>10,}")
        print(f"  Max count      : {counts.values.max():>10,}")

        if metadata is not None and "condition" in metadata.columns:
            print(f"\n  Conditions:")
            for cond, grp in metadata.groupby("condition"):
                print(f"    {cond:<30} {len(grp):>3} samples")

        print(f"{'='*55}\n")
