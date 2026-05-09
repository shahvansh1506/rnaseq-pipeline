"""
tests/test_alignment.py — Unit tests for alignment and count parsing.

Run with:
    pytest tests/test_alignment.py -v
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from src.alignment.aligner import Aligner, AlignResult
from src.parsers.count_parser import CountParser


# ══════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════

def make_count_matrix(n_genes: int = 100,
                      n_samples: int = 6,
                      seed: int = 42) -> pd.DataFrame:
    """Create a fake count matrix for testing."""
    rng      = np.random.default_rng(seed)
    counts   = rng.integers(0, 1000, size=(n_genes, n_samples))
    genes    = [f"Gene_{i}" for i in range(n_genes)]
    samples  = [f"sample_{i}" for i in range(n_samples)]
    return pd.DataFrame(counts, index=genes, columns=samples)


def make_count_csv(path: Path, n_genes: int = 100) -> Path:
    """Write a fake count matrix CSV."""
    df = make_count_matrix(n_genes=n_genes)
    df.to_csv(path)
    return path


# ══════════════════════════════════════════════════════════════
# ALIGN RESULT TESTS
# ══════════════════════════════════════════════════════════════

class TestAlignResult:

    def test_summary_good_alignment(self):
        result = AlignResult(
            sample_name        = "sample1",
            total_reads        = 1_000_000,
            aligned_once       = 900_000,
            overall_align_rate = 95.0,
        )
        assert "GOOD" in result.summary()
        assert "95.0%" in result.summary()

    def test_summary_low_alignment(self):
        result = AlignResult(
            sample_name        = "sample1",
            total_reads        = 1_000_000,
            aligned_once       = 500_000,
            overall_align_rate = 50.0,
        )
        assert "LOW" in result.summary()

    def test_parse_hisat2_stats(self, tmp_path):
        aligner = Aligner(output_dir=tmp_path)

        # Simulate HISAT2 stderr output
        fake_stderr = """
1000000 reads; of these:
  1000000 (100.00%) were unpaired; of these:
    50000 (5.00%) aligned 0 times
    900000 (90.00%) aligned exactly 1 time
    50000 (5.00%) aligned >1 times
95.00% overall alignment rate
"""
        result = aligner._parse_hisat2_stats(
            fake_stderr, "sample1", "sample1_sorted.bam"
        )

        assert result.total_reads        == 1_000_000
        assert result.aligned_once       == 900_000
        assert result.aligned_multiple   == 50_000
        assert result.unaligned          == 50_000
        assert result.overall_align_rate == 95.0

    def test_parse_hisat2_empty_stderr(self, tmp_path):
        aligner = Aligner(output_dir=tmp_path)
        result  = aligner._parse_hisat2_stats("", "sample1", "out.bam")

        # Should return empty result without crashing
        assert result.total_reads == 0
        assert result.overall_align_rate == 0.0


# ══════════════════════════════════════════════════════════════
# COUNT PARSER TESTS
# ══════════════════════════════════════════════════════════════

class TestCountParser:

    def test_load_and_filter_basic(self, tmp_path):
        csv_path = make_count_csv(tmp_path / "counts.csv", n_genes=100)
        parser   = CountParser(min_total_count=10)
        counts   = parser.load_and_filter(csv_path)

        assert isinstance(counts, pd.DataFrame)
        assert len(counts) <= 100   # Some genes may be filtered

    def test_filter_removes_low_count_genes(self):
        parser = CountParser(min_total_count=10)

        # Create matrix where half the genes have zero counts
        counts = pd.DataFrame({
            "s1": [0, 0, 100, 200],
            "s2": [0, 0, 150, 250],
        }, index=["gene1", "gene2", "gene3", "gene4"])

        filtered = parser._filter_low_counts(counts)
        assert len(filtered) == 2
        assert "gene3" in filtered.index
        assert "gene4" in filtered.index
        assert "gene1" not in filtered.index

    def test_load_file_not_found(self):
        parser = CountParser()
        with pytest.raises(FileNotFoundError):
            parser.load_and_filter(Path("nonexistent.csv"))

    def test_build_metadata_basic(self):
        parser = CountParser()
        counts = make_count_matrix(n_samples=4)
        counts.columns = ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]

        condition_map = {
            "ctrl_1" : "control",
            "ctrl_2" : "control",
            "treat_1": "treated",
            "treat_2": "treated",
        }
        metadata = parser.build_metadata(counts, condition_map)

        assert "condition" in metadata.columns
        assert len(metadata) == 4
        assert set(metadata["condition"]) == {"control", "treated"}

    def test_build_metadata_missing_sample_raises(self):
        parser = CountParser()
        counts = make_count_matrix(n_samples=2)
        counts.columns = ["s1", "s2"]

        # Only map s1, missing s2
        with pytest.raises(ValueError, match="no condition assigned"):
            parser.build_metadata(counts, {"s1": "control"})

    def test_build_metadata_from_pattern(self):
        parser = CountParser()
        counts = make_count_matrix(n_samples=4)
        counts.columns = ["basal_virgin_1", "basal_virgin_2",
                          "luminal_virgin_1", "luminal_virgin_2"]

        metadata = parser.build_metadata_from_pattern(counts)
        conditions = set(metadata["condition"].tolist())
        assert "basal_virgin" in conditions
        assert "luminal_virgin" in conditions

    def test_validate_passes(self):
        parser   = CountParser()
        counts   = make_count_matrix(n_samples=4)
        counts.columns = ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]
        metadata = parser.build_metadata(counts, {
            "ctrl_1": "control", "ctrl_2": "control",
            "treat_1": "treated", "treat_2": "treated",
        })
        assert parser.validate(counts, metadata) is True

    def test_validate_fails_on_mismatch(self):
        parser   = CountParser()
        counts   = make_count_matrix(n_samples=2)
        counts.columns = ["s1", "s2"]

        # Metadata has different sample names
        metadata = pd.DataFrame(
            {"condition": ["control", "treated"]},
            index=["wrong1", "wrong2"]
        )
        with pytest.raises(ValueError):
            parser.validate(counts, metadata)

    def test_normalize_cpm(self):
        parser = CountParser()
        counts = pd.DataFrame({
            "s1": [100, 900],
            "s2": [200, 800],
        }, index=["gene1", "gene2"])

        cpm = parser.normalize_cpm(counts)

        # s1 total = 1000, gene1 CPM = 100/1000 * 1e6 = 100,000
        assert abs(cpm.loc["gene1", "s1"] - 100_000) < 1

    def test_normalize_log_cpm_no_negatives(self):
        parser  = CountParser()
        counts  = make_count_matrix(n_genes=50, n_samples=4)
        log_cpm = parser.normalize_log_cpm(counts)

        # log2(CPM + 1) should never be negative
        assert (log_cpm >= 0).all().all()
