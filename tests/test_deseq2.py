"""
tests/test_deseq2.py — Unit tests for DESeq2 runner.

These tests validate the Python side of the DESeq2 wrapper
without actually running R (which requires DESeq2 installed).

Run with:
    pytest tests/test_deseq2.py -v
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock
from src.analysis.run_deseq2 import DESeq2Runner


# ══════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════

def make_counts(n_genes: int = 50, n_samples: int = 6) -> pd.DataFrame:
    """Create a fake count matrix."""
    rng     = np.random.default_rng(42)
    data    = rng.integers(10, 1000, size=(n_genes, n_samples))
    genes   = [f"Gene_{i}" for i in range(n_genes)]
    samples = (["ctrl_1", "ctrl_2", "ctrl_3",
                "treat_1", "treat_2", "treat_3"])
    return pd.DataFrame(data, index=genes, columns=samples)


def make_metadata() -> pd.DataFrame:
    """Create fake sample metadata."""
    return pd.DataFrame(
        {"condition": ["control"] * 3 + ["treated"] * 3},
        index=["ctrl_1", "ctrl_2", "ctrl_3",
               "treat_1", "treat_2", "treat_3"]
    )


def make_fake_results(n_genes: int = 100) -> pd.DataFrame:
    """Create fake DESeq2 results for testing."""
    rng = np.random.default_rng(42)
    return pd.DataFrame({
        "gene"           : [f"Gene_{i}" for i in range(n_genes)],
        "baseMean"       : rng.uniform(10, 1000, n_genes),
        "log2FoldChange" : rng.uniform(-4, 4, n_genes),
        "lfcSE"          : rng.uniform(0.1, 0.5, n_genes),
        "stat"           : rng.uniform(-5, 5, n_genes),
        "pvalue"         : rng.uniform(0, 1, n_genes),
        "padj"           : rng.uniform(0, 1, n_genes),
        "comparison"     : "treated_vs_control",
    })


# ══════════════════════════════════════════════════════════════
# DESEQ2 RUNNER TESTS
# ══════════════════════════════════════════════════════════════

class TestDESeq2Runner:

    def test_init_creates_output_dir(self, tmp_path):
        runner = DESeq2Runner(output_dir=tmp_path / "deseq2_out")
        assert runner.output_dir.exists()

    def test_run_saves_input_csvs(self, tmp_path):
        """Test that run() saves count matrix and metadata to CSV."""
        runner   = DESeq2Runner(output_dir=tmp_path)
        counts   = make_counts()
        metadata = make_metadata()

        # Mock _run_r_script so we don't need R installed
        with patch.object(runner, "_run_r_script"):
            # Also mock results file
            fake_results = make_fake_results()
            fake_results.to_csv(tmp_path / "de_results.csv", index=False)

            runner.run(counts, metadata, "condition")

        # Check CSV files were saved
        assert (tmp_path / "input_counts.csv").exists()
        assert (tmp_path / "input_metadata.csv").exists()

    def test_run_loads_results(self, tmp_path):
        """Test that run() correctly loads the results CSV."""
        runner       = DESeq2Runner(output_dir=tmp_path)
        counts       = make_counts()
        metadata     = make_metadata()
        fake_results = make_fake_results(50)

        with patch.object(runner, "_run_r_script"):
            fake_results.to_csv(tmp_path / "de_results.csv", index=False)
            results = runner.run(counts, metadata, "condition")

        assert isinstance(results, pd.DataFrame)
        assert len(results) == 50
        assert "gene" in results.columns

    def test_run_raises_if_no_results(self, tmp_path):
        """Test that run() raises error if R script produces no results."""
        runner   = DESeq2Runner(output_dir=tmp_path)
        counts   = make_counts()
        metadata = make_metadata()

        with patch.object(runner, "_run_r_script"):
            # Don't create de_results.csv — should raise
            with pytest.raises(RuntimeError, match="DESeq2 results not found"):
                runner.run(counts, metadata, "condition")

    def test_r_script_path_exists(self):
        """Test that de_analysis.R exists where expected."""
        runner = DESeq2Runner.__new__(DESeq2Runner)
        runner.r_script = Path("src/analysis/de_analysis.R")
        assert runner.r_script.exists(), \
            "de_analysis.R not found — make sure it's in src/analysis/"

    def test_run_r_script_raises_on_failure(self, tmp_path):
        """Test that _run_r_script raises on non-zero return code."""
        runner          = DESeq2Runner(output_dir=tmp_path)
        runner.r_script = Path("src/analysis/de_analysis.R")

        mock_result           = MagicMock()
        mock_result.returncode = 1
        mock_result.stdout     = ""
        mock_result.stderr     = "Error: something went wrong"

        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(RuntimeError, match="DESeq2 R script failed"):
                runner._run_r_script(
                    Path("counts.csv"),
                    Path("meta.csv"),
                    "condition"
                )

    def test_get_significant_genes(self, tmp_path):
        """Test filtering results to significant genes."""
        runner  = DESeq2Runner(output_dir=tmp_path)
        results = make_fake_results(200)

        # Manually set some genes to be significant
        results.loc[:9, "padj"]           = 0.001
        results.loc[:9, "log2FoldChange"] = 2.5

        sig = runner.get_significant_genes(
            results, pval_cutoff=0.05, lfc_cutoff=1.5
        )

        assert len(sig) >= 10
        assert all(sig["padj"] < 0.05)
        assert all(sig["log2FoldChange"].abs() > 1.5)

    def test_get_significant_genes_direction(self, tmp_path):
        """Test that direction (Up/Down) is correctly assigned."""
        runner  = DESeq2Runner(output_dir=tmp_path)
        results = pd.DataFrame({
            "gene"           : ["GeneA", "GeneB"],
            "log2FoldChange" : [2.0, -2.0],
            "padj"           : [0.01, 0.01],
            "comparison"     : ["treated_vs_control"] * 2,
        })

        sig = runner.get_significant_genes(
            results, pval_cutoff=0.05, lfc_cutoff=1.5
        )

        assert sig[sig["gene"] == "GeneA"]["direction"].values[0] == "Up"
        assert sig[sig["gene"] == "GeneB"]["direction"].values[0] == "Down"

    def test_run_from_files_raises_if_missing(self, tmp_path):
        """Test run_from_files raises if files don't exist."""
        runner = DESeq2Runner(output_dir=tmp_path)

        with pytest.raises(FileNotFoundError):
            runner.run_from_files(
                Path("nonexistent_counts.csv"),
                Path("nonexistent_meta.csv"),
            )
