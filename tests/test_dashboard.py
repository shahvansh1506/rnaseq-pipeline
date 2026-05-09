"""
tests/test_dashboard.py — Unit tests for validation and dashboard logic.

Run with:
    pytest tests/test_dashboard.py -v
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock
from src.validation.ncbi_validator import NCBIValidator, GeneValidation


# ══════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════

def make_de_results(n: int = 50) -> pd.DataFrame:
    """Create fake DE results DataFrame."""
    rng = np.random.default_rng(42)
    return pd.DataFrame({
        "gene"           : [f"Gene_{i}" for i in range(n)],
        "log2FoldChange" : rng.uniform(-4, 4, n),
        "padj"           : rng.uniform(0.001, 0.5, n),
        "baseMean"       : rng.uniform(10, 1000, n),
        "comparison"     : "treated_vs_control",
    })


# ══════════════════════════════════════════════════════════════
# GENE VALIDATION TESTS
# ══════════════════════════════════════════════════════════════

class TestGeneValidation:

    def test_summary_valid_gene(self):
        v = GeneValidation(
            symbol      = "Csn2",
            validated   = True,
            description = "casein beta",
        )
        assert "VALID" in v.summary()
        assert "Csn2" in v.summary()

    def test_summary_invalid_gene(self):
        v = GeneValidation(symbol="FakeGene123", validated=False)
        assert "NOT FOUND" in v.summary()

    def test_summary_truncates_long_description(self):
        v = GeneValidation(
            symbol      = "Gene1",
            validated   = True,
            description = "A" * 100,
        )
        assert "..." in v.summary()

    def test_default_not_validated(self):
        v = GeneValidation(symbol="Gene1")
        assert v.validated is False
        assert v.ncbi_id   == ""


# ══════════════════════════════════════════════════════════════
# NCBI VALIDATOR TESTS
# ══════════════════════════════════════════════════════════════

class TestNCBIValidator:

    def _make_validator(self) -> NCBIValidator:
        """Create validator with mocked Entrez."""
        validator        = NCBIValidator.__new__(NCBIValidator)
        validator.email  = "test@test.com"

        mock_entrez           = MagicMock()
        mock_entrez.email     = "test@test.com"
        validator.Entrez      = mock_entrez
        return validator

    def test_validate_found_gene(self):
        validator = self._make_validator()

        # Mock NCBI returning a gene ID
        mock_handle = MagicMock()
        validator.Entrez.esearch.return_value = mock_handle
        validator.Entrez.read.return_value    = {"IdList": ["12345"]}

        mock_fetch = MagicMock()
        mock_fetch.read.return_value = "Csn2\tcasein beta\tMus musculus"
        validator.Entrez.efetch.return_value = mock_fetch

        with patch("time.sleep"):
            result = validator._query_ncbi("Csn2", "Mus musculus")

        assert result.validated is True
        assert result.ncbi_id   == "12345"

    def test_validate_missing_gene(self):
        validator = self._make_validator()

        # Mock NCBI returning no results
        mock_handle = MagicMock()
        validator.Entrez.esearch.return_value = mock_handle
        validator.Entrez.read.return_value    = {"IdList": []}

        with patch("time.sleep"):
            result = validator._query_ncbi("FakeGene999", "Mus musculus")

        assert result.validated is False
        assert "not found" in result.error.lower()

    def test_validate_genes_returns_list(self):
        validator = self._make_validator()

        mock_handle = MagicMock()
        validator.Entrez.esearch.return_value = mock_handle
        validator.Entrez.read.return_value    = {"IdList": ["111"]}

        mock_fetch = MagicMock()
        mock_fetch.read.return_value = "Gene1\tdescription"
        validator.Entrez.efetch.return_value = mock_fetch

        with patch("time.sleep"):
            results = validator.validate_genes(
                ["Gene1", "Gene2"], "Mus musculus"
            )

        assert len(results) == 2
        assert all(isinstance(r, GeneValidation) for r in results)

    def test_save_report(self, tmp_path):
        validator   = self._make_validator()
        validations = [
            GeneValidation("Csn2",     validated=True,
                           ncbi_id="12345", description="casein beta"),
            GeneValidation("FakeGene", validated=False,
                           error="not found"),
        ]

        path = validator.save_report(
            validations, tmp_path / "report.csv"
        )

        assert path.exists()
        df = pd.read_csv(path)
        assert len(df) == 2
        assert "gene" in df.columns
        assert "validated" in df.columns

    def test_print_report_runs(self, capsys):
        validator   = self._make_validator()
        validations = [
            GeneValidation("Csn2", validated=True,
                           description="casein beta"),
            GeneValidation("FakeGene", validated=False),
        ]
        validator.print_report(validations)
        captured = capsys.readouterr()
        assert "NCBI VALIDATION REPORT" in captured.out
        assert "Csn2" in captured.out

    def test_validate_deseq2_results(self):
        validator   = self._make_validator()
        results_df  = make_de_results(20)

        # Make some genes significant
        results_df.loc[:4, "padj"] = 0.001

        mock_handle = MagicMock()
        validator.Entrez.esearch.return_value = mock_handle
        validator.Entrez.read.return_value    = {"IdList": ["111"]}

        mock_fetch = MagicMock()
        mock_fetch.read.return_value = "gene info"
        validator.Entrez.efetch.return_value = mock_fetch

        with patch("time.sleep"):
            merged = validator.validate_deseq2_results(
                results_df, max_genes=5
            )

        assert "validated" in merged.columns
        assert len(merged) == 5
