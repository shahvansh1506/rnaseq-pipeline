"""
ncbi_validator.py — Validates DE genes against NCBI database.

What this file does:
- Takes a list of gene symbols from DESeq2 results
- Queries NCBI Gene database to verify each gene exists
- Retrieves gene descriptions and aliases
- Flags any genes that couldn't be validated
- Saves a validation report

Why validate?
- Confirms your DE genes are real, annotated genes
- Catches annotation errors or outdated gene symbols
- Adds biological context (gene descriptions) to your results

Usage:
    validator = NCBIValidator(email="your@email.com")
    report    = validator.validate_genes(
        gene_list = ["Csn1s2a", "Wap", "Gapdh"],
        organism  = "Mus musculus"
    )
"""

import time
import logging
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field

from src.config import ENTREZ_EMAIL, RESULTS_DIR, LOG_FORMAT, LOG_LEVEL

# ── Setup logger ──────────────────────────────────────────────
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════
# DATA CLASS — validation result for one gene
# ══════════════════════════════════════════════════════════════

@dataclass
class GeneValidation:
    """Stores NCBI validation result for a single gene."""
    symbol      : str
    validated   : bool  = False
    ncbi_id     : str   = ""
    description : str   = ""
    organism    : str   = ""
    aliases     : list  = field(default_factory=list)
    error       : str   = ""

    def summary(self) -> str:
        status = "✅ VALID" if self.validated else "❌ NOT FOUND"
        desc   = self.description[:50] + "..." \
                 if len(self.description) > 50 else self.description
        return f"{status} | {self.symbol:<15} | {desc}"


# ══════════════════════════════════════════════════════════════
# MAIN CLASS
# ══════════════════════════════════════════════════════════════

class NCBIValidator:
    """
    Validates gene symbols against the NCBI Gene database.

    Uses Biopython's Entrez module to query NCBI.
    NCBI requires an email address for API access.

    Usage:
        validator = NCBIValidator(email="your@email.com")
        report    = validator.validate_genes(["Csn2", "Wap", "Gapdh"])
    """

    # NCBI rate limit: max 3 requests/second without API key
    # We wait 0.4s between requests to be safe
    REQUEST_DELAY = 0.4

    def __init__(self, email: str = None):
        """
        Args:
            email: Your email for NCBI Entrez API access.
                   Required by NCBI terms of service.
        """
        self.email = email or ENTREZ_EMAIL

        # Import here so missing biopython gives a clear error
        try:
            from Bio import Entrez
            Entrez.email = self.email
            self.Entrez  = Entrez
            logger.info(f"NCBI Entrez configured with email: {self.email}")
        except ImportError:
            raise ImportError(
                "Biopython not installed. Run: pip install biopython"
            )

    # ──────────────────────────────────────────────────────────
    # MAIN VALIDATION METHOD
    # ──────────────────────────────────────────────────────────

    def validate_genes(self, gene_list : list,
                       organism   : str = "Mus musculus") -> list:
        """
        Validate a list of gene symbols against NCBI.

        Args:
            gene_list: List of gene symbols e.g. ["Csn2", "Wap"]
            organism : Scientific name e.g. "Mus musculus" or
                       "Homo sapiens"

        Returns:
            List of GeneValidation objects
        """
        logger.info(f"Validating {len(gene_list)} genes "
                    f"in {organism}...")

        results = []
        for i, gene in enumerate(gene_list, 1):
            logger.info(f"  Validating {i}/{len(gene_list)}: {gene}")

            validation = self._query_ncbi(gene, organism)
            results.append(validation)

            # Respect NCBI rate limit
            time.sleep(self.REQUEST_DELAY)

        validated = sum(1 for r in results if r.validated)
        logger.info(f"Validation complete: "
                    f"{validated}/{len(results)} genes found")

        return results

    def validate_deseq2_results(self, results_df  : pd.DataFrame,
                                 organism     : str = "Mus musculus",
                                 gene_col     : str = "gene",
                                 max_genes    : int = 50) -> pd.DataFrame:
        """
        Validate top DE genes from DESeq2 results.

        Takes the most significant genes and validates them
        against NCBI, then merges validation info back.

        Args:
            results_df: DESeq2 results DataFrame
            organism  : Scientific organism name
            gene_col  : Column containing gene symbols
            max_genes : Max number of genes to validate
                        (keep low to avoid rate limiting)

        Returns:
            DataFrame with validation info added
        """
        # Get top significant genes
        sig = results_df[
            results_df["padj"] < 0.05
        ].nsmallest(max_genes, "padj")

        gene_list   = sig[gene_col].tolist()
        validations = self.validate_genes(gene_list, organism)

        # Build validation DataFrame
        val_df = pd.DataFrame([{
            "gene"       : v.symbol,
            "validated"  : v.validated,
            "ncbi_id"    : v.ncbi_id,
            "description": v.description,
            "aliases"    : ", ".join(v.aliases[:3]),
        } for v in validations])

        # Merge back with DE results
        merged = sig.merge(val_df, on="gene", how="left")
        return merged

    # ──────────────────────────────────────────────────────────
    # NCBI QUERY
    # ──────────────────────────────────────────────────────────

    def _query_ncbi(self, gene_symbol: str,
                    organism    : str) -> GeneValidation:
        """
        Query NCBI Gene database for a single gene symbol.

        Args:
            gene_symbol: Gene symbol e.g. "Csn2"
            organism   : Scientific name e.g. "Mus musculus"

        Returns:
            GeneValidation with results
        """
        validation = GeneValidation(symbol=gene_symbol)

        try:
            # Search NCBI Gene database
            search_term = (f"{gene_symbol}[Gene Name] AND "
                           f"{organism}[Organism]")

            handle  = self.Entrez.esearch(db="gene", term=search_term)
            record  = self.Entrez.read(handle)
            handle.close()

            id_list = record.get("IdList", [])

            if not id_list:
                validation.error = "Gene not found in NCBI"
                return validation

            # Fetch details for the first matching gene
            gene_id = id_list[0]
            handle  = self.Entrez.efetch(
                db     = "gene",
                id     = gene_id,
                rettype= "gene_table",
                retmode= "text"
            )
            content = handle.read()
            handle.close()

            # Mark as validated and store basic info
            validation.validated = True
            validation.ncbi_id   = gene_id
            validation.organism  = organism

            # Parse description from content
            lines = content.split("\n") if isinstance(content, str) \
                    else content.decode().split("\n")

            for line in lines[:10]:
                if gene_symbol.lower() in line.lower() and len(line) > 10:
                    validation.description = line.strip()[:200]
                    break

            if not validation.description:
                validation.description = f"NCBI Gene ID: {gene_id}"

        except Exception as e:
            validation.error = str(e)
            logger.warning(f"NCBI query failed for {gene_symbol}: {e}")

        return validation

    # ──────────────────────────────────────────────────────────
    # REPORTING
    # ──────────────────────────────────────────────────────────

    def print_report(self, validations: list) -> None:
        """Print validation report to console."""
        validated = [v for v in validations if v.validated]
        not_found = [v for v in validations if not v.validated]

        print(f"\n{'='*65}")
        print(f"  NCBI VALIDATION REPORT")
        print(f"  Total: {len(validations)}  |  "
              f"Valid: {len(validated)}  |  "
              f"Not found: {len(not_found)}")
        print(f"{'='*65}")

        for v in validations:
            print(f"  {v.summary()}")
            if v.error:
                print(f"    Error: {v.error}")

        if not_found:
            print(f"\n  Genes not found in NCBI:")
            for v in not_found:
                print(f"    - {v.symbol}: {v.error}")

        print(f"{'='*65}\n")

    def save_report(self, validations : list,
                    output_path  : Path = None) -> Path:
        """
        Save validation report to CSV.

        Args:
            validations: List of GeneValidation objects
            output_path: Where to save. Auto-named if None.

        Returns:
            Path to saved CSV
        """
        output_path = output_path or RESULTS_DIR / "gene_validation.csv"

        rows = [{
            "gene"       : v.symbol,
            "validated"  : v.validated,
            "ncbi_id"    : v.ncbi_id,
            "description": v.description,
            "organism"   : v.organism,
            "error"      : v.error,
        } for v in validations]

        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        logger.info(f"Validation report saved: {output_path}")
        return output_path
