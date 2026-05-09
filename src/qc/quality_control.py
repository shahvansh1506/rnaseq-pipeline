"""
quality_control.py — Quality control for raw RNA-seq reads.

What this file does:
- Runs FastQC on raw FASTQ files to generate QC reports
- Parses FastQC output to extract key QC metrics
- Runs MultiQC to combine reports from all samples
- Flags samples that fail QC thresholds

Usage:
    qc = QualityController()
    qc.run_fastqc("data/raw/sample.fastq.gz")
    report = qc.parse_fastqc_zip("results/qc/sample_fastqc.zip")
    qc.run_multiqc("results/qc/")
"""

import subprocess
import zipfile
import logging
import csv
from pathlib import Path
from dataclasses import dataclass, field

from src.config import (
    FASTQC_PATH, MULTIQC_PATH,
    RESULTS_DIR, THREADS, LOG_FORMAT, LOG_LEVEL
)

# ── Setup logger ──────────────────────────────────────────────
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════
# DATA CLASS — holds QC results for one sample
# ══════════════════════════════════════════════════════════════

@dataclass
class QCReport:
    """
    Stores parsed FastQC metrics for a single sample.
    FastQC grades each metric as PASS, WARN, or FAIL.
    """
    sample_name             : str
    total_sequences         : int  = 0
    poor_quality            : int  = 0
    sequence_length         : str  = ""
    gc_percent              : int  = 0
    encoding                : str  = ""

    # Per-module PASS/WARN/FAIL grades from FastQC
    basic_statistics        : str  = "UNKNOWN"
    per_base_quality        : str  = "UNKNOWN"
    per_sequence_quality    : str  = "UNKNOWN"
    per_base_gc_content     : str  = "UNKNOWN"
    adapter_content         : str  = "UNKNOWN"
    sequence_duplication    : str  = "UNKNOWN"

    # Overall verdict
    passes_qc               : bool = True
    warnings                : list = field(default_factory=list)

    def summary(self) -> str:
        """Return a one-line summary string."""
        status = "PASS" if self.passes_qc else "FAIL"
        return (f"[{status}] {self.sample_name:<30} | "
                f"Reads: {self.total_sequences:>10,} | "
                f"GC: {self.gc_percent:>3}% | "
                f"Quality: {self.per_base_quality}")


# ══════════════════════════════════════════════════════════════
# MAIN CLASS
# ══════════════════════════════════════════════════════════════

class QualityController:
    """
    Runs FastQC and MultiQC on FASTQ files and parses results.

    Usage:
        qc = QualityController(output_dir=Path("results/qc"))
        qc.run_fastqc(Path("data/raw/sample.fastq.gz"))
        report = qc.parse_fastqc_zip(Path("results/qc/sample_fastqc.zip"))
        qc.run_multiqc(Path("results/qc"))
    """

    def __init__(self, output_dir: Path = None):
        self.output_dir = output_dir or RESULTS_DIR / "qc"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.reports = []

    # ──────────────────────────────────────────────────────────
    # FASTQC
    # ──────────────────────────────────────────────────────────

    def run_fastqc(self, fastq_path: Path,
                   threads: int = THREADS) -> Path:
        """
        Run FastQC on a single FASTQ file.

        FastQC generates two output files:
        - sample_fastqc.html  -> human readable report
        - sample_fastqc.zip   -> machine readable data we parse

        Args:
            fastq_path: Path to FASTQ or FASTQ.gz file
            threads   : Number of CPU threads to use

        Returns:
            Path to the output directory
        """
        fastq_path = Path(fastq_path)

        if not fastq_path.exists():
            raise FileNotFoundError(f"FASTQ not found: {fastq_path}")

        logger.info(f"Running FastQC on: {fastq_path.name}")

        cmd = [
            FASTQC_PATH,
            str(fastq_path),
            "--outdir", str(self.output_dir),
            "--threads", str(threads),
            "--quiet",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"FastQC failed:\n{result.stderr}")
            raise RuntimeError(f"FastQC failed for {fastq_path.name}")

        logger.info(f"FastQC complete. Reports saved to: {self.output_dir}")
        return self.output_dir

    def run_fastqc_batch(self, fastq_dir: Path) -> list:
        """
        Run FastQC on all FASTQ files in a directory.

        Args:
            fastq_dir: Directory containing FASTQ files

        Returns:
            List of QCReport objects for all samples
        """
        fastq_dir   = Path(fastq_dir)
        fastq_files = (list(fastq_dir.glob("*.fastq.gz")) +
                       list(fastq_dir.glob("*.fastq")))

        if not fastq_files:
            logger.warning(f"No FASTQ files found in {fastq_dir}")
            return []

        logger.info(f"Found {len(fastq_files)} FASTQ files to process")

        reports = []
        for i, fastq_file in enumerate(fastq_files, 1):
            logger.info(f"Processing {i}/{len(fastq_files)}: {fastq_file.name}")
            try:
                self.run_fastqc(fastq_file)
                zip_file = self.output_dir / f"{fastq_file.stem}_fastqc.zip"
                if zip_file.exists():
                    report = self.parse_fastqc_zip(zip_file)
                    reports.append(report)
            except Exception as e:
                logger.error(f"Failed to process {fastq_file.name}: {e}")

        self.reports = reports
        return reports

    # ──────────────────────────────────────────────────────────
    # PARSING FASTQC OUTPUT
    # ──────────────────────────────────────────────────────────

    def parse_fastqc_zip(self, zip_path: Path) -> QCReport:
        """
        Parse FastQC zip output to extract QC metrics.

        Args:
            zip_path: Path to sample_fastqc.zip

        Returns:
            QCReport with all metrics filled in
        """
        zip_path    = Path(zip_path)
        sample_name = zip_path.stem.replace("_fastqc", "")

        logger.info(f"Parsing FastQC report: {zip_path.name}")

        with zipfile.ZipFile(zip_path, "r") as zf:
            data_files = [f for f in zf.namelist()
                          if f.endswith("fastqc_data.txt")]

            if not data_files:
                raise ValueError(f"No fastqc_data.txt in {zip_path}")

            with zf.open(data_files[0]) as f:
                content = f.read().decode("utf-8")

        return self._parse_fastqc_data(content, sample_name)

    def _parse_fastqc_data(self, content: str,
                            sample_name: str) -> QCReport:
        """
        Parse the raw text content of fastqc_data.txt.

        The file has sections like:
            >>Basic Statistics    PASS
            #Measure    Value
            Total Sequences    1000000
            >>END_MODULE

        Args:
            content    : Full text of fastqc_data.txt
            sample_name: Name of the sample

        Returns:
            Populated QCReport object
        """
        report      = QCReport(sample_name=sample_name)
        lines       = content.split("\n")
        current_mod = None

        # Map FastQC module names to QCReport fields
        module_map = {
            "Basic Statistics"           : "basic_statistics",
            "Per base sequence quality"  : "per_base_quality",
            "Per sequence quality scores": "per_sequence_quality",
            "Per sequence GC content"    : "per_base_gc_content",
            "Adapter Content"            : "adapter_content",
            "Sequence Duplication Levels": "sequence_duplication",
        }

        for line in lines:
            line = line.strip()

            # Detect section headers e.g. ">>Basic Statistics\tPASS"
            if line.startswith(">>") and not line.startswith(">>END"):
                parts = line[2:].split("\t")
                if len(parts) == 2:
                    module_name, grade = parts[0], parts[1]
                    current_mod        = module_name

                    if module_name in module_map:
                        setattr(report, module_map[module_name], grade)

                        if grade == "WARN":
                            report.warnings.append(f"WARN: {module_name}")
                        elif grade == "FAIL":
                            report.warnings.append(f"FAIL: {module_name}")
                            report.passes_qc = False

            # Parse key metrics from Basic Statistics section
            elif current_mod == "Basic Statistics" and "\t" in line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    key, val = parts[0], parts[1]
                    if key == "Total Sequences":
                        report.total_sequences = int(val)
                    elif key == "Sequences flagged as poor quality":
                        report.poor_quality = int(val)
                    elif key == "Sequence length":
                        report.sequence_length = val
                    elif key == "%GC":
                        report.gc_percent = int(val)
                    elif key == "Encoding":
                        report.encoding = val

        logger.info(f"Parsed {sample_name}: "
                    f"{'PASS' if report.passes_qc else 'FAIL'}")
        return report

    # ──────────────────────────────────────────────────────────
    # MULTIQC
    # ──────────────────────────────────────────────────────────

    def run_multiqc(self, search_dir: Path = None) -> Path:
        """
        Run MultiQC to combine all FastQC reports into one.

        Args:
            search_dir: Directory with FastQC reports.

        Returns:
            Path to the MultiQC HTML report
        """
        search_dir  = search_dir or self.output_dir
        multiqc_out = self.output_dir / "multiqc"
        multiqc_out.mkdir(exist_ok=True)

        logger.info(f"Running MultiQC on: {search_dir}")

        cmd = [
            MULTIQC_PATH,
            str(search_dir),
            "--outdir", str(multiqc_out),
            "--filename", "multiqc_report",
            "--force",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"MultiQC failed:\n{result.stderr}")
            raise RuntimeError("MultiQC failed")

        report_path = multiqc_out / "multiqc_report.html"
        logger.info(f"MultiQC report: {report_path}")
        return report_path

    # ──────────────────────────────────────────────────────────
    # REPORTING
    # ──────────────────────────────────────────────────────────

    def print_summary(self) -> None:
        """Print a summary table of all QC reports."""
        if not self.reports:
            print("No reports loaded.")
            return

        passed = sum(1 for r in self.reports if r.passes_qc)
        failed = len(self.reports) - passed

        print(f"\n{'='*70}")
        print(f"  QC SUMMARY  |  {len(self.reports)} samples  |  "
              f"Passed: {passed}  |  Failed: {failed}")
        print(f"{'='*70}")
        for report in self.reports:
            print(f"  {report.summary()}")
            for warning in report.warnings:
                print(f"      -> {warning}")
        print(f"{'='*70}\n")

    def get_failed_samples(self) -> list:
        """Return list of sample names that failed QC."""
        return [r.sample_name for r in self.reports if not r.passes_qc]

    def save_summary_csv(self, output_path: Path = None) -> Path:
        """Save QC summary to a CSV file."""
        output_path = output_path or self.output_dir / "qc_summary.csv"

        fields = [
            "sample_name", "total_sequences", "gc_percent",
            "sequence_length", "per_base_quality",
            "adapter_content", "passes_qc", "warnings"
        ]

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for report in self.reports:
                writer.writerow({
                    "sample_name"     : report.sample_name,
                    "total_sequences" : report.total_sequences,
                    "gc_percent"      : report.gc_percent,
                    "sequence_length" : report.sequence_length,
                    "per_base_quality": report.per_base_quality,
                    "adapter_content" : report.adapter_content,
                    "passes_qc"       : report.passes_qc,
                    "warnings"        : "; ".join(report.warnings),
                })

        logger.info(f"QC summary saved: {output_path}")
        return output_path