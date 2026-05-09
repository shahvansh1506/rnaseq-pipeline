"""
tests/test_qc.py — Unit tests for quality control and trimmer modules.

Run with:
    pytest tests/test_qc.py -v
"""

import gzip
import pytest
from pathlib import Path
from src.qc.quality_control import QualityController, QCReport
from src.qc.trimmer import Trimmer, TrimResult


# ══════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════

def make_fastq(path: Path, num_reads: int = 5,
               quality: str = "I") -> Path:
    """Create a test FASTQ.gz file with given quality character."""
    seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"  # 52bp
    with gzip.open(path, "wt") as f:
        for i in range(num_reads):
            f.write(f"@read_{i+1}\n")
            f.write(f"{seq}\n")
            f.write(f"+\n")
            f.write(f"{quality * len(seq)}\n")
    return path


FAKE_FASTQC_DATA = """\
##FastQC\t0.11.9
>>Basic Statistics\tPASS
#Measure\tValue
Filename\tsample.fastq.gz
File type\tConventional base calls
Encoding\tSanger / Illumina 1.9
Total Sequences\t1000000
Sequences flagged as poor quality\t0
Sequence length\t150
%GC\t48
>>END_MODULE
>>Per base sequence quality\tPASS
>>END_MODULE
>>Per sequence quality scores\tWARN
>>END_MODULE
>>Per sequence GC content\tPASS
>>END_MODULE
>>Adapter Content\tFAIL
>>END_MODULE
>>Sequence Duplication Levels\tPASS
>>END_MODULE
"""


# ══════════════════════════════════════════════════════════════
# QC REPORT TESTS
# ══════════════════════════════════════════════════════════════

class TestQCReport:

    def test_default_passes_qc(self):
        report = QCReport(sample_name="test")
        assert report.passes_qc is True

    def test_summary_contains_sample_name(self):
        report = QCReport(sample_name="my_sample", passes_qc=True)
        assert "my_sample" in report.summary()

    def test_summary_pass_status(self):
        report = QCReport(sample_name="s1", passes_qc=True)
        assert "PASS" in report.summary()

    def test_summary_fail_status(self):
        report = QCReport(sample_name="s1", passes_qc=False)
        assert "FAIL" in report.summary()


# ══════════════════════════════════════════════════════════════
# QUALITY CONTROLLER TESTS
# ══════════════════════════════════════════════════════════════

class TestQualityController:

    def test_parse_fastqc_data_basic_stats(self, tmp_path):
        qc = QualityController(output_dir=tmp_path)
        report = qc._parse_fastqc_data(FAKE_FASTQC_DATA, "sample")

        assert report.total_sequences == 1000000
        assert report.gc_percent      == 48
        assert report.sequence_length == "150"
        assert report.encoding        == "Sanger / Illumina 1.9"

    def test_parse_fastqc_data_grades(self, tmp_path):
        qc = QualityController(output_dir=tmp_path)
        report = qc._parse_fastqc_data(FAKE_FASTQC_DATA, "sample")

        assert report.basic_statistics     == "PASS"
        assert report.per_base_quality     == "PASS"
        assert report.per_sequence_quality == "WARN"
        assert report.adapter_content      == "FAIL"

    def test_parse_fastqc_data_fail_sets_passes_qc_false(self, tmp_path):
        qc = QualityController(output_dir=tmp_path)
        report = qc._parse_fastqc_data(FAKE_FASTQC_DATA, "sample")
        # FAIL in adapter content should mark sample as failed
        assert report.passes_qc is False

    def test_parse_fastqc_data_warnings_collected(self, tmp_path):
        qc = QualityController(output_dir=tmp_path)
        report = qc._parse_fastqc_data(FAKE_FASTQC_DATA, "sample")
        # Should have WARN for per_sequence_quality and FAIL for adapter
        assert len(report.warnings) == 2

    def test_get_failed_samples(self, tmp_path):
        qc = QualityController(output_dir=tmp_path)
        qc.reports = [
            QCReport(sample_name="good", passes_qc=True),
            QCReport(sample_name="bad",  passes_qc=False),
        ]
        failed = qc.get_failed_samples()
        assert failed == ["bad"]

    def test_save_summary_csv(self, tmp_path):
        qc = QualityController(output_dir=tmp_path)
        qc.reports = [
            QCReport(sample_name="s1", total_sequences=1000,
                     gc_percent=50, passes_qc=True),
        ]
        csv_path = qc.save_summary_csv(tmp_path / "summary.csv")
        assert csv_path.exists()

        content = csv_path.read_text()
        assert "s1" in content
        assert "1000" in content


# ══════════════════════════════════════════════════════════════
# TRIMMER TESTS
# ══════════════════════════════════════════════════════════════

class TestTrimmer:

    def test_quality_trim_3prime_removes_low_qual(self):
        from src.parsers.fastq_parser import FastqRecord
        trimmer = Trimmer(min_quality=20)

        # Last 3 bases have quality '!' = Phred 0 → should be trimmed
        record = FastqRecord(
            header   = "r1",
            sequence = "ATCGATCGATCG",
            plus     = "+",
            quality  = "IIIIIIIII!!!"   # last 3 are bad
        )
        trimmed = trimmer._quality_trim_3prime(record)
        assert len(trimmed.sequence) == 9
        assert trimmed.sequence == "ATCGATCGA"

    def test_quality_trim_3prime_keeps_good_read(self):
        from src.parsers.fastq_parser import FastqRecord
        trimmer = Trimmer(min_quality=20)

        record = FastqRecord("r1", "ATCGATCG", "+", "IIIIIIII")
        trimmed = trimmer._quality_trim_3prime(record)
        assert trimmed.sequence == "ATCGATCG"   # unchanged

    def test_trim_adapters_removes_adapter(self):
        from src.parsers.fastq_parser import FastqRecord
        trimmer = Trimmer()

        # Read contains start of TruSeq adapter
        adapter_start = "AGATCGGAAGAGCACACGTC"
        record = FastqRecord(
            header   = "r1",
            sequence = "ATCGATCG" + adapter_start,
            plus     = "+",
            quality  = "I" * (8 + len(adapter_start))
        )
        trimmed = trimmer._trim_adapters(record)
        assert trimmed.sequence == "ATCGATCG"

    def test_python_fallback_trims_correctly(self, tmp_path):
        # Create input with mixed quality reads
        input_f  = make_fastq(tmp_path / "input.fastq.gz",
                               num_reads=10, quality="I")
        output_f = tmp_path / "output.fastq.gz"

        trimmer = Trimmer(
            output_dir  = tmp_path,
            min_quality = 20,
            min_length  = 30,
        )
        result = trimmer.trim_python_fallback(input_f, output_f)

        assert result.input_reads == 10
        assert result.surviving_reads <= result.input_reads
        assert isinstance(result.survival_pct, float)
        assert output_f.exists()

    def test_trim_result_summary(self):
        result = TrimResult(
            sample_name     = "test",
            input_reads     = 1000,
            surviving_reads = 900,
            dropped_reads   = 100,
            survival_pct    = 90.0,
        )
        summary = result.summary()
        assert "900" in summary
        assert "1,000" in summary
        assert "90.0%" in summary
