"""
tests/test_parsers.py — Unit tests for geo_fetcher and fastq_parser.

Run with:
    pytest tests/test_parsers.py -v
"""

import pytest
import gzip
import tempfile
from pathlib import Path
from src.parsers.fastq_parser import (
    parse_fastq,
    get_fastq_stats,
    FastqRecord,
)


# ══════════════════════════════════════════════════════════════
# HELPERS — create tiny test FASTQ files
# ══════════════════════════════════════════════════════════════

def make_fastq(path: Path, compressed: bool = False) -> Path:
    """Write a minimal 2-read FASTQ file for testing."""
    content = (
        "@read1 test\n"
        "ATCGATCG\n"
        "+\n"
        "IIIIIIII\n"   # Phred 40 — excellent quality
        "@read2 test\n"
        "GCGCGCGC\n"
        "+\n"
        "!!!!!!!!\n"   # Phred 0  — terrible quality
    )
    opener = gzip.open if compressed else open
    with opener(path, "wt") as f:
        f.write(content)
    return path


# ══════════════════════════════════════════════════════════════
# FASTQ RECORD TESTS
# ══════════════════════════════════════════════════════════════

class TestFastqRecord:

    def test_length(self):
        record = FastqRecord("r1", "ATCGATCG", "+", "IIIIIIII")
        assert record.length == 8

    def test_gc_content_all_gc(self):
        record = FastqRecord("r1", "GCGCGCGC", "+", "IIIIIIII")
        assert record.gc_content == 100.0

    def test_gc_content_no_gc(self):
        record = FastqRecord("r1", "ATATATATAT", "+", "IIIIIIIII!")
        assert record.gc_content == 0.0

    def test_gc_content_mixed(self):
        record = FastqRecord("r1", "ATGC", "+", "IIII")
        assert record.gc_content == 50.0

    def test_mean_quality_high(self):
        # 'I' = ASCII 73, Phred = 73 - 33 = 40
        record = FastqRecord("r1", "ATCG", "+", "IIII")
        assert record.mean_quality == 40.0

    def test_mean_quality_low(self):
        # '!' = ASCII 33, Phred = 33 - 33 = 0
        record = FastqRecord("r1", "ATCG", "+", "!!!!")
        assert record.mean_quality == 0.0

    def test_is_good_quality_passes(self):
        record = FastqRecord("r1", "A" * 60, "+", "I" * 60)
        assert record.is_good_quality(min_quality=20, min_length=50) is True

    def test_is_good_quality_too_short(self):
        record = FastqRecord("r1", "ATCG", "+", "IIII")
        assert record.is_good_quality(min_quality=20, min_length=50) is False

    def test_is_good_quality_low_qual(self):
        record = FastqRecord("r1", "A" * 60, "+", "!" * 60)
        assert record.is_good_quality(min_quality=20, min_length=50) is False


# ══════════════════════════════════════════════════════════════
# PARSER TESTS
# ══════════════════════════════════════════════════════════════

class TestParseFastq:

    def test_parse_plain_fastq(self, tmp_path):
        path = make_fastq(tmp_path / "test.fastq")
        records = list(parse_fastq(path))
        assert len(records) == 2

    def test_parse_gzipped_fastq(self, tmp_path):
        path = make_fastq(tmp_path / "test.fastq.gz", compressed=True)
        records = list(parse_fastq(path))
        assert len(records) == 2

    def test_record_fields(self, tmp_path):
        path = make_fastq(tmp_path / "test.fastq")
        records = list(parse_fastq(path))
        assert records[0].header   == "read1 test"
        assert records[0].sequence == "ATCGATCG"
        assert records[0].quality  == "IIIIIIII"

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            list(parse_fastq(Path("nonexistent.fastq")))


# ══════════════════════════════════════════════════════════════
# STATS TESTS
# ══════════════════════════════════════════════════════════════

class TestFastqStats:

    def test_total_reads(self, tmp_path):
        path = make_fastq(tmp_path / "test.fastq")
        stats = get_fastq_stats(path)
        assert stats["total_reads"] == 2

    def test_total_bases(self, tmp_path):
        path = make_fastq(tmp_path / "test.fastq")
        stats = get_fastq_stats(path)
        assert stats["total_bases"] == 16   # 8 bases x 2 reads

    def test_gc_content(self, tmp_path):
        path = make_fastq(tmp_path / "test.fastq")
        stats = get_fastq_stats(path)
        # read1: ATCGATCG = 50% GC, read2: GCGCGCGC = 100% GC → 75% overall
        assert stats["gc_content_pct"] == 75.0

    def test_low_quality_detected(self, tmp_path):
        path = make_fastq(tmp_path / "test.fastq")
        # read2 has quality '!' = Phred 0 → should be flagged
        stats = get_fastq_stats(path, min_quality=20)
        assert stats["low_quality_reads"] == 1
        assert stats["low_quality_pct"]   == 50.0
