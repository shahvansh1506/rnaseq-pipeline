"""
config.py — Global configuration for the RNA-seq pipeline.

All paths, parameters, and tool settings live here.
Import this module at the top of any other module:
    from src.config import RAW_DIR, THREADS, PVALUE_CUTOFF
"""

from pathlib import Path

# ═══════════════════════════════════════════════════════════════
# DIRECTORY STRUCTURE
# ═══════════════════════════════════════════════════════════════

BASE_DIR      = Path(__file__).resolve().parent.parent

# Input/output data directories
DATA_DIR      = BASE_DIR / "data"
RAW_DIR       = DATA_DIR / "raw"           # Raw FASTQ files from GEO
PROCESSED_DIR = DATA_DIR / "processed"     # Trimmed reads
ALIGNED_DIR   = DATA_DIR / "aligned"       # BAM files from HISAT2
COUNTS_DIR    = DATA_DIR / "counts"        # Count matrices
REF_DIR       = DATA_DIR / "reference"     # Genome FASTA + GTF

# Output directories
RESULTS_DIR   = BASE_DIR / "results"       # DE results, plots
REPORTS_DIR   = BASE_DIR / "reports"       # R Markdown reports

# Auto-create all directories on import
for _dir in [RAW_DIR, PROCESSED_DIR, ALIGNED_DIR, COUNTS_DIR,
             REF_DIR, RESULTS_DIR, REPORTS_DIR]:
    _dir.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════
# PIPELINE PARAMETERS
# ═══════════════════════════════════════════════════════════════

# Compute resources
THREADS = 8                  # CPU threads for alignment & counting

# Read quality filtering (Trimmomatic)
MIN_QUALITY = 20             # Minimum Phred quality score to keep
MIN_LENGTH  = 50             # Minimum read length after trimming (bp)

# Differential expression thresholds (DESeq2)
PVALUE_CUTOFF  = 0.05        # Adjusted p-value significance cutoff
LOG2FC_CUTOFF  = 1.5         # Minimum absolute log2 fold change

# Heatmap: number of top DE genes to display
TOP_GENES_HEATMAP = 50

# ═══════════════════════════════════════════════════════════════
# EXTERNAL TOOL PATHS
# ═══════════════════════════════════════════════════════════════
# These assume tools are installed in your conda environment.
# Change these if you have custom installation paths.

HISAT2_PATH          = "hisat2"
HISAT2_BUILD_PATH    = "hisat2-build"
FASTQC_PATH          = "fastqc"
TRIMMOMATIC_PATH     = "trimmomatic"
FEATURECOUNTS_PATH   = "featureCounts"
MULTIQC_PATH         = "multiqc"
RSCRIPT_PATH         = "Rscript"
SAMTOOLS_PATH        = "samtools"

# ═══════════════════════════════════════════════════════════════
# DATABASE
# ═══════════════════════════════════════════════════════════════

DB_PATH = BASE_DIR / "pipeline.db"     # SQLite database for metadata

# ═══════════════════════════════════════════════════════════════
# NCBI / ENTREZ SETTINGS
# ═══════════════════════════════════════════════════════════════
# Required by Biopython's Entrez module for gene validation.
# Replace with your own email address.

ENTREZ_EMAIL = "shahvansh1506@gmail.com"

# ═══════════════════════════════════════════════════════════════
# LOGGING
# ═══════════════════════════════════════════════════════════════

LOG_LEVEL  = "INFO"          # Options: DEBUG, INFO, WARNING, ERROR
LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
