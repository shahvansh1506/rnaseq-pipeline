# RNA-seq Analysis Pipeline

A full end-to-end RNA-seq pipeline built with Python and R — takes raw sequencing reads and turns them into differential expression results with interactive visualizations.

Built this because most RNA-seq tutorials stop halfway. This one doesn't.

---

## What it does

Raw sequencing data is messy — gigabytes of compressed reads, quality issues, alignment headaches, and statistical complexity. This pipeline handles all of it:

1. **Downloads** public RNA-seq datasets directly from NCBI GEO
2. **Parses** raw FASTQ files without loading them fully into memory
3. **Quality checks** reads with FastQC and flags problematic samples
4. **Trims** low-quality bases and adapter sequences
5. **Aligns** reads to a reference genome using HISAT2
6. **Counts** reads per gene using featureCounts
7. **Finds DE genes** using DESeq2 in R with proper statistical modelling
8. **Validates** results against the NCBI Gene database
9. **Visualizes** everything in an interactive Streamlit dashboard

---

## The biology behind it

RNA-seq measures gene expression — which genes are switched on in a cell and by how much. By comparing expression across conditions (e.g. healthy vs diseased, treated vs untreated), we can identify which genes are driving biological differences.

The dataset used for development is **GSE60450** — a mammary gland study comparing basal and luminal cell populations across virgin, pregnant and lactating mice. Classic benchmark dataset in the bioinformatics community.

---

## Tech stack

| Layer | Tools |
|---|---|
| Language | Python 3.11, R 4.3 |
| DE Analysis | DESeq2, edgeR |
| Alignment | HISAT2, featureCounts |
| Quality Control | FastQC, MultiQC, Trimmomatic |
| Data | Biopython, GEOparse, pandas, numpy |
| Dashboard | Streamlit, Plotly |
| Database | NCBI Entrez API, SQLAlchemy |
| Testing | pytest (65 tests) |

---

## Project structure

```
rnaseq_pipeline/
├── src/
│   ├── config.py               # all paths and parameters in one place
│   ├── parsers/
│   │   ├── geo_fetcher.py      # downloads datasets from NCBI GEO
│   │   ├── fastq_parser.py     # memory-efficient FASTQ reader
│   │   └── count_parser.py     # loads and cleans count matrices
│   ├── qc/
│   │   ├── quality_control.py  # FastQC wrapper + report parser
│   │   └── trimmer.py          # Trimmomatic wrapper + Python fallback
│   ├── alignment/
│   │   ├── aligner.py          # HISAT2 wrapper
│   │   └── quantifier.py       # featureCounts wrapper
│   ├── analysis/
│   │   ├── de_analysis.R       # DESeq2 script (volcano, heatmap, PCA)
│   │   └── run_deseq2.py       # Python wrapper that calls the R script
│   └── validation/
│       └── ncbi_validator.py   # checks DE genes against NCBI
├── dashboard/
│   └── app.py                  # Streamlit dashboard
├── tests/                      # 65 unit tests across all modules
├── environment.yml             # conda environment
└── requirements.txt
```

---

## Setup

### Requirements
- Anaconda / Miniconda
- R 4.3+
- Git

### Installation

```bash
# clone the repo
git clone https://github.com/shahvansh1506/rnaseq-pipeline.git
cd rnaseq-pipeline

# create and activate environment
conda env create -f environment.yml
conda activate rnaseq_pipeline

# install R packages (takes ~15 mins first time)
Rscript install_r_packages.R
```

Then open `src/config.py` and set your email:
```python
ENTREZ_EMAIL = "your@email.com"
```

---

## Running the pipeline

```bash
# test each module works
python demo_module2.py    # FASTQ parser
python demo_module3.py    # QC and trimming
python demo_module4.py    # alignment and quantification
python demo_module5.py    # DESeq2 differential expression

# launch the dashboard
streamlit run dashboard/app.py

# run all tests
pytest tests/ -v
```

---

## Key design decisions

**Generator-based FASTQ parser** — standard file reading loads the whole file into RAM. With FASTQ files hitting 50GB+, that's not viable. The parser yields one read at a time, keeping memory usage flat regardless of file size.

**Python wrapper around R** — DESeq2 is genuinely the best tool for RNA-seq differential expression, but it lives in R. Rather than rewriting it, the pipeline calls it via subprocess and handles all the file I/O from Python. Best of both worlds.

**Python fallback trimmer** — Trimmomatic does not run on Windows. Rather than blocking Windows users, the pipeline detects this and falls back to a pure Python trimmer that does the same job.

**Modular architecture** — each step reads from and writes to files. This means you can re-run any single step without rerunning the whole thing, and swap out any tool (e.g. STAR instead of HISAT2) with minimal changes.

---

## Test results

```
65 passed in 12.3s
```

Tests cover: FASTQ parsing, GC content calculation, quality trimming,
adapter removal, count matrix filtering, CPM normalization, metadata
validation, DESeq2 result parsing, NCBI validation, and dashboard logic.

---

## Output files

After running the full pipeline:

```
results/deseq2/
├── de_results.csv       # full table of all genes with stats
├── volcano_plot.png     # log2FC vs significance
├── pca_plot.png         # sample clustering
└── heatmap_top50.png    # top 50 DE genes across samples
```

---

## What I learned building this

- How RNA-seq data actually flows from sequencer to biological insight
- Why generators matter when files are larger than your RAM
- How DESeq2's negative binomial model handles biological replication
- How to wrap R scripts cleanly in a Python pipeline
- How to write tests for bioinformatics code that depends on external tools

---

## Roadmap

- [ ] Snakemake workflow for full automation
- [ ] STAR aligner support
- [ ] Gene ontology enrichment (clusterProfiler)
- [ ] Multi-comparison support in the dashboard
- [ ] Docker container for reproducibility

---

## License

MIT
