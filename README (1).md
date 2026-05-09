# 🧬 RNA-seq Full Analysis Pipeline

A production-grade, end-to-end RNA-seq analysis pipeline built with Python and R.

## What this pipeline does

1. Downloads RNA-seq datasets from NCBI GEO
2. Parses raw FASTQ files memory-efficiently using generators
3. Quality controls raw reads (FastQC + MultiQC)
4. Trims low-quality bases and adapters (Trimmomatic / Python fallback)
5. Aligns reads to a reference genome (HISAT2)
6. Quantifies gene expression (featureCounts)
7. Analyses differential expression (DESeq2 in R)
8. Visualizes results (volcano plot, heatmap, PCA plot)
9. Validates DE genes against NCBI Gene database
10. Serves results in an interactive Streamlit dashboard

## Setup

### 1. Create conda environment
```bash
conda env create -f environment.yml
conda activate rnaseq_pipeline
```

### 2. Install R packages
```bash
Rscript install_r_packages.R
```

### 3. Update config.py
```python
ENTREZ_EMAIL = "your_email@example.com"
THREADS      = 8
```

## Running the pipeline

```bash
python demo_module2.py    # FASTQ parser
python demo_module3.py    # QC + trimming
python demo_module4.py    # Alignment + counting
python demo_module5.py    # DESeq2 analysis
streamlit run dashboard/app.py   # Launch dashboard
pytest tests/ -v          # Run all tests
```
