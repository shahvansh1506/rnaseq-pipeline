# install_r_packages.R
# Run this script ONCE to install all R packages needed for DESeq2 analysis.
#
# Run from Anaconda Prompt:
#     Rscript install_r_packages.R

cat("Installing R packages for RNA-seq pipeline...\n")
cat("This may take 10-20 minutes the first time.\n\n")

# ── 1. Install BiocManager (gateway to Bioconductor packages) ─
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

# ── 2. Install Bioconductor packages ─────────────────────────
bioc_packages <- c(
  "DESeq2",          # Differential expression analysis
  "apeglm"           # Log fold change shrinkage for DESeq2
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(paste(pkg, "already installed\n"))
  }
}

# ── 3. Install CRAN packages ──────────────────────────────────
cran_packages <- c(
  "ggplot2",         # Plotting
  "pheatmap",        # Heatmaps
  "ggrepel",         # Non-overlapping labels on plots
  "dplyr",           # Data manipulation
  "RColorBrewer",    # Color palettes
  "tibble"           # Modern data frames
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, repos = "https://cran.r-project.org")
  } else {
    cat(paste(pkg, "already installed\n"))
  }
}

cat("\nAll packages installed successfully!\n")
cat("You can now run the DESeq2 analysis.\n")
