# de_analysis.R
# ============================================================
# Differential Expression Analysis using DESeq2
#
# What this script does:
# 1. Loads count matrix and sample metadata
# 2. Builds a DESeq2 object
# 3. Runs differential expression analysis
# 4. Saves results table (CSV)
# 5. Generates volcano plot (PNG)
# 6. Generates PCA plot (PNG)
# 7. Generates heatmap of top 50 DE genes (PNG)
#
# Called from Python via run_deseq2.py:
#     Rscript de_analysis.R <counts> <metadata> <output_dir> <condition_col>
#
# Or run directly:
#     Rscript de_analysis.R data/counts/count_matrix.csv \
#                           data/counts/metadata.csv \
#                           results/ condition
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(ggrepel)
  library(RColorBrewer)
  library(tibble)
})

# ── 1. Parse command line arguments ──────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript de_analysis.R <counts_csv> <metadata_csv>",
      "<output_dir> <condition_column>\n")
  quit(status = 1)
}

counts_file   <- args[1]
metadata_file <- args[2]
output_dir    <- args[3]
condition_col <- args[4]

cat("============================================================\n")
cat("  DESeq2 Differential Expression Analysis\n")
cat("============================================================\n")
cat(paste("Counts file  :", counts_file, "\n"))
cat(paste("Metadata file:", metadata_file, "\n"))
cat(paste("Output dir   :", output_dir, "\n"))
cat(paste("Condition col:", condition_col, "\n\n"))

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── 2. Load data ──────────────────────────────────────────────
cat("Loading data...\n")

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
meta   <- read.csv(metadata_file, row.names = 1)

cat(paste("  Genes  :", nrow(counts), "\n"))
cat(paste("  Samples:", ncol(counts), "\n"))

# Make sure sample order matches between counts and metadata
# DESeq2 requires these to be in the same order
shared_samples <- intersect(colnames(counts), rownames(meta))

if (length(shared_samples) == 0) {
  stop("No matching samples between count matrix and metadata!")
}

counts <- counts[, shared_samples]
meta   <- meta[shared_samples, , drop = FALSE]

cat(paste("  Matched samples:", length(shared_samples), "\n\n"))

# ── 3. Build DESeq2 object ────────────────────────────────────
cat("Building DESeq2 dataset...\n")

# Convert condition column to factor
meta[[condition_col]] <- factor(meta[[condition_col]])

# Build the DESeq2 dataset
# design = which column defines the groups to compare
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),             # Must be integers
  colData   = meta,
  design    = as.formula(paste("~", condition_col))
)

# ── 4. Pre-filter low count genes ────────────────────────────
# Remove genes with fewer than 10 total counts
# This speeds up analysis without losing meaningful genes
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat(paste("  Genes after filtering:", nrow(dds), "\n\n"))

# ── 5. Run DESeq2 ─────────────────────────────────────────────
cat("Running DESeq2 (this may take a few minutes)...\n")
dds     <- DESeq(dds)
cat("DESeq2 complete!\n\n")

# Get names of all conditions for pairwise comparisons
conditions <- levels(meta[[condition_col]])
cat(paste("Conditions found:", paste(conditions, collapse = ", "), "\n\n"))

# ── 6. Extract results ────────────────────────────────────────
# Compare first two conditions as the main comparison
# e.g. condition2 vs condition1
cat("Extracting results...\n")

results_list <- list()

# Loop through all pairwise comparisons
for (i in 1:(length(conditions) - 1)) {
  for (j in (i + 1):length(conditions)) {

    cond_a <- conditions[i]
    cond_b <- conditions[j]
    comparison_name <- paste0(cond_b, "_vs_", cond_a)

    cat(paste("  Comparing:", cond_b, "vs", cond_a, "\n"))

    res <- results(
      dds,
      contrast = c(condition_col, cond_b, cond_a),
      alpha    = 0.05
    )

    # Convert to data frame and add gene column
    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene") %>%
      arrange(padj) %>%
      mutate(comparison = comparison_name)

    results_list[[comparison_name]] <- res_df
  }
}

# Save all results to one CSV
all_results <- bind_rows(results_list)
results_path <- file.path(output_dir, "de_results.csv")
write.csv(all_results, results_path, row.names = FALSE)
cat(paste("\nResults saved:", results_path, "\n"))

# Use first comparison for plots
main_results    <- results_list[[1]]
main_comparison <- names(results_list)[1]
cat(paste("Generating plots for:", main_comparison, "\n\n"))

# ── 7. Volcano plot ───────────────────────────────────────────
cat("Generating volcano plot...\n")

# Label significance categories
plot_data <- main_results %>%
  filter(!is.na(padj)) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange >  1.5 ~ "Up",
      padj < 0.05 & log2FoldChange < -1.5 ~ "Down",
      TRUE                                ~ "NS"
    ),
    neg_log10_padj = -log10(padj)
  )

# Label top 10 significant genes by fold change
top_genes <- plot_data %>%
  filter(significance != "NS") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10)

volcano_plot <- ggplot(
  plot_data,
  aes(x = log2FoldChange, y = neg_log10_padj, color = significance)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data    = top_genes,
    aes(label = gene),
    size    = 3,
    color   = "black",
    max.overlaps = 20
  ) +
  scale_color_manual(
    values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")
  ) +
  geom_vline(xintercept = c(-1.5, 1.5),
             linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title    = paste("Volcano Plot:", main_comparison),
    subtitle = paste(
      sum(plot_data$significance == "Up"),   "upregulated |",
      sum(plot_data$significance == "Down"), "downregulated"
    ),
    x        = "Log2 Fold Change",
    y        = "-Log10 Adjusted P-value",
    color    = "Expression"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

volcano_path <- file.path(output_dir, "volcano_plot.png")
ggsave(volcano_path, volcano_plot, width = 8, height = 6, dpi = 300)
cat(paste("  Saved:", volcano_path, "\n"))

# ── 8. PCA plot ───────────────────────────────────────────────
cat("Generating PCA plot...\n")

# Variance stabilizing transformation for visualization
# (NOT used for DE testing — only for PCA and heatmap)
vsd      <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = condition_col, returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, color = "black") +
  labs(
    title    = "PCA — Sample Clustering",
    subtitle = "Samples should cluster by condition",
    x        = paste0("PC1: ", pct_var[1], "% variance"),
    y        = paste0("PC2: ", pct_var[2], "% variance"),
    color    = "Condition"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

pca_path <- file.path(output_dir, "pca_plot.png")
ggsave(pca_path, pca_plot, width = 7, height = 6, dpi = 300)
cat(paste("  Saved:", pca_path, "\n"))

# ── 9. Heatmap of top 50 DE genes ────────────────────────────
cat("Generating heatmap...\n")

# Get top 50 most significantly DE genes
top50 <- main_results %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)

if (length(top50) >= 2) {

  # Get VST expression values for top genes
  mat <- assay(vsd)[top50, ]

  # Center each gene (subtract its mean across samples)
  # This makes the heatmap show relative expression
  mat <- mat - rowMeans(mat)

  # Sample annotation bar at top of heatmap
  annotation_col <- as.data.frame(
    colData(dds)[, condition_col, drop = FALSE]
  )

  heatmap_path <- file.path(output_dir, "heatmap_top50.png")
  pheatmap(
    mat,
    annotation_col  = annotation_col,
    show_rownames   = TRUE,
    show_colnames   = TRUE,
    scale           = "row",
    color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    fontsize_row    = 8,
    fontsize_col    = 10,
    filename        = heatmap_path,
    width           = 10,
    height          = 12
  )
  cat(paste("  Saved:", heatmap_path, "\n"))

} else {
  cat("  Not enough significant genes for heatmap (need >= 2)\n")
}

# ── 10. Summary ───────────────────────────────────────────────
cat("\n============================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("============================================================\n")

for (comp_name in names(results_list)) {
  res      <- results_list[[comp_name]]
  n_sig    <- sum(!is.na(res$padj) & res$padj < 0.05, na.rm = TRUE)
  n_up     <- sum(!is.na(res$padj) & res$padj < 0.05 &
                    res$log2FoldChange > 1.5, na.rm = TRUE)
  n_down   <- sum(!is.na(res$padj) & res$padj < 0.05 &
                    res$log2FoldChange < -1.5, na.rm = TRUE)

  cat(paste0("  ", comp_name, "\n"))
  cat(paste0("    Significant genes : ", n_sig, "\n"))
  cat(paste0("    Upregulated       : ", n_up, "\n"))
  cat(paste0("    Downregulated     : ", n_down, "\n\n"))
}

cat(paste("Output files saved to:", output_dir, "\n"))
cat("============================================================\n")
