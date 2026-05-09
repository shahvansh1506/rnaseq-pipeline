"""
app.py — Interactive Streamlit dashboard for RNA-seq results.

Run with:
    streamlit run dashboard/app.py

What this dashboard shows:
- Summary metrics (total genes, up/down regulated)
- Interactive volcano plot
- PCA plot image
- Heatmap image
- Filterable results table
- Download button for results CSV
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import RESULTS_DIR

# ── Page configuration ────────────────────────────────────────
st.set_page_config(
    page_title = "RNA-seq Analysis Dashboard",
    page_icon  = "🧬",
    layout     = "wide",
)

# ── Custom CSS ────────────────────────────────────────────────
st.markdown("""
<style>
    .metric-card {
        background-color: #f0f2f6;
        border-radius: 10px;
        padding: 15px;
        text-align: center;
    }
    .stMetric { text-align: center; }
</style>
""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════

@st.cache_data
def load_results(path: str) -> pd.DataFrame:
    """Load DE results from CSV (cached for performance)."""
    return pd.read_csv(path)


def find_results_file() -> Path:
    """Find the DESeq2 results CSV."""
    candidates = [
        RESULTS_DIR / "deseq2" / "de_results.csv",
        Path("results/deseq2/de_results.csv"),
        Path("de_results.csv"),
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


# ══════════════════════════════════════════════════════════════
# MAIN DASHBOARD
# ══════════════════════════════════════════════════════════════

def main():

    # ── Header ────────────────────────────────────────────────
    st.title("🧬 RNA-seq Full Analysis Pipeline")
    st.markdown(
        "Interactive exploration of differential expression results. "
        "Use the sidebar to filter results."
    )
    st.divider()

    # ── Load data ─────────────────────────────────────────────
    results_path = find_results_file()

    if results_path is None:
        st.warning(
            "⚠️ No results found. Please run the pipeline first:\n\n"
            "```python demo_module5.py```"
        )
        st.stop()

    df = load_results(str(results_path))
    st.success(f"✅ Loaded results from: `{results_path}`")

    # ── Sidebar filters ───────────────────────────────────────
    st.sidebar.header("🔧 Filters")
    st.sidebar.markdown("Adjust thresholds to filter results")

    pval_thresh = st.sidebar.slider(
        "Adjusted P-value cutoff",
        min_value = 0.001,
        max_value = 0.1,
        value     = 0.05,
        step      = 0.001,
        format    = "%.3f"
    )

    lfc_thresh = st.sidebar.slider(
        "Log2 Fold Change cutoff",
        min_value = 0.5,
        max_value = 4.0,
        value     = 1.5,
        step      = 0.1,
    )

    # Filter by comparison if multiple exist
    if "comparison" in df.columns:
        comparisons = df["comparison"].unique().tolist()
        selected_comparison = st.sidebar.selectbox(
            "Comparison",
            options = comparisons,
        )
        df_filtered = df[df["comparison"] == selected_comparison].copy()
    else:
        df_filtered = df.copy()

    # ── Classify genes ────────────────────────────────────────
    df_filtered = df_filtered.dropna(subset=["padj", "log2FoldChange"])

    df_filtered["significance"] = "NS"
    df_filtered.loc[
        (df_filtered["padj"] < pval_thresh) &
        (df_filtered["log2FoldChange"] > lfc_thresh),
        "significance"
    ] = "Up"
    df_filtered.loc[
        (df_filtered["padj"] < pval_thresh) &
        (df_filtered["log2FoldChange"] < -lfc_thresh),
        "significance"
    ] = "Down"

    up   = df_filtered[df_filtered["significance"] == "Up"]
    down = df_filtered[df_filtered["significance"] == "Down"]
    sig  = pd.concat([up, down])

    # ── Summary metrics ───────────────────────────────────────
    st.subheader("📊 Summary")
    col1, col2, col3, col4 = st.columns(4)

    col1.metric(
        "Total Genes Tested",
        f"{len(df_filtered):,}"
    )
    col2.metric(
        "⬆️ Upregulated",
        f"{len(up):,}",
        delta = f"padj<{pval_thresh}, LFC>{lfc_thresh}"
    )
    col3.metric(
        "⬇️ Downregulated",
        f"{len(down):,}",
        delta = f"padj<{pval_thresh}, LFC<-{lfc_thresh}"
    )
    col4.metric(
        "Total Significant",
        f"{len(sig):,}",
        delta = f"{len(sig)/len(df_filtered)*100:.1f}% of genes"
    )

    st.divider()

    # ── Volcano plot ──────────────────────────────────────────
    st.subheader("🌋 Volcano Plot")

    df_filtered["-log10(padj)"] = -np.log10(
        df_filtered["padj"].clip(lower=1e-300)
    )

    # Label top 15 significant genes
    top_labels = sig.nsmallest(15, "padj") if len(sig) > 0 else pd.DataFrame()

    fig_volcano = px.scatter(
        df_filtered,
        x          = "log2FoldChange",
        y          = "-log10(padj)",
        color      = "significance",
        hover_name = "gene" if "gene" in df_filtered.columns else None,
        hover_data = {
            "log2FoldChange": ":.3f",
            "padj"          : ":.2e",
            "significance"  : False,
        },
        color_discrete_map = {
            "Up"  : "#E41A1C",
            "Down": "#377EB8",
            "NS"  : "#AAAAAA"
        },
        opacity = 0.6,
        labels  = {
            "log2FoldChange": "Log2 Fold Change",
            "-log10(padj)"  : "-Log10 Adjusted P-value",
        },
    )

    # Add threshold lines
    fig_volcano.add_vline(
        x=lfc_thresh,   line_dash="dash",
        line_color="black", opacity=0.4
    )
    fig_volcano.add_vline(
        x=-lfc_thresh,  line_dash="dash",
        line_color="black", opacity=0.4
    )
    fig_volcano.add_hline(
        y=-np.log10(pval_thresh), line_dash="dash",
        line_color="black", opacity=0.4
    )

    # Add gene labels for top hits
    if len(top_labels) > 0 and "gene" in top_labels.columns:
        fig_volcano.add_trace(go.Scatter(
            x    = top_labels["log2FoldChange"],
            y    = -np.log10(top_labels["padj"].clip(lower=1e-300)),
            mode = "text",
            text = top_labels["gene"],
            textposition = "top center",
            textfont     = dict(size=9, color="black"),
            showlegend   = False,
        ))

    fig_volcano.update_layout(
        height      = 500,
        showlegend  = True,
    )

    st.plotly_chart(fig_volcano, use_container_width=True)

    # ── Plots row ─────────────────────────────────────────────
    st.subheader("📈 Additional Plots")
    col_pca, col_heat = st.columns(2)

    deseq2_dir = RESULTS_DIR / "deseq2"

    with col_pca:
        pca_path = deseq2_dir / "pca_plot.png"
        if pca_path.exists():
            st.image(str(pca_path), caption="PCA Plot — Sample Clustering",
                     use_container_width=True)
        else:
            st.info("PCA plot not found. Run demo_module5.py first.")

    with col_heat:
        heat_path = deseq2_dir / "heatmap_top50.png"
        if heat_path.exists():
            st.image(str(heat_path),
                     caption="Heatmap — Top 50 DE Genes",
                     use_container_width=True)
        else:
            st.info("Heatmap not found. Run demo_module5.py first.")

    st.divider()

    # ── Results table ─────────────────────────────────────────
    st.subheader("📋 Significant Genes Table")

    if len(sig) > 0:
        display_cols = [c for c in
                        ["gene", "log2FoldChange", "padj",
                         "baseMean", "significance", "comparison"]
                        if c in sig.columns]

        display_df = sig[display_cols].sort_values("padj").reset_index(drop=True)
        display_df["log2FoldChange"] = display_df["log2FoldChange"].round(3)
        display_df["padj"]           = display_df["padj"].apply(
            lambda x: f"{x:.2e}"
        )

        st.dataframe(
            display_df,
            use_container_width = True,
            height              = 400,
        )
    else:
        st.info("No significant genes with current filter settings. "
                "Try relaxing the thresholds in the sidebar.")

    # ── Download button ───────────────────────────────────────
    st.divider()
    st.subheader("⬇️ Download Results")

    col_dl1, col_dl2 = st.columns(2)

    with col_dl1:
        st.download_button(
            label     = "📥 Download Significant Genes CSV",
            data      = sig.to_csv(index=False) if len(sig) > 0 else "",
            file_name = "significant_genes.csv",
            mime      = "text/csv",
            disabled  = len(sig) == 0,
        )

    with col_dl2:
        st.download_button(
            label     = "📥 Download All Results CSV",
            data      = df_filtered.to_csv(index=False),
            file_name = "all_de_results.csv",
            mime      = "text/csv",
        )

    # ── Footer ────────────────────────────────────────────────
    st.divider()
    st.markdown(
        "Built with 🧬 **RNA-seq Full Analysis Pipeline** | "
        "Python + R + DESeq2 + Streamlit"
    )


if __name__ == "__main__":
    main()
