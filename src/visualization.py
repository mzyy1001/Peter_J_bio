"""
Visualization module: generates publication-quality figures for the paper.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import networkx as nx
from pathlib import Path
from itertools import combinations

FIGURES_DIR = Path(__file__).parent.parent / "results" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Color scheme for hallmark categories
CATEGORY_COLORS = {
    "Primary": "#E74C3C",
    "Antagonistic": "#F39C12",
    "Integrative": "#3498DB",
}

HALLMARK_COLORS = {
    "Genomic Instability": "#E74C3C",
    "Telomere Attrition": "#C0392B",
    "Epigenetic Alterations": "#E67E22",
    "Loss of Proteostasis": "#D35400",
    "Deregulated Nutrient Sensing": "#F1C40F",
    "Mitochondrial Dysfunction": "#E67E22",
    "Cellular Senescence": "#F39C12",
    "Stem Cell Exhaustion": "#3498DB",
    "Altered Intercellular Communication": "#2980B9",
}

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "sans-serif",
})


def plot_hallmark_correlation_heatmap(corr_matrix, p_matrix=None, title="",
                                      filename="hallmark_correlation_heatmap.pdf"):
    """Plot correlation heatmap with significance stars."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Short labels
    from .hallmark_genes import HALLMARK_SHORT
    short = [HALLMARK_SHORT.get(c, c) for c in corr_matrix.columns]

    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)

    sns.heatmap(corr_matrix.values, mask=mask, annot=True, fmt=".2f",
                cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                xticklabels=short, yticklabels=short,
                square=True, linewidths=0.5, ax=ax,
                cbar_kws={"label": "Spearman Correlation"})

    # Add significance stars
    if p_matrix is not None:
        for i in range(len(corr_matrix)):
            for j in range(i):
                p = p_matrix.iloc[i, j]
                if p < 0.001:
                    star = "***"
                elif p < 0.01:
                    star = "**"
                elif p < 0.05:
                    star = "*"
                else:
                    continue
                ax.text(j + 0.5, i + 0.75, star, ha="center", va="center",
                        fontsize=8, color="black")

    ax.set_title(title or "Hallmark Activity Score Correlations", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_network_graph(G, title="", filename="hallmark_network.pdf",
                       node_metrics=None, layout="spring"):
    """Plot hallmark interaction network."""
    fig, ax = plt.subplots(figsize=(12, 10))

    if layout == "spring":
        pos = nx.spring_layout(G, k=2.5, iterations=100, seed=42, weight="weight")
    elif layout == "circular":
        pos = nx.circular_layout(G)
    else:
        pos = nx.kamada_kawai_layout(G, weight="weight")

    # Node sizes based on degree or eigenvector centrality
    if node_metrics is not None and "eigenvector" in node_metrics.columns:
        node_sizes = [1500 + 3000 * node_metrics.loc[n, "eigenvector"]
                      for n in G.nodes()]
    else:
        node_sizes = [2000 for _ in G.nodes()]

    # Node colors by category
    from .hallmark_genes import HALLMARK_CATEGORIES
    cat_map = {}
    for cat, hmarks in HALLMARK_CATEGORIES.items():
        for h in hmarks:
            cat_map[h] = cat
    node_colors = [CATEGORY_COLORS.get(cat_map.get(n, ""), "#95A5A6")
                   for n in G.nodes()]

    # Edge widths and colors
    if G.number_of_edges() > 0:
        weights = [G[u][v].get("weight", 0.5) for u, v in G.edges()]
        max_w = max(weights) if weights else 1
        edge_widths = [1 + 4 * w / max_w for w in weights]

        edge_colors = []
        for u, v in G.edges():
            corr = G[u][v].get("correlation", G[u][v].get("partial_corr", 0))
            edge_colors.append("#E74C3C" if corr > 0 else "#3498DB")
    else:
        edge_widths = []
        edge_colors = []

    # Draw
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths,
                           edge_color=edge_colors, alpha=0.6)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                           node_color=node_colors, edgecolors="black",
                           linewidths=1.5, alpha=0.9)

    # Labels
    from .hallmark_genes import HALLMARK_SHORT
    label_map = {n: HALLMARK_SHORT.get(n, n) for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, label_map, ax=ax, font_size=10,
                           font_weight="bold")

    # Edge weight labels
    edge_labels = {}
    for u, v in G.edges():
        w = G[u][v].get("correlation", G[u][v].get("partial_corr",
              G[u][v].get("weight", 0)))
        edge_labels[(u, v)] = f"{w:.2f}"
    nx.draw_networkx_edge_labels(G, pos, edge_labels, ax=ax, font_size=7)

    # Legend
    legend_patches = [mpatches.Patch(color=c, label=cat)
                      for cat, c in CATEGORY_COLORS.items()]
    ax.legend(handles=legend_patches, loc="upper left", fontsize=9)

    ax.set_title(title or "Hallmark Interaction Network", fontweight="bold", fontsize=14)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_age_correlation_barplot(hallmark_scores, metadata,
                                 filename="age_hallmark_correlations.pdf"):
    """Bar plot of age-hallmark correlations."""
    from scipy import stats as sp_stats
    from .hallmark_genes import HALLMARK_SHORT, HALLMARK_CATEGORIES

    ages = metadata.set_index("sample_id").loc[hallmark_scores.index, "age"]
    hallmarks = hallmark_scores.columns.tolist()

    corrs, pvals = [], []
    for h in hallmarks:
        r, p = sp_stats.spearmanr(ages, hallmark_scores[h])
        corrs.append(r)
        pvals.append(p)

    df = pd.DataFrame({"hallmark": hallmarks, "correlation": corrs, "pvalue": pvals})
    df["short"] = df["hallmark"].map(HALLMARK_SHORT)
    df["significant"] = df["pvalue"] < 0.05

    cat_map = {}
    for cat, hs in HALLMARK_CATEGORIES.items():
        for h in hs:
            cat_map[h] = cat
    df["category"] = df["hallmark"].map(cat_map)
    df["color"] = df["category"].map(CATEGORY_COLORS)
    df = df.sort_values("correlation", ascending=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.barh(df["short"], df["correlation"], color=df["color"],
                   edgecolor="black", linewidth=0.5)

    # Add significance markers
    for i, (_, row) in enumerate(df.iterrows()):
        if row["pvalue"] < 0.001:
            ax.text(row["correlation"] + 0.02, i, "***", va="center", fontsize=9)
        elif row["pvalue"] < 0.01:
            ax.text(row["correlation"] + 0.02, i, "**", va="center", fontsize=9)
        elif row["pvalue"] < 0.05:
            ax.text(row["correlation"] + 0.02, i, "*", va="center", fontsize=9)

    ax.set_xlabel("Spearman Correlation with Age")
    ax.set_title("Age-Hallmark Activity Correlations", fontweight="bold")
    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")

    legend_patches = [mpatches.Patch(color=c, label=cat)
                      for cat, c in CATEGORY_COLORS.items()]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=9)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_feature_importance(importances_df, title="", filename="feature_importance.pdf"):
    """Plot feature importance from ML models."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    from .hallmark_genes import HALLMARK_SHORT

    for ax, col, name in zip(axes, ["rf_importance", "gb_importance"],
                              ["Random Forest", "Gradient Boosting"]):
        df = importances_df.sort_values(col, ascending=True)
        df["short"] = df["hallmark"].map(HALLMARK_SHORT)
        ax.barh(df["short"], df[col], color="#2ECC71", edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Feature Importance")
        ax.set_title(name, fontweight="bold")

    fig.suptitle(title or "Hallmark Importance for Age Group Prediction",
                 fontweight="bold", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_pca_scatter(pca_results, filename="pca_scatter.pdf"):
    """Plot PCA scatter colored by age and tissue."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    scores = pca_results["scores"]
    var_ratio = pca_results["explained_variance_ratio"]

    # Color by age
    sc1 = axes[0].scatter(scores["PC1"], scores["PC2"], c=scores["age"],
                          cmap="RdYlBu_r", s=30, alpha=0.7, edgecolors="black",
                          linewidth=0.3)
    axes[0].set_xlabel(f"PC1 ({var_ratio[0]*100:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({var_ratio[1]*100:.1f}%)")
    axes[0].set_title("Colored by Age", fontweight="bold")
    plt.colorbar(sc1, ax=axes[0], label="Age")

    # Color by tissue
    tissues = scores["tissue"].unique()
    palette = sns.color_palette("Set2", len(tissues))
    for i, t in enumerate(sorted(tissues)):
        mask = scores["tissue"] == t
        axes[1].scatter(scores.loc[mask, "PC1"], scores.loc[mask, "PC2"],
                       c=[palette[i]], label=t, s=30, alpha=0.7,
                       edgecolors="black", linewidth=0.3)
    axes[1].set_xlabel(f"PC1 ({var_ratio[0]*100:.1f}%)")
    axes[1].set_ylabel(f"PC2 ({var_ratio[1]*100:.1f}%)")
    axes[1].set_title("Colored by Tissue", fontweight="bold")
    axes[1].legend(fontsize=8)

    fig.suptitle("PCA of Hallmark Activity Scores", fontweight="bold", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_biological_age(age_results, filename="biological_age.pdf"):
    """Plot predicted vs actual age with age acceleration."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    actual = age_results["actual_ages"]
    predicted = age_results["predicted_ages"]
    accel = age_results["age_acceleration"]

    # Predicted vs actual
    axes[0].scatter(actual, predicted, c=accel["age_acceleration"],
                    cmap="RdBu_r", s=25, alpha=0.7, edgecolors="black", linewidth=0.3)
    lims = [min(actual.min(), predicted.min()) - 5,
            max(actual.max(), predicted.max()) + 5]
    axes[0].plot(lims, lims, "k--", linewidth=1, alpha=0.5)
    axes[0].set_xlabel("Chronological Age")
    axes[0].set_ylabel("Predicted Biological Age")
    axes[0].set_title("Biological Age Prediction", fontweight="bold")
    r2 = age_results["gb_r2_cv"]
    mae = age_results["gb_mae"]
    axes[0].text(0.05, 0.95, f"R² = {r2:.3f}\nMAE = {mae:.1f} yr",
                 transform=axes[0].transAxes, va="top", fontsize=10,
                 bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    # Age acceleration distribution
    axes[1].hist(accel["age_acceleration"], bins=25, color="#3498DB",
                 edgecolor="black", alpha=0.7)
    axes[1].axvline(0, color="red", linestyle="--", linewidth=1.5)
    axes[1].set_xlabel("Age Acceleration (years)")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Biological Age Acceleration Distribution", fontweight="bold")
    axes[1].text(0.05, 0.95,
                 f"Mean = {accel['age_acceleration'].mean():.2f}\n"
                 f"SD = {accel['age_acceleration'].std():.2f}",
                 transform=axes[1].transAxes, va="top", fontsize=10,
                 bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_cross_hallmark_prediction(cross_pred_df, filename="cross_hallmark_prediction.pdf"):
    """Plot cross-hallmark prediction R² values."""
    from .hallmark_genes import HALLMARK_SHORT

    fig, ax = plt.subplots(figsize=(8, 6))
    df = cross_pred_df.sort_values("rf_r2_cv", ascending=True).copy()
    df["short"] = df["target_hallmark"].map(HALLMARK_SHORT)

    colors = ["#E74C3C" if r > 0.5 else "#F39C12" if r > 0.3 else "#3498DB"
              for r in df["rf_r2_cv"]]
    ax.barh(df["short"], df["rf_r2_cv"], color=colors,
            edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Cross-validated R² (RF)")
    ax.set_title("Cross-Hallmark Predictability\n(each hallmark predicted from all others)",
                 fontweight="bold")
    ax.axvline(0.5, color="red", linestyle="--", alpha=0.5, label="R²=0.5")
    ax.legend()

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_interaction_changes(interaction_df, filename="interaction_age_changes.pdf"):
    """Plot how hallmark interactions change with age."""
    from .hallmark_genes import HALLMARK_SHORT

    # Use FDR-corrected significance if available, otherwise fall back to uncorrected
    sig_col = "significant_fdr" if "significant_fdr" in interaction_df.columns else "significant"
    sig = interaction_df[interaction_df[sig_col]].copy()
    if len(sig) == 0:
        sig = interaction_df.nlargest(10, "corr_change").copy()

    sig["label"] = sig.apply(
        lambda r: f"{HALLMARK_SHORT.get(r['hallmark_1'], r['hallmark_1'])}-"
                  f"{HALLMARK_SHORT.get(r['hallmark_2'], r['hallmark_2'])}",
        axis=1
    )
    sig = sig.sort_values("corr_change")

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ["#E74C3C" if c > 0 else "#3498DB" for c in sig["corr_change"]]
    ax.barh(sig["label"], sig["corr_change"], color=colors,
            edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Correlation Change (Old - Young)")
    ax.set_title("Age-Dependent Changes in Hallmark Interactions", fontweight="bold")
    ax.axvline(0, color="black", linewidth=0.5)

    for i, (_, row) in enumerate(sig.iterrows()):
        # Use FDR-corrected p-values for significance stars if available
        p_val = row["pvalue_fdr"] if "pvalue_fdr" in row.index else row["p_value"]
        stars = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        if stars:
            x_pos = row["corr_change"] + (0.01 if row["corr_change"] > 0 else -0.01)
            ax.text(x_pos, i, stars, va="center", fontsize=9,
                    ha="left" if row["corr_change"] > 0 else "right")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_partial_correlation_network(G, pcorr_df, filename="partial_corr_network.pdf"):
    """Plot partial correlation network (direct relationships only)."""
    plot_network_graph(G, title="Partial Correlation Network\n(Direct Hallmark Relationships)",
                       filename=filename, layout="kamada_kawai")


def plot_variance_explained(pca_results, filename="variance_explained.pdf"):
    """Plot PCA variance explained."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    var = pca_results["explained_variance_ratio"]
    cum_var = pca_results["cumulative_variance"]

    axes[0].bar(range(1, len(var)+1), var, color="#3498DB",
                edgecolor="black", linewidth=0.5)
    axes[0].set_xlabel("Principal Component")
    axes[0].set_ylabel("Variance Explained")
    axes[0].set_title("Scree Plot", fontweight="bold")

    axes[1].plot(range(1, len(cum_var)+1), cum_var, "o-", color="#E74C3C")
    axes[1].axhline(0.9, color="gray", linestyle="--", alpha=0.5)
    axes[1].set_xlabel("Number of Components")
    axes[1].set_ylabel("Cumulative Variance")
    axes[1].set_title("Cumulative Variance Explained", fontweight="bold")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_causal_graph(G, filename="causal_network.pdf"):
    """Plot the estimated causal/directed network."""
    from .hallmark_genes import HALLMARK_SHORT, HALLMARK_CATEGORIES

    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(G, k=3, iterations=100, seed=42)

    cat_map = {}
    for cat, hs in HALLMARK_CATEGORIES.items():
        for h in hs:
            cat_map[h] = cat

    node_colors = [CATEGORY_COLORS.get(cat_map.get(n, ""), "#95A5A6") for n in G.nodes()]
    label_map = {n: HALLMARK_SHORT.get(n, n) for n in G.nodes()}

    # Separate directed and undirected edges
    directed = [(u, v) for u, v, d in G.edges(data=True) if d.get("type") == "directed"]
    undirected = [(u, v) for u, v, d in G.edges(data=True) if d.get("type") == "undirected"]

    nx.draw_networkx_edges(G, pos, edgelist=directed, ax=ax,
                           edge_color="#E74C3C", width=2, alpha=0.7,
                           arrows=True, arrowsize=20,
                           connectionstyle="arc3,rad=0.1")
    nx.draw_networkx_edges(G, pos, edgelist=undirected, ax=ax,
                           edge_color="#95A5A6", width=1.5, alpha=0.5,
                           style="dashed")

    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=2500,
                           node_color=node_colors, edgecolors="black",
                           linewidths=1.5, alpha=0.9)
    nx.draw_networkx_labels(G, pos, label_map, ax=ax, font_size=10,
                           font_weight="bold")

    legend_patches = [
        mpatches.Patch(color="#E74C3C", label="Directed (causal)"),
        mpatches.Patch(color="#95A5A6", label="Undirected"),
    ] + [mpatches.Patch(color=c, label=cat) for cat, c in CATEGORY_COLORS.items()]
    ax.legend(handles=legend_patches, loc="upper left", fontsize=9)

    ax.set_title("Estimated Causal Structure (PC Algorithm)", fontweight="bold", fontsize=14)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_scoring_method_comparison(corr_matrices, method_names,
                                   filename="fig14_scoring_comparison.pdf"):
    """Plot correlation heatmaps for different scoring methods side-by-side.

    Parameters
    ----------
    corr_matrices : list of pd.DataFrame
        Correlation matrices (one per scoring method).
    method_names : list of str
        Names of the scoring methods.
    filename : str
        Output filename.
    """
    from .hallmark_genes import HALLMARK_SHORT

    n = len(corr_matrices)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
    if n == 1:
        axes = [axes]

    short = [HALLMARK_SHORT.get(c, c) for c in corr_matrices[0].columns]

    # Compute pairwise similarity scores for the title
    similarities = []
    for i in range(n):
        for j in range(i + 1, n):
            mask = np.tril(np.ones(corr_matrices[i].shape[0], dtype=bool), k=-1)
            vals_i = corr_matrices[i].values[mask]
            vals_j = corr_matrices[j].values[mask]
            from scipy.stats import pearsonr
            r, _ = pearsonr(vals_i, vals_j)
            similarities.append(f"{method_names[i]} vs {method_names[j]}: r={r:.3f}")

    for idx, (corr_mat, name) in enumerate(zip(corr_matrices, method_names)):
        ax = axes[idx]
        sns.heatmap(corr_mat.values, annot=True, fmt=".2f",
                    cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                    xticklabels=short, yticklabels=short,
                    square=True, linewidths=0.5, ax=ax,
                    cbar_kws={"label": "Correlation"})
        ax.set_title(name, fontweight="bold")

    sim_text = "Pairwise similarity: " + "; ".join(similarities)
    fig.suptitle("Scoring Method Comparison\n" + sim_text,
                 fontweight="bold", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_shared_gene_sensitivity(original_corr, reduced_corr,
                                  filename="fig15_shared_gene_sensitivity.pdf"):
    """Plot original vs shared-genes-removed correlation heatmaps.

    Parameters
    ----------
    original_corr : pd.DataFrame
        Original hallmark correlation matrix.
    reduced_corr : pd.DataFrame
        Correlation matrix after removing shared genes.
    filename : str
        Output filename.
    """
    from .hallmark_genes import HALLMARK_SHORT

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    short_orig = [HALLMARK_SHORT.get(c, c) for c in original_corr.columns]
    short_red = [HALLMARK_SHORT.get(c, c) for c in reduced_corr.columns]

    sns.heatmap(original_corr.values, annot=True, fmt=".2f",
                cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                xticklabels=short_orig, yticklabels=short_orig,
                square=True, linewidths=0.5, ax=axes[0],
                cbar_kws={"label": "Correlation"})
    axes[0].set_title("Original", fontweight="bold")

    sns.heatmap(reduced_corr.values, annot=True, fmt=".2f",
                cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                xticklabels=short_red, yticklabels=short_red,
                square=True, linewidths=0.5, ax=axes[1],
                cbar_kws={"label": "Correlation"})
    axes[1].set_title("Shared Genes Removed", fontweight="bold")

    # Compute element-wise correlation between the two matrices
    mask = np.tril(np.ones(original_corr.shape[0], dtype=bool), k=-1)
    vals_orig = original_corr.values[mask]
    vals_red = reduced_corr.values[mask]
    from scipy.stats import pearsonr
    r, _ = pearsonr(vals_orig, vals_red)

    fig.text(0.5, 0.01, f"Element-wise correlation between matrices: r = {r:.4f}",
             ha="center", fontsize=11,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    fig.suptitle("Shared Gene Sensitivity Analysis", fontweight="bold", fontsize=13)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_tissue_adjustment_comparison(unadjusted_corr, adjusted_corr,
                                      filename="fig16_tissue_adjustment.pdf"):
    """Plot unadjusted vs tissue-adjusted correlation heatmaps.

    Parameters
    ----------
    unadjusted_corr : pd.DataFrame
        Unadjusted hallmark correlation matrix.
    adjusted_corr : pd.DataFrame
        Tissue-adjusted hallmark correlation matrix.
    filename : str
        Output filename.
    """
    from .hallmark_genes import HALLMARK_SHORT

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    short_unadj = [HALLMARK_SHORT.get(c, c) for c in unadjusted_corr.columns]
    short_adj = [HALLMARK_SHORT.get(c, c) for c in adjusted_corr.columns]

    sns.heatmap(unadjusted_corr.values, annot=True, fmt=".2f",
                cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                xticklabels=short_unadj, yticklabels=short_unadj,
                square=True, linewidths=0.5, ax=axes[0],
                cbar_kws={"label": "Correlation"})
    axes[0].set_title("Unadjusted", fontweight="bold")

    sns.heatmap(adjusted_corr.values, annot=True, fmt=".2f",
                cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                xticklabels=short_adj, yticklabels=short_adj,
                square=True, linewidths=0.5, ax=axes[1],
                cbar_kws={"label": "Correlation"})
    axes[1].set_title("Tissue-Adjusted", fontweight="bold")

    # Compute similarity metric
    mask = np.tril(np.ones(unadjusted_corr.shape[0], dtype=bool), k=-1)
    vals_unadj = unadjusted_corr.values[mask]
    vals_adj = adjusted_corr.values[mask]
    from scipy.stats import pearsonr
    r, _ = pearsonr(vals_unadj, vals_adj)

    fig.text(0.5, 0.01,
             f"Similarity between unadjusted and adjusted: r = {r:.4f}",
             ha="center", fontsize=11,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    fig.suptitle("Tissue Adjustment Comparison", fontweight="bold", fontsize=13)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    fig.savefig(FIGURES_DIR / filename)
    plt.close(fig)
    print(f"  Saved: {filename}")
