"""
Sensitivity analyses for the hallmarks-of-aging interaction network.

Implements three reviewer-requested analyses:
1. Shared-gene removal sensitivity
2. Alternative scoring method comparison (ssGSEA, mean z-score, pathway PCA)
3. Tissue-covariate adjustment

Each analysis recomputes hallmark activity scores and the Spearman correlation
network, then quantifies how much the network changes relative to the baseline.
"""

import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations
from pathlib import Path
from sklearn.decomposition import PCA

from src.hallmark_genes import HALLMARK_GENES, get_gene_hallmark_map
from src.network_analysis import compute_hallmark_scores, build_correlation_network


RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _upper_triangle(corr_matrix):
    """Extract the upper-triangle values (excluding diagonal) from a square DataFrame."""
    idx = np.triu_indices_from(corr_matrix.values, k=1)
    return corr_matrix.values[idx]


def _compare_matrices(mat_a, mat_b, label_a="A", label_b="B"):
    """Element-wise Pearson r between the upper triangles of two correlation matrices."""
    a = _upper_triangle(mat_a)
    b = _upper_triangle(mat_b)
    r, p = stats.pearsonr(a, b)
    return {"comparison": f"{label_a}_vs_{label_b}", "pearson_r": r, "pvalue": p}


# ---------------------------------------------------------------------------
# 1. Shared-gene sensitivity analysis
# ---------------------------------------------------------------------------

def shared_gene_sensitivity(expr_df, hallmark_genes=None):
    """
    Remove genes that appear in 2+ hallmarks, recompute scores and the
    Spearman correlation network, and compare to the original network.

    Parameters
    ----------
    expr_df : pd.DataFrame
        Samples x genes expression matrix.
    hallmark_genes : dict or None
        Hallmark gene sets. Defaults to HALLMARK_GENES.

    Returns
    -------
    dict with keys: original_corr, reduced_corr, correlation_similarity,
    n_removed_genes, affected_hallmarks, removed_genes, reduced_gene_sets.
    """
    if hallmark_genes is None:
        hallmark_genes = HALLMARK_GENES

    print("[Sensitivity 1] Shared-gene removal analysis")
    print("  Identifying multi-hallmark genes ...")

    gene_map = get_gene_hallmark_map()
    shared_genes = {g for g, hmarks in gene_map.items() if len(hmarks) >= 2}
    affected = sorted({h for g in shared_genes for h in gene_map[g]})

    print(f"  Found {len(shared_genes)} shared genes across {len(affected)} hallmarks")

    # Build reduced gene sets
    reduced_genes = {}
    for hallmark, genes in hallmark_genes.items():
        reduced_genes[hallmark] = [g for g in genes if g not in shared_genes]

    empty_hallmarks = [h for h, g in reduced_genes.items() if len(g) == 0]
    if empty_hallmarks:
        print(f"  WARNING: hallmarks with zero genes after removal: {empty_hallmarks}")

    # Original scores and correlation
    print("  Computing original hallmark scores ...")
    orig_scores = compute_hallmark_scores(expr_df, hallmark_genes, method="ssgsea")
    _, orig_corr, _ = build_correlation_network(orig_scores, method="spearman",
                                                 threshold=0.0)

    # Reduced scores and correlation
    print("  Computing reduced hallmark scores (shared genes removed) ...")
    red_scores = compute_hallmark_scores(expr_df, reduced_genes, method="ssgsea")
    _, red_corr, _ = build_correlation_network(red_scores, method="spearman",
                                                threshold=0.0)

    # Compare
    sim = _compare_matrices(orig_corr, red_corr, "original", "reduced")
    print(f"  Correlation similarity (Pearson r of upper triangles): {sim['pearson_r']:.4f}")

    return {
        "original_corr": orig_corr,
        "reduced_corr": red_corr,
        "correlation_similarity": sim["pearson_r"],
        "similarity_pvalue": sim["pvalue"],
        "n_removed_genes": len(shared_genes),
        "removed_genes": sorted(shared_genes),
        "affected_hallmarks": affected,
        "reduced_gene_sets": reduced_genes,
    }


# ---------------------------------------------------------------------------
# 2. Alternative scoring method comparison
# ---------------------------------------------------------------------------

def _score_mean_zscore(expr_df, hallmark_genes):
    """Mean z-score scoring: average of z-scored gene expression per gene set."""
    return compute_hallmark_scores(expr_df, hallmark_genes, method="mean")


def _score_ssgsea(expr_df, hallmark_genes):
    """Rank-based ssGSEA scoring (existing implementation)."""
    return compute_hallmark_scores(expr_df, hallmark_genes, method="ssgsea")


def _score_pathway_pca(expr_df, hallmark_genes):
    """First principal component of each gene set as the hallmark activity score."""
    expr_z = (expr_df - expr_df.mean()) / expr_df.std()
    hallmarks = list(hallmark_genes.keys())
    scores = pd.DataFrame(index=expr_df.index, columns=hallmarks, dtype=float)

    for h in hallmarks:
        genes = [g for g in hallmark_genes[h] if g in expr_z.columns]
        if len(genes) < 2:
            # Not enough genes for PCA; fall back to mean
            scores[h] = expr_z[genes].mean(axis=1) if genes else 0.0
            continue
        subset = expr_z[genes].dropna(axis=1)
        if subset.shape[1] < 2:
            scores[h] = subset.mean(axis=1) if subset.shape[1] else 0.0
            continue
        pca = PCA(n_components=1, random_state=42)
        pc1 = pca.fit_transform(subset.values)
        # Anchor PC1 sign: correlate with mean z-score of hallmark genes;
        # flip if negative so PCA scores point in the same direction as mean expression.
        mean_z = subset.mean(axis=1)
        if np.corrcoef(pc1[:, 0], mean_z)[0, 1] < 0:
            pc1 = -pc1
        scores[h] = pc1[:, 0]

    return scores.astype(float)


def alternative_scoring_comparison(expr_df, hallmark_genes=None):
    """
    Compute hallmark scores with three methods (ssGSEA, mean z-score,
    pathway PCA), build Spearman correlation networks, and compare all
    pairs of correlation matrices.

    Parameters
    ----------
    expr_df : pd.DataFrame
        Samples x genes expression matrix.
    hallmark_genes : dict or None
        Hallmark gene sets. Defaults to HALLMARK_GENES.

    Returns
    -------
    dict with keys: corr_ssgsea, corr_meanz, corr_pca, pairwise_similarity
    (list of dicts with comparison, pearson_r, pvalue).
    """
    if hallmark_genes is None:
        hallmark_genes = HALLMARK_GENES

    print("[Sensitivity 2] Alternative scoring method comparison")

    methods = {
        "ssgsea": _score_ssgsea,
        "meanz": _score_mean_zscore,
        "pca": _score_pathway_pca,
    }

    corr_matrices = {}
    for name, func in methods.items():
        print(f"  Computing scores with method: {name} ...")
        scores = func(expr_df, hallmark_genes)
        _, corr, _ = build_correlation_network(scores, method="spearman",
                                                threshold=0.0)
        corr_matrices[name] = corr

    # Pairwise comparison of correlation matrices
    print("  Comparing correlation matrices pairwise ...")
    names = list(corr_matrices.keys())
    pairwise = []
    for a, b in combinations(names, 2):
        sim = _compare_matrices(corr_matrices[a], corr_matrices[b], a, b)
        pairwise.append(sim)
        print(f"    {a} vs {b}: Pearson r = {sim['pearson_r']:.4f}")

    return {
        "corr_ssgsea": corr_matrices["ssgsea"],
        "corr_meanz": corr_matrices["meanz"],
        "corr_pca": corr_matrices["pca"],
        "pairwise_similarity": pairwise,
    }


# ---------------------------------------------------------------------------
# 3. Tissue-covariate adjustment
# ---------------------------------------------------------------------------

def tissue_covariate_adjustment(expr_df, metadata, tissue_column="tissue",
                                hallmark_genes=None):
    """
    Regress out tissue effects from each gene using OLS with tissue dummy
    variables, recompute hallmark scores from residuals, and compare the
    adjusted network to the unadjusted network.

    Parameters
    ----------
    expr_df : pd.DataFrame
        Samples x genes expression matrix. Index must align with metadata.
    metadata : pd.DataFrame
        Must contain a column specified by *tissue_column* and have the same
        index (or a superset) as expr_df.
    tissue_column : str
        Column name in metadata that holds tissue labels.
    hallmark_genes : dict or None
        Hallmark gene sets. Defaults to HALLMARK_GENES.

    Returns
    -------
    dict with keys: adjusted_scores, adjusted_corr, unadjusted_corr,
    similarity, similarity_pvalue, tissue_residual_expr.
    """
    if hallmark_genes is None:
        hallmark_genes = HALLMARK_GENES

    print("[Sensitivity 3] Tissue-covariate adjustment")

    # Align metadata and expression
    # Handle case where metadata has sample_id as column (not index)
    if "sample_id" in metadata.columns and metadata.index.name != "sample_id":
        meta_indexed = metadata.set_index("sample_id")
    else:
        meta_indexed = metadata
    common = expr_df.index.intersection(meta_indexed.index)
    if len(common) == 0:
        raise ValueError("No overlapping samples between expr_df and metadata.")
    expr_aligned = expr_df.loc[common]
    meta_aligned = meta_indexed.loc[common]

    tissues = meta_aligned[tissue_column]
    unique_tissues = tissues.unique()
    print(f"  Found {len(unique_tissues)} tissue types, {len(common)} samples")

    if len(unique_tissues) < 2:
        print("  WARNING: fewer than 2 tissue types; skipping adjustment.")
        scores = compute_hallmark_scores(expr_aligned, hallmark_genes, method="ssgsea")
        _, corr, _ = build_correlation_network(scores, method="spearman", threshold=0.0)
        return {
            "adjusted_scores": scores,
            "adjusted_corr": corr,
            "unadjusted_corr": corr,
            "similarity": 1.0,
            "similarity_pvalue": 0.0,
            "tissue_residual_expr": expr_aligned,
        }

    # Build design matrix (tissue dummies, drop first to avoid collinearity)
    print("  Building OLS design matrix ...")
    dummies = pd.get_dummies(tissues, drop_first=True, dtype=float)
    # Add intercept
    dummies.insert(0, "_intercept", 1.0)
    X = dummies.values  # (n_samples, n_dummies)

    # Regress out tissue from every gene
    print("  Regressing out tissue effects from each gene ...")
    expr_vals = expr_aligned.values  # (n_samples, n_genes)
    # OLS: beta = (X'X)^-1 X' Y,  residuals = Y - X beta
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ X.T @ expr_vals
    residuals = expr_vals - X @ beta

    expr_residual = pd.DataFrame(residuals, index=expr_aligned.index,
                                 columns=expr_aligned.columns)

    # Unadjusted scores and network
    print("  Computing unadjusted hallmark scores ...")
    unadj_scores = compute_hallmark_scores(expr_aligned, hallmark_genes, method="ssgsea")
    _, unadj_corr, _ = build_correlation_network(unadj_scores, method="spearman",
                                                  threshold=0.0)

    # Adjusted scores and network
    print("  Computing tissue-adjusted hallmark scores ...")
    adj_scores = compute_hallmark_scores(expr_residual, hallmark_genes, method="ssgsea")
    _, adj_corr, _ = build_correlation_network(adj_scores, method="spearman",
                                                threshold=0.0)

    sim = _compare_matrices(unadj_corr, adj_corr, "unadjusted", "adjusted")
    print(f"  Adjusted vs unadjusted similarity (Pearson r): {sim['pearson_r']:.4f}")

    return {
        "adjusted_scores": adj_scores,
        "adjusted_corr": adj_corr,
        "unadjusted_corr": unadj_corr,
        "similarity": sim["pearson_r"],
        "similarity_pvalue": sim["pvalue"],
        "tissue_residual_expr": expr_residual,
    }


# ---------------------------------------------------------------------------
# 4. Run all sensitivity analyses
# ---------------------------------------------------------------------------

def run_all_sensitivity(expr_df, metadata=None, tissue_column="tissue",
                        hallmark_genes=None, output_dir=None):
    """
    Execute all three sensitivity analyses and save key results to CSV.

    Parameters
    ----------
    expr_df : pd.DataFrame
        Samples x genes expression matrix.
    metadata : pd.DataFrame or None
        Sample metadata. Required for tissue adjustment (must contain the
        tissue column and share index with expr_df).
    tissue_column : str
        Column in metadata with tissue labels.
    hallmark_genes : dict or None
        Hallmark gene sets. Defaults to HALLMARK_GENES.
    output_dir : str or Path or None
        Directory for output files. Defaults to project results/.

    Returns
    -------
    dict with keys 'shared_gene', 'scoring_methods', 'tissue_adjustment',
    each containing the respective analysis result dict.
    """
    if hallmark_genes is None:
        hallmark_genes = HALLMARK_GENES
    if output_dir is None:
        output_dir = RESULTS_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    # --- 1. Shared-gene sensitivity ---
    print("=" * 60)
    res1 = shared_gene_sensitivity(expr_df, hallmark_genes)
    results["shared_gene"] = res1
    res1["original_corr"].to_csv(output_dir / "sensitivity_shared_gene_original_corr.csv")
    res1["reduced_corr"].to_csv(output_dir / "sensitivity_shared_gene_reduced_corr.csv")

    # --- 2. Alternative scoring methods ---
    print("=" * 60)
    res2 = alternative_scoring_comparison(expr_df, hallmark_genes)
    results["scoring_methods"] = res2
    res2["corr_ssgsea"].to_csv(output_dir / "sensitivity_scoring_corr_ssgsea.csv")
    res2["corr_meanz"].to_csv(output_dir / "sensitivity_scoring_corr_meanz.csv")
    res2["corr_pca"].to_csv(output_dir / "sensitivity_scoring_corr_pca.csv")
    pd.DataFrame(res2["pairwise_similarity"]).to_csv(
        output_dir / "sensitivity_scoring_pairwise_similarity.csv", index=False
    )

    # --- 3. Tissue-covariate adjustment ---
    if metadata is not None and tissue_column in metadata.columns:
        print("=" * 60)
        res3 = tissue_covariate_adjustment(expr_df, metadata, tissue_column,
                                           hallmark_genes)
        results["tissue_adjustment"] = res3
        res3["adjusted_corr"].to_csv(
            output_dir / "sensitivity_tissue_adjusted_corr.csv"
        )
        res3["unadjusted_corr"].to_csv(
            output_dir / "sensitivity_tissue_unadjusted_corr.csv"
        )
    else:
        print("=" * 60)
        print("[Sensitivity 3] Skipped: no metadata or tissue column not found.")
        results["tissue_adjustment"] = None

    # --- Summary ---
    print("=" * 60)
    summary = {
        "analysis": [],
        "key_metric": [],
        "value": [],
    }
    summary["analysis"].append("shared_gene_removal")
    summary["key_metric"].append("correlation_similarity")
    summary["value"].append(res1["correlation_similarity"])

    summary["analysis"].append("shared_gene_removal")
    summary["key_metric"].append("n_removed_genes")
    summary["value"].append(res1["n_removed_genes"])

    for pair in res2["pairwise_similarity"]:
        summary["analysis"].append("scoring_comparison")
        summary["key_metric"].append(pair["comparison"])
        summary["value"].append(pair["pearson_r"])

    if results["tissue_adjustment"] is not None:
        summary["analysis"].append("tissue_adjustment")
        summary["key_metric"].append("adjusted_vs_unadjusted_similarity")
        summary["value"].append(results["tissue_adjustment"]["similarity"])

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_dir / "sensitivity_summary.csv", index=False)
    print(f"\nSensitivity summary saved to {output_dir / 'sensitivity_summary.csv'}")
    print(summary_df.to_string(index=False))

    return results
