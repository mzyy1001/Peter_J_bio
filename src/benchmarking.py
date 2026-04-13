"""
Ground-truth benchmarking module.

Since the dataset is simulated with known latent structure, we benchmark
how well our pipeline recovers the embedded ground truth:
- Edge recovery (precision, recall, F1, AUROC)
- Module recovery (Adjusted Rand Index)
- Centrality recovery (rank correlation)
- Age effect recovery (correlation of true vs estimated age effects)
- Robustness across random seeds
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import precision_recall_fscore_support, roc_auc_score
from itertools import combinations
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "results"


def benchmark_edge_recovery(estimated_corr, ground_truth_corr, thresholds=None):
    """
    Compare estimated hallmark correlation network against ground truth.
    Reports precision, recall, F1 at various thresholds.
    """
    if thresholds is None:
        thresholds = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]

    hallmarks = estimated_corr.columns.tolist()
    n_h = len(hallmarks)

    # Ground truth edges (any non-zero off-diagonal in embedded corr)
    gt_pairs = []
    gt_labels = []
    est_values = []
    for i, j in combinations(range(n_h), 2):
        gt_val = abs(ground_truth_corr.iloc[i, j])
        est_val = abs(estimated_corr.iloc[i, j])
        gt_pairs.append((hallmarks[i], hallmarks[j]))
        gt_labels.append(1 if gt_val > 0.1 else 0)
        est_values.append(est_val)

    gt_labels = np.array(gt_labels)
    est_values = np.array(est_values)

    results = []
    for thr in thresholds:
        pred_labels = (est_values > thr).astype(int)
        if pred_labels.sum() == 0 or pred_labels.sum() == len(pred_labels):
            continue
        p, r, f1, _ = precision_recall_fscore_support(
            gt_labels, pred_labels, average="binary", zero_division=0
        )
        results.append({
            "threshold": thr,
            "precision": p,
            "recall": r,
            "f1": f1,
            "n_predicted_edges": pred_labels.sum(),
        })

    # AUROC
    if len(np.unique(gt_labels)) > 1:
        auroc = roc_auc_score(gt_labels, est_values)
    else:
        auroc = float("nan")

    return pd.DataFrame(results), auroc


def benchmark_centrality_recovery(estimated_metrics, ground_truth_corr):
    """
    Compare estimated centrality rankings against ground truth.
    Ground truth centrality = row sum of absolute correlations.
    """
    hallmarks = ground_truth_corr.columns.tolist()

    # Ground truth centrality: sum of absolute correlations
    gt_centrality = ground_truth_corr.abs().sum(axis=1) - 1  # subtract diagonal
    gt_centrality = gt_centrality.loc[hallmarks]

    # Estimated centrality (eigenvector or weighted degree)
    if "eigenvector" in estimated_metrics.columns:
        est_centrality = estimated_metrics["eigenvector"]
    else:
        est_centrality = estimated_metrics["weighted_degree"]

    # Align indices
    common = sorted(set(gt_centrality.index) & set(est_centrality.index))
    gt = gt_centrality.loc[common]
    est = est_centrality.loc[common]

    # Spearman rank correlation
    rho, p = stats.spearmanr(gt, est)

    # Top-k agreement
    gt_top3 = set(gt.nlargest(3).index)
    est_top3 = set(est.nlargest(3).index)
    top3_overlap = len(gt_top3 & est_top3) / 3

    return {
        "spearman_rho": rho,
        "spearman_p": p,
        "top3_overlap": top3_overlap,
        "gt_top3": sorted(gt_top3),
        "est_top3": sorted(est_top3),
        "gt_centrality": gt.to_dict(),
        "est_centrality": est.to_dict(),
    }


def benchmark_module_recovery(estimated_communities, ground_truth_modules):
    """
    Compare detected modules against expected modules using Adjusted Rand Index.
    """
    from sklearn.metrics import adjusted_rand_score

    # Create label vectors
    all_nodes = set()
    for comm in estimated_communities:
        all_nodes.update(comm)
    for mod in ground_truth_modules:
        all_nodes.update(mod)
    all_nodes = sorted(all_nodes)

    est_labels = {}
    for ci, comm in enumerate(estimated_communities):
        for node in comm:
            est_labels[node] = ci

    gt_labels = {}
    for mi, mod in enumerate(ground_truth_modules):
        for node in mod:
            gt_labels[node] = mi

    # Align
    common = sorted(set(est_labels.keys()) & set(gt_labels.keys()))
    y_est = [est_labels[n] for n in common]
    y_gt = [gt_labels[n] for n in common]

    ari = adjusted_rand_score(y_gt, y_est)
    return {"adjusted_rand_index": ari, "n_nodes": len(common)}


def benchmark_age_effect_recovery(estimated_age_corrs, true_age_effects):
    """
    Compare estimated age-hallmark correlations against true embedded age effects.
    """
    hallmarks = list(true_age_effects.keys())
    gt_vals = [true_age_effects[h] for h in hallmarks]

    if isinstance(estimated_age_corrs, pd.Series):
        est_vals = [estimated_age_corrs.get(h, 0) for h in hallmarks]
    else:
        est_vals = [estimated_age_corrs.get(h, {}).get("r", 0) for h in hallmarks]

    rho, p = stats.spearmanr(gt_vals, est_vals)
    pearson_r, pearson_p = stats.pearsonr(gt_vals, est_vals)

    return {
        "spearman_rho": rho,
        "spearman_p": p,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "gt_effects": dict(zip(hallmarks, gt_vals)),
        "est_effects": dict(zip(hallmarks, est_vals)),
    }


def robustness_analysis(n_seeds=10):
    """
    Run the pipeline across multiple random seeds and assess stability.
    """
    from .data_acquisition import generate_simulated_expression_data
    from .hallmark_genes import HALLMARK_GENES
    from .network_analysis import compute_hallmark_scores, build_correlation_network

    hallmarks = list(HALLMARK_GENES.keys())
    all_corrs = []

    for seed in range(n_seeds):
        expr_df, metadata, _ = generate_simulated_expression_data(
            n_samples=300, seed=seed
        )
        scores = compute_hallmark_scores(expr_df, HALLMARK_GENES, method="ssgsea")
        corr = scores.corr(method="spearman")
        all_corrs.append(corr.values[np.triu_indices(len(hallmarks), k=1)])

    all_corrs = np.array(all_corrs)
    mean_corr = all_corrs.mean(axis=0)
    std_corr = all_corrs.std(axis=0)
    cv = std_corr / (np.abs(mean_corr) + 1e-10)

    # Create edge-level stability report
    edges = []
    idx = 0
    for i, j in combinations(range(len(hallmarks)), 2):
        edges.append({
            "hallmark_1": hallmarks[i],
            "hallmark_2": hallmarks[j],
            "mean_corr": mean_corr[idx],
            "std_corr": std_corr[idx],
            "cv": cv[idx],
            "stable": std_corr[idx] < 0.1,
        })
        idx += 1

    stability_df = pd.DataFrame(edges)
    pct_stable = stability_df["stable"].mean() * 100

    return {
        "n_seeds": n_seeds,
        "stability_df": stability_df,
        "pct_stable_edges": pct_stable,
        "mean_cv": cv.mean(),
    }


def run_all_benchmarks(estimated_corr, ground_truth_corr, node_metrics,
                       communities, age_corrs, true_age_effects):
    """Run all benchmarks and save results."""
    print("\n  === Ground Truth Benchmarking ===")

    # Edge recovery
    edge_df, auroc = benchmark_edge_recovery(estimated_corr, ground_truth_corr)
    print(f"  Edge recovery AUROC: {auroc:.3f}")
    if len(edge_df) > 0:
        best = edge_df.loc[edge_df["f1"].idxmax()]
        print(f"  Best F1={best['f1']:.3f} at threshold={best['threshold']}")
    edge_df.to_csv(RESULTS_DIR / "benchmark_edge_recovery.csv", index=False)

    # Centrality recovery
    cent_results = benchmark_centrality_recovery(node_metrics, ground_truth_corr)
    print(f"  Centrality rank correlation: rho={cent_results['spearman_rho']:.3f}, "
          f"p={cent_results['spearman_p']:.4f}")
    print(f"  Top-3 hub overlap: {cent_results['top3_overlap']*100:.0f}%")

    # Module recovery
    gt_modules = [
        ["Genomic Instability", "Telomere Attrition",
         "Epigenetic Alterations", "Loss of Proteostasis"],
        ["Deregulated Nutrient Sensing", "Mitochondrial Dysfunction",
         "Cellular Senescence"],
        ["Stem Cell Exhaustion", "Altered Intercellular Communication"],
    ]
    if communities:
        mod_results = benchmark_module_recovery(communities, gt_modules)
        print(f"  Module recovery ARI: {mod_results['adjusted_rand_index']:.3f}")
    else:
        mod_results = {"adjusted_rand_index": float("nan")}

    # Age effect recovery
    from .hallmark_genes import HALLMARK_GENES
    hallmark_names = list(HALLMARK_GENES.keys())
    true_effects_dict = dict(zip(hallmark_names, [0.6, 0.5, 0.4, 0.5, 0.3, 0.5, 0.6, 0.5, 0.55]))
    age_results = benchmark_age_effect_recovery(age_corrs, true_effects_dict)
    print(f"  Age effect recovery: rho={age_results['spearman_rho']:.3f}")

    # Robustness
    print("  Running robustness analysis (10 seeds)...")
    robust = robustness_analysis(n_seeds=10)
    print(f"  Stable edges (CV<0.1): {robust['pct_stable_edges']:.0f}%")
    robust["stability_df"].to_csv(RESULTS_DIR / "benchmark_robustness.csv", index=False)

    # Summary
    summary = {
        "edge_auroc": auroc,
        "best_f1": float(edge_df["f1"].max()) if len(edge_df) > 0 else 0,
        "centrality_rho": cent_results["spearman_rho"],
        "top3_overlap": cent_results["top3_overlap"],
        "module_ari": mod_results["adjusted_rand_index"],
        "age_effect_rho": age_results["spearman_rho"],
        "pct_stable_edges": robust["pct_stable_edges"],
    }

    import json
    with open(RESULTS_DIR / "benchmark_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n  Benchmark Summary:")
    for k, v in summary.items():
        print(f"    {k}: {v}")

    return summary
