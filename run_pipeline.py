#!/usr/bin/env python3
"""
Main pipeline: Hallmarks of Aging Interaction Analysis

Runs the complete bioinformatics analysis:
1. Data acquisition & preprocessing
2. Hallmark activity score computation
3. Network construction (correlation, partial correlation, causal)
4. Machine learning models (age prediction, cross-hallmark prediction)
5. Visualization
6. Results summary export
"""

import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

# Project paths
PROJECT_DIR = Path(__file__).parent
sys.path.insert(0, str(PROJECT_DIR))

from src.hallmark_genes import (
    HALLMARK_GENES, HALLMARK_CATEGORIES, HALLMARK_SHORT,
    get_all_genes, get_gene_hallmark_map, get_shared_genes
)
from src.data_acquisition import generate_simulated_expression_data, download_genage_genes
from src.network_analysis import (
    build_gene_overlap_network, compute_hallmark_scores,
    build_correlation_network, partial_correlation_network,
    compute_network_metrics, age_stratified_networks,
    causal_inference_pc
)
from src.ml_models import (
    age_group_classifier, biological_age_estimator,
    cross_hallmark_prediction, hallmark_pca_analysis,
    tissue_specific_analysis, interaction_strength_model
)
from src.benchmarking import run_all_benchmarks
from src.sensitivity import run_all_sensitivity
from src.visualization import (
    plot_hallmark_correlation_heatmap, plot_network_graph,
    plot_age_correlation_barplot, plot_feature_importance,
    plot_pca_scatter, plot_biological_age,
    plot_cross_hallmark_prediction, plot_interaction_changes,
    plot_partial_correlation_network, plot_variance_explained,
    plot_causal_graph,
    plot_scoring_method_comparison,
    plot_shared_gene_sensitivity,
    plot_tissue_adjustment_comparison
)

RESULTS_DIR = PROJECT_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)
DATA_DIR = PROJECT_DIR / "data"


def main():
    print("=" * 70)
    print("HALLMARKS OF AGING INTERACTION ANALYSIS PIPELINE")
    print("=" * 70)

    # =====================================================
    # STEP 1: Data Acquisition
    # =====================================================
    print("\n[1/6] Data Acquisition & Preprocessing")
    print("-" * 40)

    # Download GenAge
    genage = download_genage_genes()

    # Generate simulated expression data (300 samples, 5 tissues)
    expr_df, metadata, gt_scores = generate_simulated_expression_data(n_samples=300, seed=42)

    # Shared genes analysis
    shared = get_shared_genes()
    print(f"\n  Cross-hallmark genes (in 2+ hallmarks): {len(shared)}")
    for gene, hallmarks in sorted(shared.items(), key=lambda x: -len(x[1]))[:10]:
        print(f"    {gene}: {', '.join(HALLMARK_SHORT[h] for h in hallmarks)}")

    # =====================================================
    # STEP 2: Hallmark Activity Scoring
    # =====================================================
    print("\n[2/6] Computing Hallmark Activity Scores")
    print("-" * 40)

    hallmark_scores = compute_hallmark_scores(expr_df, HALLMARK_GENES, method="ssgsea")
    print(f"  Score matrix: {hallmark_scores.shape}")
    print(f"  Score ranges:")
    for h in hallmark_scores.columns:
        print(f"    {HALLMARK_SHORT[h]}: [{hallmark_scores[h].min():.3f}, {hallmark_scores[h].max():.3f}]")

    hallmark_scores.to_csv(RESULTS_DIR / "computed_hallmark_scores.csv")

    # =====================================================
    # STEP 3: Network Analysis
    # =====================================================
    print("\n[3/6] Network Analysis")
    print("-" * 40)

    # 3a. Gene overlap network
    print("  3a. Gene overlap network...")
    overlap_G = build_gene_overlap_network(HALLMARK_GENES)
    overlap_metrics = compute_network_metrics(overlap_G)
    print(f"      Edges: {overlap_metrics['n_edges']}, Density: {overlap_metrics['density']:.3f}")

    # 3b. Correlation network
    print("  3b. Correlation network...")
    corr_G, corr_matrix, p_matrix = build_correlation_network(
        hallmark_scores, method="spearman", threshold=0.1
    )
    corr_metrics = compute_network_metrics(corr_G)
    print(f"      Significant edges: {corr_metrics['n_edges']}")
    corr_matrix.to_csv(RESULTS_DIR / "hallmark_correlation_matrix.csv")
    p_matrix.to_csv(RESULTS_DIR / "hallmark_pvalue_matrix.csv")

    # 3c. Partial correlation network
    print("  3c. Partial correlation network...")
    pcorr_G, pcorr_df = partial_correlation_network(hallmark_scores, threshold=0.08)
    pcorr_metrics = compute_network_metrics(pcorr_G)
    print(f"      Direct edges: {pcorr_metrics['n_edges']}")
    pcorr_df.to_csv(RESULTS_DIR / "partial_correlation_matrix.csv")

    # 3d. Age-stratified networks
    print("  3d. Age-stratified networks...")
    age_networks, rewiring_df = age_stratified_networks(hallmark_scores, metadata)
    for label, net_data in age_networks.items():
        n_edges = net_data["graph"].number_of_edges()
        print(f"      {label} (n={net_data['n_samples']}): {n_edges} edges")
    if rewiring_df is not None and len(rewiring_df) > 0:
        sig_col = "significant_fdr" if "significant_fdr" in rewiring_df.columns else "significant"
        sig_rewire = rewiring_df[rewiring_df[sig_col] == True]
        print(f"      Significant rewiring edges (FDR<0.05): {len(sig_rewire)}")
        rewiring_df.to_csv(RESULTS_DIR / "age_rewiring_fdr.csv", index=False)

    # 3e. Causal inference
    print("  3e. Causal structure learning (PC algorithm)...")
    causal_G, adj_matrix = causal_inference_pc(hallmark_scores, alpha=0.05)
    n_directed = sum(1 for _, _, d in causal_G.edges(data=True) if d.get("type") == "directed")
    print(f"      Directed edges: {n_directed}, Total edges: {causal_G.number_of_edges()}")

    # Node-level metrics
    if "node_metrics" in corr_metrics:
        node_df = corr_metrics["node_metrics"]
        print("\n  Network centrality (correlation network):")
        for h in node_df.index:
            print(f"    {HALLMARK_SHORT.get(h, h)}: degree={node_df.loc[h, 'weighted_degree']:.2f}, "
                  f"betweenness={node_df.loc[h, 'betweenness']:.3f}, "
                  f"eigenvector={node_df.loc[h, 'eigenvector']:.3f}")
        node_df.to_csv(RESULTS_DIR / "network_node_metrics.csv")

    # =====================================================
    # STEP 4: Machine Learning Models
    # =====================================================
    print("\n[4/6] Machine Learning Models")
    print("-" * 40)

    # 4a. Age group classifier
    print("  4a. Age group classification...")
    clf_results = age_group_classifier(hallmark_scores, metadata)
    print(f"      RF Accuracy: {clf_results['rf_accuracy']:.3f} +/- {clf_results['rf_std']:.3f}")
    print(f"      GB Accuracy: {clf_results['gb_accuracy']:.3f} +/- {clf_results['gb_std']:.3f}")
    print(f"      LR Accuracy: {clf_results['lr_accuracy']:.3f} +/- {clf_results['lr_std']:.3f}")
    clf_results["feature_importances"].to_csv(RESULTS_DIR / "age_clf_importances.csv", index=False)

    # 4b. Biological age estimator
    print("  4b. Biological age estimation...")
    age_results = biological_age_estimator(hallmark_scores, metadata)
    print(f"      GB R² (CV): {age_results['gb_r2_cv']:.3f} +/- {age_results['gb_r2_cv_std']:.3f}")
    print(f"      GB MAE: {age_results['gb_mae']:.1f} years")
    print(f"      EN R²: {age_results['en_r2']:.3f}, MAE: {age_results['en_mae']:.1f} years")
    age_results["age_acceleration"].to_csv(RESULTS_DIR / "age_acceleration.csv", index=False)
    age_results["en_coefficients"].to_csv(RESULTS_DIR / "age_en_coefficients.csv", index=False)
    age_results["gb_importances"].to_csv(RESULTS_DIR / "age_gb_importances.csv", index=False)

    # 4c. Cross-hallmark prediction
    print("  4c. Cross-hallmark prediction...")
    cross_pred = cross_hallmark_prediction(hallmark_scores)
    print("      Hallmark predictability (R² from other hallmarks):")
    for _, row in cross_pred.iterrows():
        print(f"        {HALLMARK_SHORT[row['target_hallmark']]}: "
              f"R²={row['rf_r2_cv']:.3f} (top predictor: {HALLMARK_SHORT[row['top_predictor']]})")
    cross_pred.to_csv(RESULTS_DIR / "cross_hallmark_prediction.csv", index=False)

    # 4d. PCA analysis
    print("  4d. PCA / dimensionality reduction...")
    pca_results = hallmark_pca_analysis(hallmark_scores, metadata)
    print(f"      PC1: {pca_results['explained_variance_ratio'][0]*100:.1f}%, "
          f"PC2: {pca_results['explained_variance_ratio'][1]*100:.1f}%, "
          f"PC3: {pca_results['explained_variance_ratio'][2]*100:.1f}%")
    pca_results["components"].to_csv(RESULTS_DIR / "pca_loadings.csv")

    # 4e. Tissue-specific analysis
    print("  4e. Tissue-specific analysis...")
    tissue_results = tissue_specific_analysis(hallmark_scores, metadata)
    for tissue, res in tissue_results.items():
        top_h = res["age_correlations"]["r"].abs().idxmax()
        top_r = res["age_correlations"].loc[top_h, "r"]
        print(f"      {tissue} (n={res['n_samples']}): strongest age-corr = "
              f"{HALLMARK_SHORT.get(top_h, top_h)} (r={top_r:.3f})")

    # 4f. Interaction strength model
    print("  4f. Interaction strength & age-dependent changes...")
    mi_df, interaction_df = interaction_strength_model(hallmark_scores, metadata)
    mi_df.to_csv(RESULTS_DIR / "mutual_information_matrix.csv")
    interaction_df.to_csv(RESULTS_DIR / "interaction_age_changes.csv", index=False)
    sig_changes = interaction_df[interaction_df["significant"]]
    print(f"      Significant age-dependent changes: {len(sig_changes)} pairs")

    # =====================================================
    # STEP 5: Visualization
    # =====================================================
    print("\n[5/6] Generating Figures")
    print("-" * 40)

    node_met = corr_metrics.get("node_metrics", None)

    plot_hallmark_correlation_heatmap(corr_matrix, p_matrix,
        title="Spearman Correlation Between Hallmark Activity Scores",
        filename="fig1_hallmark_correlation_heatmap.pdf")

    plot_network_graph(corr_G, title="Hallmark Correlation Network",
        filename="fig2_correlation_network.pdf", node_metrics=node_met)

    plot_partial_correlation_network(pcorr_G, pcorr_df,
        filename="fig3_partial_correlation_network.pdf")

    plot_age_correlation_barplot(hallmark_scores, metadata,
        filename="fig4_age_hallmark_correlations.pdf")

    plot_feature_importance(clf_results["feature_importances"],
        title="Hallmark Importance for Age Group Prediction",
        filename="fig5_feature_importance.pdf")

    plot_pca_scatter(pca_results, filename="fig6_pca_scatter.pdf")

    plot_variance_explained(pca_results, filename="fig7_variance_explained.pdf")

    plot_biological_age(age_results, filename="fig8_biological_age.pdf")

    plot_cross_hallmark_prediction(cross_pred,
        filename="fig9_cross_hallmark_prediction.pdf")

    plot_interaction_changes(interaction_df,
        filename="fig10_interaction_age_changes.pdf")

    plot_causal_graph(causal_G, filename="fig11_causal_network.pdf")

    # Age-stratified correlation heatmaps
    for label, net_data in age_networks.items():
        plot_hallmark_correlation_heatmap(
            net_data["corr"], net_data["pval"],
            title=f"Hallmark Correlations - {label} (n={net_data['n_samples']})",
            filename=f"fig12_{label.lower()}_correlation_heatmap.pdf"
        )

    # Gene overlap network
    plot_network_graph(overlap_G,
        title="Gene Overlap Network Between Hallmarks",
        filename="fig13_gene_overlap_network.pdf", layout="circular")

    # =====================================================
    # STEP 5b: Ground Truth Benchmarking
    # =====================================================
    print("\n[5b/7] Ground Truth Benchmarking")
    print("-" * 40)

    # Load ground truth correlation matrix
    # Prefer observed latent correlations (includes age+tissue effects) if available
    gt_observed_path = DATA_DIR / "ground_truth_observed_latent_correlations.csv"
    gt_base_path = DATA_DIR / "ground_truth_base_hallmark_correlations.csv"
    if gt_observed_path.exists():
        gt_corr = pd.read_csv(gt_observed_path, index_col=0)
        print("  Using observed latent correlations as benchmark target")
    else:
        gt_corr = pd.read_csv(gt_base_path, index_col=0)
        print("  Using base latent correlations as benchmark target")

    # Age correlations for benchmarking
    from scipy import stats as sp_stats
    ages = metadata.set_index("sample_id").loc[hallmark_scores.index, "age"]
    age_corrs_dict = {}
    for h in hallmark_scores.columns:
        r, p = sp_stats.spearmanr(ages, hallmark_scores[h])
        age_corrs_dict[h] = r

    benchmark_summary = run_all_benchmarks(
        estimated_corr=corr_matrix,
        ground_truth_corr=gt_corr,
        node_metrics=corr_metrics.get("node_metrics", pd.DataFrame()),
        communities=corr_metrics.get("communities", []),
        age_corrs=pd.Series(age_corrs_dict),
        true_age_effects=dict(zip(
            list(HALLMARK_GENES.keys()),
            [0.6, 0.5, 0.4, 0.5, 0.3, 0.5, 0.6, 0.5, 0.55]
        )),
    )

    # =====================================================
    # STEP 5c: Sensitivity Analyses
    # =====================================================
    print("\n[5c/7] Sensitivity Analyses")
    print("-" * 40)

    sensitivity_results = run_all_sensitivity(
        expr_df, metadata=metadata, tissue_column="tissue",
        hallmark_genes=HALLMARK_GENES, output_dir=str(RESULTS_DIR)
    )

    # Extract sub-results
    sg = sensitivity_results["shared_gene"]
    sm = sensitivity_results["scoring_methods"]
    ta = sensitivity_results.get("tissue_adjustment")

    # Plot scoring method comparison
    scoring_corrs = [sm["corr_ssgsea"], sm["corr_meanz"], sm["corr_pca"]]
    scoring_names = ["ssGSEA", "Mean Z-score", "Pathway PCA"]
    plot_scoring_method_comparison(scoring_corrs, scoring_names)

    # Plot shared gene sensitivity
    plot_shared_gene_sensitivity(sg["original_corr"], sg["reduced_corr"])

    # Plot tissue adjustment comparison
    if ta is not None:
        plot_tissue_adjustment_comparison(ta["unadjusted_corr"], ta["adjusted_corr"])

    # Print summary
    print("\n  Sensitivity analysis summary:")
    print(f"    Shared gene removal: {sg['n_removed_genes']} genes removed, "
          f"correlation similarity = {sg['correlation_similarity']:.4f}")
    for pair in sm["pairwise_similarity"]:
        print(f"    Scoring: {pair['comparison']}: r = {pair['pearson_r']:.4f}")
    if ta is not None:
        print(f"    Tissue adjustment similarity: {ta['similarity']:.4f}")

    # =====================================================
    # STEP 6: Export Summary
    # =====================================================
    print("\n[6/7] Exporting Results Summary")
    print("-" * 40)

    summary = {
        "dataset": {
            "n_samples": len(metadata),
            "n_genes": expr_df.shape[1],
            "n_hallmarks": len(HALLMARK_GENES),
            "tissues": metadata["tissue"].value_counts().to_dict(),
            "age_range": f"{metadata['age'].min():.0f}-{metadata['age'].max():.0f}",
        },
        "shared_genes": {g: [HALLMARK_SHORT[h] for h in hs]
                         for g, hs in sorted(shared.items(), key=lambda x: -len(x[1]))[:15]},
        "network": {
            "correlation_edges": corr_metrics["n_edges"],
            "partial_corr_edges": pcorr_metrics["n_edges"],
            "density": corr_metrics["density"],
            "communities": corr_metrics.get("communities", []),
        },
        "ml_results": {
            "age_classification": {
                "rf_accuracy": f"{clf_results['rf_accuracy']:.3f}",
                "gb_accuracy": f"{clf_results['gb_accuracy']:.3f}",
                "lr_accuracy": f"{clf_results['lr_accuracy']:.3f}",
            },
            "biological_age": {
                "gb_r2": f"{age_results['gb_r2_cv']:.3f}",
                "gb_mae": f"{age_results['gb_mae']:.1f}",
                "en_r2": f"{age_results['en_r2']:.3f}",
            },
        },
        "key_findings": [],
    }

    # Key findings
    # Top hub hallmarks
    if "node_metrics" in corr_metrics:
        top_hub = corr_metrics["node_metrics"]["eigenvector"].idxmax()
        summary["key_findings"].append(
            f"Hub hallmark: {top_hub} (highest eigenvector centrality)"
        )

    # Strongest correlation pair
    flat_corr = []
    from itertools import combinations
    for h1, h2 in combinations(corr_matrix.columns, 2):
        flat_corr.append((h1, h2, corr_matrix.loc[h1, h2]))
    flat_corr.sort(key=lambda x: abs(x[2]), reverse=True)
    h1, h2, r = flat_corr[0]
    summary["key_findings"].append(
        f"Strongest correlation: {HALLMARK_SHORT[h1]}-{HALLMARK_SHORT[h2]} (r={r:.3f})"
    )

    # Most predictable hallmark
    top_pred = cross_pred.iloc[0]
    summary["key_findings"].append(
        f"Most predictable hallmark: {HALLMARK_SHORT[top_pred['target_hallmark']]} "
        f"(R²={top_pred['rf_r2_cv']:.3f})"
    )

    with open(RESULTS_DIR / "analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Provenance
    import platform
    provenance = {
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "seed": 42,
        "n_samples": 300,
        "n_genes": expr_df.shape[1],
        "pipeline_parameters": {
            "correlation_method": "spearman",
            "bonferroni_alpha": 0.05,
            "partial_corr_threshold": 0.15,
            "pc_alpha": 0.05,
            "ml_cv_folds": 5,
        },
    }
    try:
        import sklearn, scipy, networkx, matplotlib
        provenance["package_versions"] = {
            "scikit-learn": sklearn.__version__,
            "scipy": scipy.__version__,
            "networkx": networkx.__version__,
            "matplotlib": matplotlib.__version__,
        }
    except Exception:
        pass
    with open(RESULTS_DIR / "provenance.json", "w") as f:
        json.dump(provenance, f, indent=2)

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"\nResults directory: {RESULTS_DIR}")
    print(f"Figures: {RESULTS_DIR / 'figures'}")
    print(f"\nKey findings:")
    for finding in summary["key_findings"]:
        print(f"  - {finding}")

    return summary


if __name__ == "__main__":
    main()
