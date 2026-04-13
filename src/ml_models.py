"""
Machine learning models for hallmark interaction analysis.

Models:
1. Random Forest classifier for aging stage prediction from hallmark scores
2. Elastic Net for hallmark-hallmark predictive relationships
3. Multi-output regression for cross-hallmark prediction
4. Gradient boosting for biological age estimation
5. Autoencoder for hallmark latent space
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import (
    cross_val_score, cross_val_predict, StratifiedKFold, KFold
)
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import Pipeline
from sklearn.ensemble import (
    RandomForestClassifier, RandomForestRegressor,
    GradientBoostingRegressor, GradientBoostingClassifier
)
from sklearn.linear_model import ElasticNet, ElasticNetCV, LogisticRegression
from sklearn.multioutput import MultiOutputRegressor
from sklearn.metrics import (
    accuracy_score, f1_score, r2_score, mean_absolute_error,
    classification_report, confusion_matrix, roc_auc_score
)
from sklearn.inspection import permutation_importance
from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

RESULTS_DIR = Path(__file__).parent.parent / "results"


def _bootstrap_ci(y_true, y_pred, metric_fn, n_boot=1000, alpha=0.05,
                  random_state=42):
    """Compute bootstrap confidence interval for a metric."""
    rng = np.random.RandomState(random_state)
    n = len(y_true)
    scores = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.randint(0, n, size=n)
        scores[i] = metric_fn(y_true[idx], y_pred[idx])
    lo = np.percentile(scores, 100 * alpha / 2)
    hi = np.percentile(scores, 100 * (1 - alpha / 2))
    return lo, hi


def age_group_classifier(hallmark_scores, metadata, n_groups=3):
    """
    Train a classifier to predict age group from hallmark activity scores.
    Returns feature importances showing which hallmarks are most predictive.
    """
    # Create age groups
    ages = metadata.set_index("sample_id").loc[hallmark_scores.index, "age"]
    if n_groups == 3:
        bins = [0, 40, 60, 100]
        labels_str = ["Young", "Middle", "Old"]
    else:
        bins = [0, 35, 50, 65, 100]
        labels_str = ["Young", "Early-Mid", "Late-Mid", "Old"]

    age_groups = pd.cut(ages, bins=bins, labels=labels_str)
    le = LabelEncoder()
    y = le.fit_transform(age_groups)

    X = hallmark_scores.values
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # --- Pipelines (scaler inside CV to prevent data leakage) ---
    rf_pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", RandomForestClassifier(n_estimators=200, max_depth=8,
                                       random_state=42)),
    ])
    gb_pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", GradientBoostingClassifier(n_estimators=150, max_depth=4,
                                           random_state=42)),
    ])
    lr_pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=1000, C=1.0)),
    ])

    # Cross-validated scores
    rf_scores = cross_val_score(rf_pipe, X, y, cv=cv, scoring="accuracy")
    gb_scores = cross_val_score(gb_pipe, X, y, cv=cv, scoring="accuracy")
    lr_scores = cross_val_score(lr_pipe, X, y, cv=cv, scoring="accuracy")

    # --- Baselines ---
    # Majority class baseline
    dummy_majority = DummyClassifier(strategy="most_frequent")
    majority_scores = cross_val_score(dummy_majority, X, y, cv=cv,
                                      scoring="accuracy")

    # Stratified random baseline
    dummy_strat = DummyClassifier(strategy="stratified", random_state=42)
    stratified_scores = cross_val_score(dummy_strat, X, y, cv=cv,
                                        scoring="accuracy")

    # Tissue-only baseline (if tissue column available)
    tissue_baseline_acc = None
    if "tissue" in metadata.columns:
        tissues = metadata.set_index("sample_id").loc[
            hallmark_scores.index, "tissue"
        ]
        tissue_le = LabelEncoder()
        X_tissue = tissue_le.fit_transform(tissues).reshape(-1, 1)
        tissue_pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=1000)),
        ])
        tissue_baseline_scores = cross_val_score(tissue_pipe, X_tissue, y,
                                                 cv=cv, scoring="accuracy")
        tissue_baseline_acc = tissue_baseline_scores.mean()

    # Fit final models for feature importance (on all data, clearly labelled)
    rf_pipe.fit(X, y)
    gb_pipe.fit(X, y)
    lr_pipe.fit(X, y)

    feature_names = hallmark_scores.columns.tolist()
    importances = pd.DataFrame({
        "hallmark": feature_names,
        "rf_importance": rf_pipe.named_steps["clf"].feature_importances_,
        "gb_importance": gb_pipe.named_steps["clf"].feature_importances_,
    }).sort_values("rf_importance", ascending=False)

    # Permutation importance (on full data, but measured on fitted pipeline)
    X_scaled_full = rf_pipe.named_steps["scaler"].transform(X)
    perm_imp = permutation_importance(
        rf_pipe.named_steps["clf"], X_scaled_full, y,
        n_repeats=30, random_state=42, scoring="accuracy"
    )
    importances["rf_perm_importance_mean"] = perm_imp.importances_mean
    importances["rf_perm_importance_std"] = perm_imp.importances_std

    # --- Confidence intervals via bootstrap on CV predictions ---
    y_pred_rf = cross_val_predict(rf_pipe, X, y, cv=cv)
    rf_acc_ci = _bootstrap_ci(
        y, y_pred_rf,
        lambda yt, yp: accuracy_score(yt, yp)
    )

    y_pred_gb = cross_val_predict(gb_pipe, X, y, cv=cv)
    gb_acc_ci = _bootstrap_ci(
        y, y_pred_gb,
        lambda yt, yp: accuracy_score(yt, yp)
    )

    results = {
        "rf_accuracy": rf_scores.mean(),
        "rf_std": rf_scores.std(),
        "rf_accuracy_ci95": rf_acc_ci,
        "gb_accuracy": gb_scores.mean(),
        "gb_std": gb_scores.std(),
        "gb_accuracy_ci95": gb_acc_ci,
        "lr_accuracy": lr_scores.mean(),
        "lr_std": lr_scores.std(),
        "feature_importances": importances,
        "classes": le.classes_,
        "rf_model": rf_pipe,
        "gb_model": gb_pipe,
        # Baselines
        "baseline_majority_accuracy": majority_scores.mean(),
        "baseline_stratified_accuracy": stratified_scores.mean(),
        "baseline_tissue_accuracy": tissue_baseline_acc,
    }

    # Confusion matrix from RF cross-validated predictions
    results["confusion_matrix"] = confusion_matrix(y, y_pred_rf)
    results["classification_report"] = classification_report(
        y, y_pred_rf, target_names=le.classes_, output_dict=True
    )

    return results


def biological_age_estimator(hallmark_scores, metadata):
    """
    Predict chronological age from hallmark scores (biological age model).
    The residual (predicted - actual) represents biological age acceleration.
    """
    ages = metadata.set_index("sample_id").loc[hallmark_scores.index, "age"].values
    X = hallmark_scores.values
    cv = KFold(n_splits=5, shuffle=True, random_state=42)

    # --- Pipelines (scaler inside CV to prevent data leakage) ---
    en_pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("reg", ElasticNetCV(l1_ratio=[0.1, 0.3, 0.5, 0.7, 0.9], cv=5,
                             max_iter=5000, random_state=42)),
    ])
    gb_pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("reg", GradientBoostingRegressor(n_estimators=200, max_depth=4,
                                          learning_rate=0.05,
                                          random_state=42)),
    ])
    rf_pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("reg", RandomForestRegressor(n_estimators=200, max_depth=8,
                                      random_state=42)),
    ])

    # --- Cross-validated predictions (out-of-fold) ---
    en_pred_cv = cross_val_predict(en_pipe, X, ages, cv=cv)
    gb_pred_cv = cross_val_predict(gb_pipe, X, ages, cv=cv)
    rf_pred_cv = cross_val_predict(rf_pipe, X, ages, cv=cv)

    # Cross-validated R^2 scores (per-fold)
    gb_scores = cross_val_score(gb_pipe, X, ages, cv=cv, scoring="r2")
    rf_scores = cross_val_score(rf_pipe, X, ages, cv=cv, scoring="r2")

    # --- Baselines ---
    # Predict-mean baseline
    dummy_mean = DummyRegressor(strategy="mean")
    mean_pred_cv = cross_val_predict(dummy_mean, X, ages, cv=cv)
    baseline_mean_mae = mean_absolute_error(ages, mean_pred_cv)
    baseline_mean_r2 = r2_score(ages, mean_pred_cv)

    # Tissue-only baseline (predict mean age per tissue)
    baseline_tissue_mae = None
    baseline_tissue_r2 = None
    if "tissue" in metadata.columns:
        tissues = metadata.set_index("sample_id").loc[
            hallmark_scores.index, "tissue"
        ]
        tissue_le = LabelEncoder()
        X_tissue = tissue_le.fit_transform(tissues).reshape(-1, 1)
        tissue_pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("reg", GradientBoostingRegressor(
                n_estimators=50, max_depth=2, random_state=42)),
        ])
        tissue_pred_cv = cross_val_predict(tissue_pipe, X_tissue, ages, cv=cv)
        baseline_tissue_mae = mean_absolute_error(ages, tissue_pred_cv)
        baseline_tissue_r2 = r2_score(ages, tissue_pred_cv)

    # Fit final models on all data (for coefficients / importances)
    en_pipe.fit(X, ages)
    gb_pipe.fit(X, ages)
    rf_pipe.fit(X, ages)

    # Age acceleration from cross-validated GB predictions
    age_accel = gb_pred_cv - ages

    feature_names = hallmark_scores.columns.tolist()
    en_model = en_pipe.named_steps["reg"]
    en_coefs = pd.DataFrame({
        "hallmark": feature_names,
        "coefficient": en_model.coef_,
        "abs_coef": np.abs(en_model.coef_),
    }).sort_values("abs_coef", ascending=False)

    gb_model = gb_pipe.named_steps["reg"]
    gb_importances = pd.DataFrame({
        "hallmark": feature_names,
        "importance": gb_model.feature_importances_,
    }).sort_values("importance", ascending=False)

    # Permutation importance for GB
    X_scaled_full = gb_pipe.named_steps["scaler"].transform(X)
    perm_imp = permutation_importance(
        gb_model, X_scaled_full, ages,
        n_repeats=30, random_state=42, scoring="r2"
    )
    gb_importances["perm_importance_mean"] = perm_imp.importances_mean
    gb_importances["perm_importance_std"] = perm_imp.importances_std

    # --- Confidence intervals via bootstrap on CV predictions ---
    en_r2_ci = _bootstrap_ci(
        ages, en_pred_cv, lambda yt, yp: r2_score(yt, yp)
    )
    gb_r2_ci = _bootstrap_ci(
        ages, gb_pred_cv, lambda yt, yp: r2_score(yt, yp)
    )
    gb_mae_ci = _bootstrap_ci(
        ages, gb_pred_cv, lambda yt, yp: mean_absolute_error(yt, yp)
    )

    results = {
        "en_r2": r2_score(ages, en_pred_cv),
        "en_r2_ci95": en_r2_ci,
        "en_mae": mean_absolute_error(ages, en_pred_cv),
        "en_coefficients": en_coefs,
        "gb_r2_cv": gb_scores.mean(),
        "gb_r2_cv_std": gb_scores.std(),
        "gb_r2_ci95": gb_r2_ci,
        "gb_mae": mean_absolute_error(ages, gb_pred_cv),
        "gb_mae_ci95": gb_mae_ci,
        "gb_importances": gb_importances,
        "rf_r2_cv": rf_scores.mean(),
        "rf_r2_cv_std": rf_scores.std(),
        "age_acceleration": pd.DataFrame({
            "sample_id": hallmark_scores.index,
            "chronological_age": ages,
            "predicted_age": gb_pred_cv,
            "age_acceleration": age_accel,
        }),
        "actual_ages": ages,
        "predicted_ages": gb_pred_cv,
        # Baselines
        "baseline_mean_mae": baseline_mean_mae,
        "baseline_mean_r2": baseline_mean_r2,
        "baseline_tissue_mae": baseline_tissue_mae,
        "baseline_tissue_r2": baseline_tissue_r2,
    }
    return results


def cross_hallmark_prediction(hallmark_scores):
    """
    For each hallmark, predict its score from all other hallmarks.
    High R2 indicates a hallmark is strongly determined by others.
    Uses leave-one-hallmark-out approach.
    """
    hallmarks = hallmark_scores.columns.tolist()
    results = []

    for target_h in hallmarks:
        X = hallmark_scores.drop(columns=[target_h]).values
        y = hallmark_scores[target_h].values

        # ElasticNet via Pipeline (scaler inside CV)
        en_pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("reg", ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], cv=5,
                                 max_iter=3000)),
        ])
        cv = KFold(n_splits=5, shuffle=True, random_state=42)
        en_pred_cv = cross_val_predict(en_pipe, X, y, cv=cv)
        en_r2 = r2_score(y, en_pred_cv)

        # Fit on all data for predictor weights
        en_pipe.fit(X, y)
        other_h = [h for h in hallmarks if h != target_h]
        predictor_weights = dict(zip(other_h,
                                     en_pipe.named_steps["reg"].coef_))

        # RF cross-validated via Pipeline
        rf_pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("reg", RandomForestRegressor(n_estimators=100, max_depth=6,
                                          random_state=42)),
        ])
        rf_r2 = cross_val_score(rf_pipe, X, y, cv=cv, scoring="r2").mean()

        results.append({
            "target_hallmark": target_h,
            "en_r2": en_r2,
            "rf_r2_cv": rf_r2,
            "top_predictor": max(predictor_weights,
                                 key=lambda k: abs(predictor_weights[k])),
            "top_predictor_weight": max(abs(v)
                                        for v in predictor_weights.values()),
            "predictor_weights": predictor_weights,
        })

    return pd.DataFrame(results).sort_values("rf_r2_cv", ascending=False)


def hallmark_pca_analysis(hallmark_scores, metadata):
    """
    PCA and t-SNE analysis of hallmark score space.
    """
    X = hallmark_scores.values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # PCA
    pca = PCA()
    X_pca = pca.fit_transform(X_scaled)

    pca_results = {
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumulative_variance": np.cumsum(pca.explained_variance_ratio_),
        "components": pd.DataFrame(
            pca.components_,
            columns=hallmark_scores.columns,
            index=[f"PC{i+1}" for i in range(len(pca.components_))]
        ),
        "scores": pd.DataFrame(
            X_pca[:, :3],
            columns=["PC1", "PC2", "PC3"],
            index=hallmark_scores.index,
        ),
    }

    # Add metadata
    ages = metadata.set_index("sample_id").loc[hallmark_scores.index, "age"]
    pca_results["scores"]["age"] = ages.values
    tissues = metadata.set_index("sample_id").loc[hallmark_scores.index, "tissue"]
    pca_results["scores"]["tissue"] = tissues.values

    # t-SNE
    tsne = TSNE(n_components=2, perplexity=30, random_state=42)
    X_tsne = tsne.fit_transform(X_scaled)
    pca_results["tsne"] = pd.DataFrame(
        X_tsne, columns=["tSNE1", "tSNE2"], index=hallmark_scores.index
    )
    pca_results["tsne"]["age"] = ages.values
    pca_results["tsne"]["tissue"] = tissues.values

    return pca_results


def tissue_specific_analysis(hallmark_scores, metadata):
    """
    Analyze tissue-specific hallmark patterns and interactions.
    """
    tissues = metadata.set_index("sample_id").loc[hallmark_scores.index, "tissue"]
    results = {}

    for tissue in tissues.unique():
        mask = tissues == tissue
        scores_t = hallmark_scores[mask]
        if len(scores_t) < 15:
            continue

        # Mean scores
        mean_scores = scores_t.mean()

        # Correlation within tissue
        corr = scores_t.corr(method="spearman")

        # Age correlation per hallmark
        ages_t = metadata.set_index("sample_id").loc[scores_t.index, "age"]
        age_corrs = {}
        for h in hallmark_scores.columns:
            r, p = stats.spearmanr(ages_t, scores_t[h])
            age_corrs[h] = {"r": r, "p": p}

        results[tissue] = {
            "n_samples": len(scores_t),
            "mean_scores": mean_scores,
            "correlation": corr,
            "age_correlations": pd.DataFrame(age_corrs).T,
        }

    return results


def interaction_strength_model(hallmark_scores, metadata):
    """
    Quantify pairwise interaction strengths between hallmarks using
    mutual information and maximal information coefficient approximation.
    """
    from sklearn.feature_selection import mutual_info_regression
    from itertools import combinations

    hallmarks = hallmark_scores.columns.tolist()
    n_h = len(hallmarks)

    # Mutual information matrix
    mi_matrix = np.zeros((n_h, n_h))
    for i in range(n_h):
        X = hallmark_scores.drop(columns=[hallmarks[i]]).values
        y = hallmark_scores[hallmarks[i]].values
        mi_vals = mutual_info_regression(X, y, random_state=42)
        for j_idx, j in enumerate([k for k in range(n_h) if k != i]):
            mi_matrix[i, j] = mi_vals[j_idx]

    mi_df = pd.DataFrame(mi_matrix, index=hallmarks, columns=hallmarks)

    # Age-conditioned interaction changes
    ages = metadata.set_index("sample_id").loc[hallmark_scores.index, "age"].values
    young_mask = ages < 40
    old_mask = ages >= 60

    interactions = []
    for h1, h2 in combinations(hallmarks, 2):
        if young_mask.sum() >= 10 and old_mask.sum() >= 10:
            r_young, _ = stats.spearmanr(
                hallmark_scores.loc[hallmark_scores.index[young_mask], h1],
                hallmark_scores.loc[hallmark_scores.index[young_mask], h2]
            )
            r_old, _ = stats.spearmanr(
                hallmark_scores.loc[hallmark_scores.index[old_mask], h1],
                hallmark_scores.loc[hallmark_scores.index[old_mask], h2]
            )
            # Fisher z-test for correlation difference
            z_young = np.arctanh(r_young)
            z_old = np.arctanh(r_old)
            n_y, n_o = young_mask.sum(), old_mask.sum()
            z_diff = (z_old - z_young) / np.sqrt(1/(n_y-3) + 1/(n_o-3))
            p_diff = 2 * (1 - stats.norm.cdf(abs(z_diff)))

            interactions.append({
                "hallmark_1": h1, "hallmark_2": h2,
                "corr_young": r_young, "corr_old": r_old,
                "corr_change": r_old - r_young,
                "z_statistic": z_diff, "p_value": p_diff,
                "significant": p_diff < 0.05,
            })

    interaction_df = pd.DataFrame(interactions)

    # Apply BH-FDR correction to p-values
    if len(interaction_df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, pvals_fdr, _, _ = multipletests(interaction_df["p_value"].values,
                                           method="fdr_bh")
        interaction_df["pvalue_fdr"] = pvals_fdr
        interaction_df["significant_fdr"] = pvals_fdr < 0.05

    return mi_df, interaction_df
