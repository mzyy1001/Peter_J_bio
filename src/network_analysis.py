"""
Network analysis module: constructs and analyzes hallmark interaction networks.

Approaches:
1. Gene overlap network (shared genes between hallmarks)
2. Correlation-based network (expression correlation between hallmark gene sets)
3. Protein-protein interaction (PPI) enriched network
4. Bayesian network structure learning
"""

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from itertools import combinations
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def build_gene_overlap_network(hallmark_genes):
    """
    Build a weighted network where edges represent shared genes between hallmarks.
    Weight = Jaccard similarity of gene sets.
    """
    G = nx.Graph()
    hallmarks = list(hallmark_genes.keys())
    for h in hallmarks:
        G.add_node(h, size=len(hallmark_genes[h]))

    for h1, h2 in combinations(hallmarks, 2):
        s1, s2 = set(hallmark_genes[h1]), set(hallmark_genes[h2])
        shared = s1 & s2
        if shared:
            jaccard = len(shared) / len(s1 | s2)
            G.add_edge(h1, h2, weight=jaccard, shared_genes=sorted(shared),
                       n_shared=len(shared))
    return G


def compute_hallmark_scores(expr_df, hallmark_genes, method="ssgsea"):
    """
    Compute per-sample hallmark activity scores using gene set scoring.

    Methods:
    - mean: simple mean of z-scored gene expression
    - ssgsea: simplified single-sample GSEA (rank-based)
    """
    # Z-score normalize expression
    expr_z = (expr_df - expr_df.mean()) / expr_df.std()

    hallmarks = list(hallmark_genes.keys())
    scores = pd.DataFrame(index=expr_df.index, columns=hallmarks, dtype=float)

    for h in hallmarks:
        genes = [g for g in hallmark_genes[h] if g in expr_z.columns]
        if not genes:
            scores[h] = 0.0
            continue

        if method == "mean":
            scores[h] = expr_z[genes].mean(axis=1)
        elif method == "ssgsea":
            # Simplified ssGSEA: rank genes, compute enrichment score
            ranks = expr_z.rank(axis=1)
            n_total = ranks.shape[1]
            n_set = len(genes)
            set_ranks = ranks[genes]
            # Weighted running sum approximation
            scores[h] = (set_ranks.sum(axis=1) / n_set - (n_total + 1) / 2) / n_total

    return scores.astype(float)


def build_correlation_network(hallmark_scores, method="spearman", threshold=0.1):
    """
    Build a correlation network between hallmark activity scores.
    """
    if method == "spearman":
        corr_matrix = hallmark_scores.corr(method="spearman")
    else:
        corr_matrix = hallmark_scores.corr(method="pearson")

    # Compute p-values
    n = len(hallmark_scores)
    hallmarks = hallmark_scores.columns.tolist()
    p_matrix = pd.DataFrame(np.ones((len(hallmarks), len(hallmarks))),
                            index=hallmarks, columns=hallmarks)

    for h1, h2 in combinations(hallmarks, 2):
        if method == "spearman":
            r, p = stats.spearmanr(hallmark_scores[h1], hallmark_scores[h2])
        else:
            r, p = stats.pearsonr(hallmark_scores[h1], hallmark_scores[h2])
        p_matrix.loc[h1, h2] = p
        p_matrix.loc[h2, h1] = p

    # Bonferroni correction: divide alpha by the number of pairwise tests
    n_tests = len(hallmarks) * (len(hallmarks) - 1) // 2
    bonferroni_alpha = 0.05 / n_tests if n_tests > 0 else 0.05

    G = nx.Graph()
    for h in hallmarks:
        G.add_node(h)

    for h1, h2 in combinations(hallmarks, 2):
        r = corr_matrix.loc[h1, h2]
        p = p_matrix.loc[h1, h2]
        if abs(r) > threshold and p < bonferroni_alpha:
            G.add_edge(h1, h2, weight=abs(r), correlation=r,
                       pvalue=p, significant=p < (0.001 / n_tests))

    return G, corr_matrix, p_matrix


def partial_correlation_network(hallmark_scores, threshold=0.15):
    """
    Build a partial correlation network to identify direct relationships
    (controlling for all other hallmarks).

    Uses Graphical Lasso with an aggressive alpha range to encourage sparsity,
    and a higher default threshold (0.15) to retain only meaningful edges.
    """
    from sklearn.covariance import GraphicalLassoCV

    # Standardize
    X = hallmark_scores.values
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    hallmarks = hallmark_scores.columns.tolist()

    # Graphical Lasso for sparse precision matrix
    # Use a higher alpha range to encourage genuine sparsity
    try:
        model = GraphicalLassoCV(
            alphas=np.logspace(-1, 1, 20),  # more aggressive regularization range
            cv=5, max_iter=1000
        )
        model.fit(X)
        precision = model.precision_
    except Exception:
        # Fallback: simple inverse covariance with stronger regularization
        cov = np.cov(X.T)
        precision = np.linalg.inv(cov + 0.1 * np.eye(cov.shape[0]))

    # Convert precision to partial correlations
    D = np.diag(1.0 / np.sqrt(np.diag(precision)))
    pcorr = -D @ precision @ D
    np.fill_diagonal(pcorr, 1.0)

    pcorr_df = pd.DataFrame(pcorr, index=hallmarks, columns=hallmarks)

    G = nx.Graph()
    for h in hallmarks:
        G.add_node(h)

    for h1, h2 in combinations(hallmarks, 2):
        i, j = hallmarks.index(h1), hallmarks.index(h2)
        r = pcorr[i, j]
        if abs(r) > threshold:
            G.add_edge(h1, h2, weight=abs(r), partial_corr=r)

    return G, pcorr_df


def compute_network_metrics(G):
    """Compute key network metrics."""
    metrics = {}
    metrics["n_nodes"] = G.number_of_nodes()
    metrics["n_edges"] = G.number_of_edges()
    metrics["density"] = nx.density(G)

    if G.number_of_edges() > 0:
        metrics["avg_clustering"] = nx.average_clustering(G, weight="weight")
        metrics["transitivity"] = nx.transitivity(G)

        # Node-level metrics
        degree = dict(G.degree(weight="weight"))
        betweenness = nx.betweenness_centrality(G, weight="weight")
        closeness = nx.closeness_centrality(G)
        if nx.is_connected(G):
            eigenvector = nx.eigenvector_centrality_numpy(G, weight="weight")
        else:
            eigenvector = nx.degree_centrality(G)

        node_df = pd.DataFrame({
            "weighted_degree": degree,
            "betweenness": betweenness,
            "closeness": closeness,
            "eigenvector": eigenvector,
        })
        metrics["node_metrics"] = node_df

        # Community detection
        try:
            from networkx.algorithms.community import greedy_modularity_communities
            communities = list(greedy_modularity_communities(G, weight="weight"))
            metrics["communities"] = [sorted(c) for c in communities]
            metrics["modularity"] = nx.algorithms.community.modularity(
                G, communities, weight="weight"
            )
        except Exception:
            metrics["communities"] = []
            metrics["modularity"] = 0.0

    return metrics


def _benjamini_hochberg(pvalues):
    """
    Apply Benjamini-Hochberg FDR correction to a 1-D array of p-values.
    Returns an array of adjusted p-values (same order as input).
    """
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    if n == 0:
        return pvalues
    order = np.argsort(pvalues)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    adjusted = pvalues * n / ranked
    # Enforce monotonicity (step-up procedure)
    adjusted = np.minimum(adjusted, 1.0)
    # Walk from largest rank down to enforce non-increasing adjusted p-values
    sorted_idx = np.argsort(ranked)[::-1]
    cummin = np.inf
    for idx in sorted_idx:
        cummin = min(cummin, adjusted[idx])
        adjusted[idx] = cummin
    return adjusted


def age_stratified_networks(hallmark_scores, metadata, age_bins=None,
                            fdr_alpha=0.05):
    """
    Build separate correlation networks for young, middle, and old groups.
    Test if network structure changes with age using Fisher z-tests
    with Benjamini-Hochberg FDR correction for multiple comparisons.
    """
    if age_bins is None:
        age_bins = [(20, 40, "Young"), (40, 60, "Middle"), (60, 90, "Old")]

    networks = {}
    for low, high, label in age_bins:
        mask = (metadata["age"] >= low) & (metadata["age"] < high)
        if mask.sum() < 10:
            continue
        scores_sub = hallmark_scores.loc[metadata.loc[mask, "sample_id"]]
        G, corr, pval = build_correlation_network(scores_sub, threshold=0.05)
        networks[label] = {"graph": G, "corr": corr, "pval": pval,
                           "n_samples": mask.sum()}

    # Fisher z-test for age-dependent rewiring between all pairs of age groups
    group_labels = list(networks.keys())
    hallmarks = hallmark_scores.columns.tolist()
    rewiring_results = []
    for g1, g2 in combinations(group_labels, 2):
        corr1 = networks[g1]["corr"]
        corr2 = networks[g2]["corr"]
        n1 = networks[g1]["n_samples"]
        n2 = networks[g2]["n_samples"]
        for h1, h2 in combinations(hallmarks, 2):
            r1 = corr1.loc[h1, h2]
            r2 = corr2.loc[h1, h2]
            # Fisher z-transform
            z1 = np.arctanh(np.clip(r1, -0.9999, 0.9999))
            z2 = np.arctanh(np.clip(r2, -0.9999, 0.9999))
            se = np.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))
            z_stat = (z1 - z2) / se
            p_val = 2.0 * stats.norm.sf(abs(z_stat))
            rewiring_results.append({
                "group1": g1, "group2": g2,
                "hallmark1": h1, "hallmark2": h2,
                "r1": r1, "r2": r2,
                "z_stat": z_stat, "pvalue": p_val,
            })

    if rewiring_results:
        rewiring_df = pd.DataFrame(rewiring_results)
        # Benjamini-Hochberg FDR correction
        rewiring_df["pvalue_fdr"] = _benjamini_hochberg(
            rewiring_df["pvalue"].values
        )
        rewiring_df["significant_fdr"] = rewiring_df["pvalue_fdr"] < fdr_alpha
    else:
        rewiring_df = pd.DataFrame()

    return networks, rewiring_df


def causal_inference_pc(hallmark_scores, alpha=0.05):
    """
    PC algorithm approximation for causal structure learning.
    Uses conditional independence tests to orient edges.
    """
    hallmarks = hallmark_scores.columns.tolist()
    n_h = len(hallmarks)
    n = len(hallmark_scores)

    # Start with complete undirected graph
    adj = np.ones((n_h, n_h), dtype=bool)
    np.fill_diagonal(adj, False)

    # Phase 1: Remove edges based on marginal independence
    for i, j in combinations(range(n_h), 2):
        r, p = stats.spearmanr(hallmark_scores.iloc[:, i], hallmark_scores.iloc[:, j])
        if p > alpha:
            adj[i, j] = adj[j, i] = False

    # Phase 2: Remove edges based on conditional independence (depth=1)
    for i, j in combinations(range(n_h), 2):
        if not adj[i, j]:
            continue
        for k in range(n_h):
            if k == i or k == j:
                continue
            if not (adj[i, k] or adj[j, k]):
                continue
            # Test i _||_ j | k using partial correlation
            x, y, z = (hallmark_scores.iloc[:, c].values for c in (i, j, k))
            # Residualize
            from numpy.polynomial.polynomial import polyfit, polyval
            res_x = x - np.polyval(np.polyfit(z, x, 1), z)
            res_y = y - np.polyval(np.polyfit(z, y, 1), z)
            r, p = stats.spearmanr(res_x, res_y)
            if p > alpha:
                adj[i, j] = adj[j, i] = False
                break

    # Build directed graph using v-structure detection
    G = nx.DiGraph()
    for h in hallmarks:
        G.add_node(h)

    # Orient v-structures: i -> k <- j if i-k-j and i not adj j
    oriented = set()
    for k in range(n_h):
        neighbors = [n for n in range(n_h) if adj[k, n]]
        for i, j in combinations(neighbors, 2):
            if not adj[i, j]:  # i and j not adjacent -> v-structure
                G.add_edge(hallmarks[i], hallmarks[k], type="directed")
                G.add_edge(hallmarks[j], hallmarks[k], type="directed")
                oriented.add((i, k))
                oriented.add((j, k))

    # Remaining edges as undirected (add both directions)
    for i, j in combinations(range(n_h), 2):
        if adj[i, j] and (i, j) not in oriented and (j, i) not in oriented:
            G.add_edge(hallmarks[i], hallmarks[j], type="undirected")

    return G, pd.DataFrame(adj, index=hallmarks, columns=hallmarks)


def bootstrap_correlation_ci(hallmark_scores, method="spearman",
                             n_bootstrap=2000, ci_level=0.95, seed=42):
    """
    Compute bootstrap confidence intervals for pairwise correlations
    between hallmark activity scores.

    Parameters
    ----------
    hallmark_scores : pd.DataFrame
        Samples x hallmarks matrix of activity scores.
    method : str
        'spearman' or 'pearson'.
    n_bootstrap : int
        Number of bootstrap resamples.
    ci_level : float
        Confidence level (e.g., 0.95 for 95% CI).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    ci_df : pd.DataFrame
        One row per hallmark pair with columns: hallmark1, hallmark2,
        observed_r, ci_lower, ci_upper, ci_width.
    """
    rng = np.random.RandomState(seed)
    hallmarks = hallmark_scores.columns.tolist()
    n_samples = len(hallmark_scores)
    alpha = 1.0 - ci_level
    pairs = list(combinations(hallmarks, 2))

    # Compute observed correlations
    corr_func = stats.spearmanr if method == "spearman" else stats.pearsonr
    observed = {}
    for h1, h2 in pairs:
        r, _ = corr_func(hallmark_scores[h1], hallmark_scores[h2])
        observed[(h1, h2)] = r

    # Bootstrap
    boot_corrs = {pair: [] for pair in pairs}
    for _ in range(n_bootstrap):
        idx = rng.choice(n_samples, size=n_samples, replace=True)
        sample = hallmark_scores.iloc[idx]
        for h1, h2 in pairs:
            r, _ = corr_func(sample[h1].values, sample[h2].values)
            boot_corrs[(h1, h2)].append(r)

    # Compute percentile CIs
    results = []
    for h1, h2 in pairs:
        boot = np.array(boot_corrs[(h1, h2)])
        lo = np.percentile(boot, 100 * alpha / 2)
        hi = np.percentile(boot, 100 * (1 - alpha / 2))
        results.append({
            "hallmark1": h1, "hallmark2": h2,
            "observed_r": observed[(h1, h2)],
            "ci_lower": lo, "ci_upper": hi,
            "ci_width": hi - lo,
        })
    return pd.DataFrame(results)


def network_null_model(hallmark_scores, method="spearman", threshold=0.1,
                       n_permutations=1000, seed=42):
    """
    Permutation-based null model for network metrics. Shuffles sample labels
    and recomputes the correlation network to generate an empirical null
    distribution for density, modularity, mean clustering, and transitivity.

    Parameters
    ----------
    hallmark_scores : pd.DataFrame
        Samples x hallmarks matrix of activity scores.
    method : str
        Correlation method for build_correlation_network.
    threshold : float
        Correlation magnitude threshold for edge inclusion.
    n_permutations : int
        Number of permutations.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    observed_metrics : dict
        Network metrics for the observed (unpermuted) data.
    null_distributions : dict
        Keys are metric names; values are arrays of length n_permutations.
    empirical_pvalues : dict
        Two-sided empirical p-values: fraction of permuted values at least
        as extreme as observed.
    """
    rng = np.random.RandomState(seed)

    # Observed network
    G_obs, _, _ = build_correlation_network(
        hallmark_scores, method=method, threshold=threshold
    )
    obs = compute_network_metrics(G_obs)
    observed_metrics = {
        "density": obs.get("density", 0.0),
        "modularity": obs.get("modularity", 0.0),
        "avg_clustering": obs.get("avg_clustering", 0.0),
        "transitivity": obs.get("transitivity", 0.0),
        "n_edges": obs.get("n_edges", 0),
    }

    null_distributions = {k: [] for k in observed_metrics}

    for _ in range(n_permutations):
        # Shuffle sample labels independently for each hallmark column
        shuffled = hallmark_scores.copy()
        for col in shuffled.columns:
            shuffled[col] = rng.permutation(shuffled[col].values)

        G_perm, _, _ = build_correlation_network(
            shuffled, method=method, threshold=threshold
        )
        perm_m = compute_network_metrics(G_perm)
        for k in null_distributions:
            null_distributions[k].append(perm_m.get(k, 0.0))

    # Convert to arrays
    for k in null_distributions:
        null_distributions[k] = np.array(null_distributions[k])

    # Empirical p-values (two-sided)
    empirical_pvalues = {}
    for k in observed_metrics:
        null = null_distributions[k]
        obs_val = observed_metrics[k]
        # fraction of permuted values as extreme or more extreme
        empirical_pvalues[k] = (np.sum(np.abs(null) >= abs(obs_val)) + 1) / (
            n_permutations + 1
        )

    return observed_metrics, null_distributions, empirical_pvalues
