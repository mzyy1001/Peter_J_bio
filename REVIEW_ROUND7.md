# Deep Code and Method Review - Round 7

Scope reviewed: `src/hallmark_genes.py`, `src/data_acquisition.py`, `src/network_analysis.py`, `src/ml_models.py`, `src/benchmarking.py`, `src/sensitivity.py`, `src/visualization.py`, `run_pipeline.py`, `paper/main.tex`, `results/analysis_summary.json`, and `results/benchmark_summary.json`.

## Overall Assessment

The manuscript is substantially more honest than earlier versions: it clearly states that all findings are simulation-based and reports weak benchmark dimensions instead of hiding them. The codebase is readable and mostly runnable, and several previous leakage/multiple-testing issues have been improved. However, the methodological core still has important correctness problems that affect the interpretation of the benchmark metrics.

The largest issue is that the benchmark compares recovered hallmark-score correlations against the *pre-age/pre-tissue latent correlation matrix*, while the expression data are generated after adding shared age effects and tissue offsets. This means the estimator is being evaluated against a target that is not the actual correlation structure present in the simulated observations. The reported edge AUROC of 0.675, module ARI of -0.052, and centrality rho of 0.45 are therefore useful diagnostics, but not clean estimates of recovery performance.

## Scores

**Code Correctness: 6/10**

The pipeline is coherent and modular, but several implementation choices make the simulation/benchmark internally inconsistent. The code is not broken in a trivial sense, but benchmark conclusions depend on target definitions that do not match the generated data.

**Statistical Methods: 6/10**

Cross-validation and Bonferroni/FDR corrections are improved. Remaining weaknesses are target mismatch in benchmarking, global preprocessing before ML CV, uncorrected nominal significance in the interaction figure path, heuristic causal orientation, and unclear separation between exploratory in-sample importance estimates and out-of-fold performance.

**Reproducibility: 5/10**

There is a `requirements.txt` and deterministic seeds in many places, but no root README, no lockfile or environment file, no test suite, no command-line configuration, no data manifest, and the pipeline overwrites generated data/results. GenAge download failure silently falls back to synthetic data, which is convenient but should be explicitly logged in machine-readable provenance.

## Major Findings

### 1. Benchmark ground truth does not match the simulated observed structure

In `src/data_acquisition.py`, the latent correlation matrix is created at lines 97-129 and sampled at lines 131-134. Age effects are then added at lines 136-150 and tissue offsets at lines 152-161. However, the saved ground-truth matrix at lines 199-201 is the original latent `corr_matrix`, not the correlation matrix of `hallmark_scores` after age and tissue effects.

`src/benchmarking.py` then defines true edges from this saved matrix at lines 34-43. The recovered `estimated_corr` in `run_pipeline.py` comes from `compute_hallmark_scores()` on gene expression at lines 88-112, which reflects latent correlation plus shared age trends, tissue offsets, gene-loading signs, and scoring artifacts.

This target mismatch can both penalize true recovery and reward confounded correlations. It also makes the module and centrality benchmarks hard to interpret.

Recommended code change:

```python
# src/data_acquisition.py
base_corr = pd.DataFrame(corr_matrix, index=hallmark_names, columns=hallmark_names)
base_corr.to_csv(DATA_DIR / "ground_truth_base_hallmark_correlations.csv")

observed_latent_corr = pd.DataFrame(
    hallmark_scores, columns=hallmark_names
).corr(method="spearman")
observed_latent_corr.to_csv(DATA_DIR / "ground_truth_observed_latent_correlations.csv")
```

Then decide explicitly in `run_pipeline.py` whether benchmarks target:

- **base latent biology**: residualize out age/tissue before estimating correlations, or
- **observed simulation structure**: benchmark against `ground_truth_observed_latent_correlations.csv`.

The manuscript should name which target is used.

### 2. Shared genes are simulated as single-primary-hallmark genes but scored as multi-hallmark genes

`src/data_acquisition.py` removes duplicate genes while keeping only the first hallmark label at lines 170-178. Gene expression is then generated from only that primary hallmark at lines 183-190. But `src/network_analysis.py` scores shared genes in every hallmark gene set at lines 56-71.

Example: `TP53` is scored in GI, CS, and SCE, but its expression is simulated from only the first assigned hallmark. This creates artificial asymmetric cross-talk and can distort the shared-gene sensitivity analysis.

Recommended code change: use `get_gene_hallmark_map()` during simulation and generate shared-gene expression from a weighted sum of all assigned hallmark activities:

```python
gene_map = get_gene_hallmark_map()
for gi, gene in enumerate(unique_genes):
    h_idxs = [hallmark_names.index(h) for h in gene_map[gene]]
    weights = rng.uniform(0.3, 1.0, size=len(h_idxs))
    weights = weights / weights.sum()
    latent = hallmark_scores[:, h_idxs] @ weights
    expression[:, gi] = base + loading * latent + noise
```

Also save the simulated gene-to-hallmark loadings to `data/gene_simulation_loadings.csv` for reproducibility.

### 3. Random signed gene loadings conflict with unsigned hallmark scoring

In `src/data_acquisition.py`, each gene receives a random positive or negative loading at line 187. But `compute_hallmark_scores()` in `src/network_analysis.py` uses either an unsigned mean z-score or unsigned rank-sum proxy at lines 62-71. This can invert hallmark scores depending on the random balance of gene signs.

The output summary confirms this problem: `results/analysis_summary.json` reports the strongest correlation as CS-SCE with `r=-0.572` at lines 111-114, despite the simulator embedding positive hallmark relationships.

Recommended code change:

- Store the true gene loading sign from the simulator.
- Add a direction-aware scoring option:

```python
def compute_hallmark_scores(expr_df, hallmark_genes, method="signed_mean", gene_weights=None):
    ...
    if method == "signed_mean":
        weights = pd.Series({g: gene_weights.get(g, 1.0) for g in genes})
        scores[h] = expr_z[genes].mul(weights, axis=1).mean(axis=1)
```

For real data, use curated up/down signatures or compare against established direction-aware methods such as singscore/GSVA/PLAGE.

### 4. The "ssGSEA" implementation is a rank-sum proxy, not standard ssGSEA

`src/network_analysis.py` implements `method="ssgsea"` as a normalized average rank at lines 64-71. This is closer to a simple single-sample rank-sum score than standard ssGSEA's weighted running-sum statistic.

The manuscript partially acknowledges this as "simplified ssGSEA" at `paper/main.tex` lines 134-143 and 474, which is good. But the abstract and methods still say "ssGSEA-based" at lines 46 and 134. That is acceptable only if consistently framed as a proxy.

Recommended manuscript/code alignment: rename the method to `rank_sum` or `ssgsea_proxy`, and reserve `ssgsea` for a validated implementation.

### 5. PCA scoring sensitivity is uninterpretable without sign anchoring

`src/sensitivity.py` computes pathway PCA scores at lines 130-150. PCA component signs are arbitrary, so the negative similarity values in `results/sensitivity_summary.csv` lines 4-6 and `paper/main.tex` lines 421-425 should not be interpreted as evidence that PCA gives a fundamentally different biological structure unless PC1 signs are anchored.

Recommended code change:

```python
pc1 = pca.fit_transform(subset.values)[:, 0]
mean_score = subset.mean(axis=1).values
if np.corrcoef(pc1, mean_score)[0, 1] < 0:
    pc1 = -pc1
scores[h] = pc1
```

Then rerun the scoring sensitivity analysis. This may substantially change the `ssGSEA vs pathway PCA` similarity.

### 6. Causal/PC implementation is too heuristic for the reported specificity

`src/network_analysis.py` labels `causal_inference_pc()` as a PC approximation at lines 291-351. It performs marginal Spearman tests, then depth-1 conditioning only, with OLS residualization and Spearman correlation at lines 310-327. It does not track separation sets, condition on larger sets, correct for multiple testing, or implement standard PC orientation rules beyond v-structures.

There is also a code/comment mismatch: lines 346-350 say remaining edges are added "as undirected (add both directions)", but the code adds only one directed edge in a `DiGraph` with `type="undirected"`.

The manuscript caveats are good at `paper/main.tex` lines 165-166, 299-306, and 453-455. Still, the result at lines 301-302 ("consistent with the hierarchical framework") is stronger than the implementation supports.

Recommended code change:

- Rename function to `conditional_dependence_orientation_heuristic()`.
- Save p-values and conditioning variables for each removed edge.
- Apply FDR correction or clearly report uncorrected exploratory tests.
- If claiming PC, use a tested library or implement separation-set tracking and full PC orientation rules.

### 7. Global feature construction creates mild ML leakage

`run_pipeline.py` computes hallmark scores once on all samples at line 88, using `compute_hallmark_scores()`, which z-scores expression across all samples at `src/network_analysis.py` line 51. ML models then cross-validate on these precomputed scores in `src/ml_models.py`.

The supervised target is not used in scoring, so this is not severe label leakage. But for strict predictive evaluation, test-fold expression distributions influence the training-fold feature normalization. This matters because the manuscript emphasizes leakage prevention at `paper/main.tex` lines 181-185.

Recommended code change: create a scikit-learn transformer for hallmark scoring and put it inside the CV pipeline when evaluating ML from expression-level inputs. If the intended ML task is "age prediction from already-derived cohort-level hallmark scores," then state that the CV evaluates model training only, not the full expression-to-score pipeline.

### 8. Interaction rewiring significance is inconsistent between code paths

`age_stratified_networks()` applies Benjamini-Hochberg FDR correction at `src/network_analysis.py` lines 278-284. But `interaction_strength_model()` uses nominal `p_diff < 0.05` at `src/ml_models.py` lines 489-502, and `plot_interaction_changes()` uses those nominal flags/stars at `src/visualization.py` lines 342-375.

The manuscript says no pairs survived FDR correction at `paper/main.tex` lines 366-374, but Figure 10 can still display nominal significance stars from the uncorrected `interaction_df`.

Recommended code change:

- Add BH-FDR columns to `interaction_strength_model()`.
- Use `q_value < 0.05` for stars in `plot_interaction_changes()`.
- Update the caption to specify whether stars are nominal p-values or FDR-adjusted q-values.

### 9. ML method details do not fully match the code

`paper/main.tex` line 182 says Gradient Boosting classification used learning rate 0.05. In `src/ml_models.py`, `GradientBoostingClassifier` is created at lines 82-85 without `learning_rate`, so scikit-learn's default 0.1 is used.

`paper/main.tex` lines 193-194 describe "pairwise mutual information." In `src/ml_models.py`, lines 463-471 compute mutual information by predicting each hallmark from all other hallmark columns and filling an asymmetric matrix. This is not pairwise symmetric MI.

Recommended code changes:

- Either set `learning_rate=0.05` in `GradientBoostingClassifier` or revise the manuscript to 0.1.
- For true pairwise MI, loop over hallmark pairs and call `mutual_info_regression(hallmark_scores[[h1]], hallmark_scores[h2])`, then symmetrize.

### 10. Reproducibility is incomplete

`requirements.txt` exists, but it uses broad lower bounds at lines 1-11 and there is no root `README.md`, lockfile, environment YAML, Dockerfile, CLI help, test suite, or CI. `run_pipeline.py` hard-codes `n_samples=300` and `seed=42` at line 74, and many analysis parameters are hard-coded throughout the pipeline.

`src/data_acquisition.py` silently switches to synthetic GenAge fallback on download failure at lines 28-43. That is convenient, but it means two users can get different provenance while both see a successful run.

Recommended changes:

- Add a root `README.md` with exact commands:
  - `python -m venv .venv`
  - `pip install -r requirements.txt`
  - `python run_pipeline.py`
- Add `environment.yml` or a pinned `requirements-lock.txt`.
- Add `--seed`, `--n-samples`, `--output-dir`, and `--skip-download` CLI options.
- Save `results/provenance.json` containing package versions, seed, command-line args, git commit, GenAge download status, and all key thresholds.
- Add smoke tests for `generate_simulated_expression_data()`, `compute_hallmark_scores()`, `build_correlation_network()`, and `run_all_benchmarks()`.

## Model Improvements for Current Benchmark Metrics

The current metrics are:

- edge AUROC = 0.675 (`results/benchmark_summary.json` line 2)
- module ARI = -0.052 (`results/benchmark_summary.json` line 6)
- centrality rho = 0.45 (`results/benchmark_summary.json` line 4)

Specific changes most likely to improve them:

1. **Fix the benchmark target first.** Save both base latent and observed latent correlation matrices in `src/data_acquisition.py` lines 199-206. Benchmark edge recovery against observed latent correlations if using unadjusted hallmark scores, or residualize age/tissue before benchmarking against base latent correlations. This may improve AUROC simply by evaluating against the correct target.

2. **Use signed scoring.** Replace unsigned rank-sum scoring in `src/network_analysis.py` lines 64-71 with signed mean/rank scoring using simulated gene loadings. This should reduce sign inversions such as CS-SCE `r=-0.572` and improve edge AUROC and centrality rho.

3. **Simulate shared genes from all assigned hallmarks.** Replace the primary-label-only logic in `src/data_acquisition.py` lines 170-190 with multi-hallmark loadings. This should make the score-generation and score-recovery models consistent and improve edge recovery.

4. **Benchmark multiple scoring methods directly.** Extend `run_all_benchmarks()` to accept score matrices from mean z-score, rank-sum, signed mean, anchored PCA, and latent-oracle scores. Report edge AUROC/module ARI/centrality rho by scoring method. This will identify whether the bottleneck is scoring, thresholding, or community detection.

5. **Community detection should use the full weighted matrix, not only the thresholded significance graph.** `compute_network_metrics()` runs greedy modularity on `corr_G`, whose edges are filtered in `build_correlation_network()` at lines 107-113. For module recovery, create a signed or absolute weighted graph from the full correlation matrix and benchmark communities across thresholds. This is likely to improve module ARI.

6. **Use centrality on the same object as the ground truth.** `benchmark_centrality_recovery()` defines ground-truth centrality as row sums of the full absolute correlation matrix at `src/benchmarking.py` lines 80-81, but estimated centrality uses eigenvector centrality from a thresholded network at lines 84-88. Compare weighted degree from the full estimated absolute correlation matrix to the full ground-truth weighted degree.

7. **Tune thresholds on simulation replicates, not on the evaluated seed.** `benchmark_edge_recovery()` reports best F1 across thresholds at `src/benchmarking.py` lines 48-70. Use seeds 0-8 to choose thresholds and seed 9 or held-out seeds to evaluate. This would make threshold recommendations reproducible and less circular.

8. **Increase sample size or signal-to-noise for power diagnostics.** The current `n=300` and random gene loadings produce weak age-effect recovery (`results/benchmark_summary.json` line 7). Add simulation sweeps over `n_samples`, `noise_sd`, and loading distributions. This will clarify whether metrics are limited by sample size or scoring mismatch.

## Paper-Code Alignment

Mostly aligned, but these claims need revision or code changes:

- `paper/main.tex` lines 89-96 claim curation integrates GenAge, MSigDB, and literature. `src/hallmark_genes.py` is a hard-coded dictionary at lines 6-68, and `download_genage_genes()` in `src/data_acquisition.py` lines 22-43 is not used to build or audit those sets. Add a provenance table or soften the claim to "manually curated."

- `paper/main.tex` line 182 states Gradient Boosting classifier learning rate 0.05, but `src/ml_models.py` lines 82-85 use the default 0.1.

- `paper/main.tex` lines 193-194 state pairwise mutual information, but `src/ml_models.py` lines 463-471 computes asymmetric one-vs-rest MI.

- `paper/main.tex` lines 245-246 heatmap caption marks `*p<0.05`, `**p<0.01`, `***p<0.001`, but `plot_hallmark_correlation_heatmap()` uses raw p-values from `src/visualization.py` lines 66-80, while the network uses Bonferroni filtering in `src/network_analysis.py` lines 99-113. Caption should say stars are unadjusted p-values, or the plot should use adjusted p-values.

- `paper/main.tex` line 77 says "Demonstration of the pipeline's ability to detect age-dependent network rewiring." Given no FDR-significant pairs at lines 366-374, revise to "evaluate age-dependent network rewiring."

- `paper/main.tex` lines 430-434 infer that shared-gene removal and tissue adjustment show correlations are not driven by overlap or tissue composition. In this simulation that is plausible, but because the simulator uses primary-only shared genes and modest tissue offsets, this should be phrased as "under this simulation regime."

- `paper/main.tex` line 495 says all code is open-source. The repository may be available locally, but the manuscript should include a URL, commit hash, and archived release/DOI for reproducibility.

## Final Recommendation

The paper remains acceptable as a simulation-benchmarked framework paper if the authors explicitly fix or disclose the benchmark target mismatch. Without that fix, the reported benchmark metrics are not wrong as descriptive outputs, but they are not measuring the clean recovery problem the manuscript says they measure.

Priority before submission:

1. Save and benchmark against the correct ground-truth correlation target.
2. Implement signed/direction-aware scoring or explain why random signed loadings are intentionally unrecoverable by unsigned scores.
3. Anchor PCA sensitivity signs and rerun sensitivity results.
4. Correct the manuscript-code mismatches for Gradient Boosting learning rate, mutual information, heatmap p-value stars, and rewiring language.
5. Add minimal reproducibility documentation and provenance output.
