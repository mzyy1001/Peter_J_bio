# Review Round 8

Scope reviewed: `paper/main.tex`, `REVIEW_ROUND7.md`, `results/benchmark_summary.json`, `results/sensitivity_summary.csv`, `results/provenance.json`, `README.md`, plus targeted checks of `run_pipeline.py`, `src/data_acquisition.py`, `src/benchmarking.py`, `src/network_analysis.py`, `src/ml_models.py`, `src/sensitivity.py`, and `src/visualization.py`.

## Score

**Overall: 7/10**

Suggested sub-scores:

- **Code correctness: 7/10**
- **Statistical methods: 7/10**
- **Reproducibility: 6/10**

This is a real improvement over Round 7. The authors fixed the most serious benchmark-target problem, made shared-gene simulation consistent with multi-hallmark scoring, anchored PCA signs, added FDR columns for interaction rewiring, added README/provenance files, and updated the paper numbers. The paper is now much closer to a defensible simulation-framework manuscript.

However, I would not call it submission-ready yet. The current manuscript still over-interprets weak benchmark results, contains a few paper-code/simulation mismatches, and leaves important methodological weaknesses unresolved. The corrected benchmark is also less flattering than the previous one: **AUROC = 0.454**, **module ARI = -0.169**, and **edge stability = 17%** should be treated as central limitations, not as a mostly successful validation.

## What Was Fixed From Round 7

1. **Benchmark target mismatch: fixed in code.**  
   `src/data_acquisition.py` now saves both `ground_truth_base_hallmark_correlations.csv` and `ground_truth_observed_latent_correlations.csv`, and `run_pipeline.py` uses the observed latent correlation matrix when available. This directly addresses the biggest Round 7 issue.

2. **Shared-gene simulation: substantially fixed.**  
   Shared genes are now simulated from the mean of all assigned hallmark scores rather than only a first/primary hallmark. `gene_simulation_loadings.csv` is also saved, which improves auditability.

3. **PCA sensitivity sign anchoring: fixed.**  
   `src/sensitivity.py` now flips PC1 to align with mean z-score direction. The PCA similarity values are now interpretable.

4. **FDR correction for interaction strength model: mostly fixed.**  
   `interaction_strength_model()` now adds `pvalue_fdr` and `significant_fdr`, and the plotting function uses FDR-corrected significance if present.

5. **Some paper-code mismatches: fixed.**  
   Gradient boosting learning rate and MI wording are better aligned than in Round 7. Heatmap captions now say stars are unadjusted p-values.

6. **Reproducibility scaffolding: improved.**  
   `README.md`, `.gitignore`, and `results/provenance.json` now exist. This raises reproducibility from weak to adequate, though still not strong.

## Major Remaining Issues

### 1. The manuscript overstates edge recovery despite AUROC below chance

The benchmark now reports:

- edge AUROC = **0.454**
- best F1 = **0.764**
- centrality rho = **0.433**
- module ARI = **-0.169**
- stable edges = **16.7%**

The paper still frames this as "moderate edge detection" in the abstract and emphasizes the good F1. That is too generous. AUROC below 0.5 means the continuous edge scores rank true edges worse than false edges under the chosen ground truth. Best F1 is less persuasive because it is selected post hoc over thresholds and can be inflated when the observed-latent ground truth labels many of the 36 possible pairs as positive.

Recommended revision: describe edge recovery as **threshold-dependent and weak by ranking metrics**, not moderate. The abstract should say something like: "high best-threshold F1 but weak AUROC."

### 2. The paper still has a ground-truth wording mismatch

The code benchmarks against `ground_truth_observed_latent_correlations.csv`, which includes age and tissue effects. But Table 2's caption says true edges are defined by `|rho_embedded| > 0.1` in the latent correlation matrix. That sounds like the base pre-age/pre-tissue matrix.

This is smaller than the Round 7 bug because the code is now doing the better thing, but the paper must name the benchmark target precisely:

- base latent correlation matrix: pre age/tissue
- observed latent correlation matrix: post age/tissue, pre gene-expression/noise/scoring

The benchmark table should say it uses the **observed latent hallmark correlations after age and tissue effects**, if that is the intended target.

### 3. Random signed gene loadings remain inconsistent with unsigned scoring

The simulator still assigns random positive or negative gene loadings (`src/data_acquisition.py`, loading sampled with random sign), while `compute_hallmark_scores()` remains unsigned mean/rank scoring. This is now clearly visible in `results/analysis_summary.json`: the strongest recovered correlation is reported as **GI-MD r = -0.541**, while the latent target has a positive GI-MD association.

This is probably the main reason the corrected AUROC is poor. The current paper acknowledges missing directionality as a limitation, but the benchmark diagnosis should be sharper: the simulation intentionally creates up/down mixed signatures that the primary scoring method cannot recover.

Recommended fix: either implement a signed scoring benchmark using `gene_simulation_loadings.csv`, or explicitly state that the current primary score is an unsigned pathway activity proxy being tested under a difficult mixed-direction regime.

### 4. The PC/causal orientation section is still too strong

The methods caveat is better, but the Results still say the PC approximation identified orientations "consistent with the hierarchical framework." The implementation remains a depth-1 conditional-dependence heuristic with no separation-set tracking, no full PC orientation rules, and no multiple-testing correction across conditional independence tests.

This can stay as an exploratory visualization, but the Results should not imply validation of hierarchy or information flow. Rename the code/function and paper language away from "PC algorithm" unless a real PC implementation is used.

### 5. Some biological-category claims are plainly wrong

The abstract says "integrative hallmarks (MD, DNS) were most predictable" and "primary hallmarks (SCE) retained the most independent variation." Under the nine-hallmark framework used in the paper:

- MD and DNS are antagonistic hallmarks, not integrative.
- SCE is integrative, not primary.

The same conceptual error reappears in the conclusion, which claims the framework distinguishes integrative hallmarks from primary hallmarks via cross-hallmark prediction. The actual result is more modest: in this simulation, MD and DNS were most predictable, while SCE was least predictable. Do not map that cleanly onto the primary/antagonistic/integrative hierarchy.

### 6. The manuscript claims age-dependent rewiring was embedded, but the simulator does not appear to embed it

The Discussion says the simulation embedded differential correlation structure between age groups and the pipeline failed to detect it. In `src/data_acquisition.py`, I see age main effects and tissue offsets, but no age-varying hallmark covariance matrix or age-specific interaction terms. Therefore, the "failure" is not a power failure to detect embedded rewiring; it is a negative analysis in a simulation where rewiring was apparently not explicitly simulated.

Recommended revision: say no FDR-significant age-dependent correlation changes were detected, and avoid claiming that age-dependent rewiring was embedded unless the simulator is changed to include it.

### 7. ML claims remain a little too favorable

The paper says all age-group classifiers exceeded both baselines, but Gradient Boosting achieved **43.3%**, equal to the majority-class and stratified-random baselines reported in the same paragraph. This should be corrected.

The phrase "biological age estimation" is still somewhat overloaded. The paper now correctly says simulated residuals are not biological age acceleration, but the section and figure names still lean on biological-age framing. For a simulation-only paper, "chronological age prediction from hallmark scores" is cleaner.

### 8. Benchmark metrics are not fully methodologically clean

Several Round 7 statistical caveats remain:

- Best F1 is chosen on the evaluated seed rather than selected on training/simulation replicates and evaluated on held-out replicates.
- Centrality benchmark compares ground-truth weighted degree from a full matrix to estimated eigenvector centrality from a thresholded network.
- Module recovery benchmarks against the three hallmark categories, although the observed latent correlation matrix after age/tissue effects may not naturally preserve that partition.
- Robustness analysis repeatedly writes generated data files during each seed run, which is not ideal for provenance and can overwrite the main seed's intermediate data.

These are not fatal for a simulation-framework paper, but they should be disclosed or tightened before publication.

## New Issues Introduced or Exposed by the Fixes

1. **Observed-target benchmarking changed the class balance.**  
   Because age effects raise many off-diagonal observed correlations above 0.1, the "true edge" definition is now dense. This makes best F1 easier to interpret incorrectly. Report the number of positive/negative benchmark edges.

2. **The simulator now saves multiple ground-truth matrices, but the old `data/ground_truth_hallmark_correlations.csv` remains.**  
   This stale file could confuse users. Remove it or document that it is deprecated.

3. **Provenance is still incomplete.**  
   It records seed, sample count, some package versions, and thresholds, but not git commit, command-line invocation, GenAge download/fallback status, operating directory, complete dependency versions, or the exact ground-truth target used.

4. **README is minimal.**  
   It is enough to run the pipeline, but does not explain that all outputs are simulation-only, that running robustness overwrites data files, or how to reproduce the exact manuscript figures and PDF.

## Publishability Assessment

This is now publishable in principle as a **simulation-benchmarked software/framework note**, not as a biological discovery paper. The manuscript has become much more honest, and the code is coherent enough for readers to inspect and run.

Before submission, I would require these changes:

1. Reframe edge recovery as weak/mixed, not moderate, because AUROC is below chance and stability is low.
2. Correct the benchmark-target wording to specify observed latent correlations after age/tissue effects.
3. Fix the hallmark-category errors around MD, DNS, and SCE.
4. Remove or correct the claim that age-dependent rewiring was embedded.
5. Explicitly diagnose unsigned scoring vs random signed gene loadings as a major source of recovery failure.
6. Tone down the PC orientation Results section to exploratory conditional-dependence visualization only.
7. Add a short benchmark class-balance table or at least report true-positive edge count under the observed target.

With those revisions, I would score it **8/10** as a transparent simulation study. Without them, the current version is a **7/10**: substantially improved and technically interesting, but still too rhetorically strong relative to its own benchmark results.

