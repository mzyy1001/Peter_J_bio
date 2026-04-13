# Review Round 2

## 1. Overall Assessment

**Recommendation: Major Revision**  
**Score: 5/10**

The revision is a substantial improvement over Round 1. The authors have now explicitly reframed the work as a simulation-benchmarked computational framework, state that all results come from simulated data, remove much of the original biological-discovery framing, add ground-truth benchmarking code, correct several machine-learning evaluation problems, add baselines and confidence intervals, make the Graphical Lasso implementation sparser, and add FDR correction for age-rewiring tests. These changes directly address many of the most serious Round 1 concerns.

However, the revised manuscript is not yet acceptable. The main remaining problem is no longer hidden simulation, but rather internal inconsistency and overstatement of benchmark success. The manuscript repeatedly claims strong or reliable recovery, but `results/benchmark_summary.json` reports only modest edge recovery (AUROC = 0.675, best F1 = 0.653), weak centrality recovery (rho = 0.45), poor top-hub overlap (0.333), essentially failed module recovery (ARI = -0.052), weak age-effect recovery (rho = 0.201), and low edge stability (22.2%). These values do not support claims of "high fidelity," "faithful" recovery, "consistently detected" modules, or robust reliability.

There are also major manuscript/result mismatches. For example, `results/analysis_summary.json` reports biological-age MAE = 15.9 years and GB R2 = 0.103, while the manuscript still claims MAE approximately 5 years. The same JSON reports 17 correlation edges, partial-correlation edges = 12, and density = 0.472, while the manuscript reports 22 correlation edges, partial-correlation edges = 36, and density = 0.611. The FDR-corrected rewiring file reports no FDR-significant rewiring edges, while the manuscript still discusses four altered pairs. These are not minor polish issues; they undermine trust in the revision and must be fixed before the paper can be judged on its scientific merits.

In short, the authors have moved the paper into a potentially viable category: a simulation recovery and benchmarking study. But the current revision still overclaims, uses placeholders instead of reporting benchmark values in the main table, and contains stale or contradictory results. I would not reject outright because the trajectory is positive and many Round 1 concerns were genuinely addressed, but the manuscript requires another major revision.

## 2. Round 1 Issues Adequately Addressed

1. **Simulation framing and transparency.**  
   This is the largest improvement. The title, abstract, Methods, and Discussion now explicitly identify the study as simulation-benchmarked. The Methods include a "Critical distinction" paragraph stating that all reported results are from simulated data and are not independent biological discoveries. This directly addresses the Round 1 fatal concern that simulated data were presented as empirical multi-tissue transcriptomics.

2. **Repositioning as method validation rather than empirical discovery.**  
   The paper now describes the analysis as a framework validation step before application to real cohorts. The Discussion repeatedly states that tissue-specific signatures, hub findings, and network structure reflect embedded simulation design. This is a major improvement over the original version.

3. **Ground-truth benchmarking added.**  
   The code now includes benchmarking for edge recovery, centrality recovery, module recovery, age-effect recovery, and robustness across seeds. This is exactly the type of analysis requested in Round 1. The problem is that the benchmark outcomes are weak and not honestly reflected in the manuscript, but the analysis category itself has been added.

4. **Machine-learning leakage substantially improved.**  
   `src/ml_models.py` now uses scikit-learn `Pipeline` objects with scaling inside cross-validation and computes age-prediction metrics from out-of-fold predictions via `cross_val_predict`. This addresses the Round 1 concerns about in-sample biological-age MAE and preprocessing leakage.

5. **Baselines added.**  
   The code now includes majority-class and stratified-random baselines for classification, and mean-age and tissue-only baselines for regression. This is a meaningful improvement over comparison against chance only.

6. **Bootstrap confidence intervals added in code.**  
   The ML code includes bootstrap confidence intervals for classification accuracy and biological-age metrics, and the network code includes bootstrap correlation CIs. The current manuscript should report these values more fully, but the implementation direction is appropriate.

7. **Graphical Lasso implementation improved.**  
   The code now uses a more aggressive alpha range and a thresholded partial-correlation graph. `analysis_summary.json` reports 12 partial-correlation edges, which is much more plausible than the complete 36-edge "sparse" network criticized in Round 1.

8. **Multiple-testing correction improved in code.**  
   The correlation-network code now applies Bonferroni correction by using `0.05 / n_tests`, and the age-rewiring code now computes Benjamini-Hochberg adjusted p-values. This addresses two concrete Round 1 technical issues, although the manuscript still does not consistently reflect the corrected results.

9. **Causal language reduced in the Methods and Discussion.**  
   The Methods now call the PC result "Conditional-Dependence Orientation" and explicitly state that the implementation is approximate, conditions only on single variables, and should not be interpreted as causal inference. This is a substantial improvement.

10. **Limitations section improved.**  
    The revised Limitations section now acknowledges simulated data only, simplified ssGSEA, informal gene-set curation, shared-gene sensitivity concerns, and the approximate nature of the PC algorithm.

## 3. Round 1 Issues Remaining Unaddressed

1. **Overclaiming remains, now in benchmarking language rather than biological-discovery language.**  
   The original overstatement has not fully disappeared; it has shifted from "we discovered biological aging architecture" to "the framework reliably recovers known structure." The reported benchmark values do not justify this. Module ARI = -0.052 indicates worse-than-random or no meaningful agreement with the embedded module labels. Edge stability = 22.2% is not robust. Centrality rho = 0.45 is modest. These should be described as mixed or limited recovery, not high-fidelity validation.

2. **Biological interpretations are still too strong in several Results subsections.**  
   The tissue-specific section still says brain SCE suggests pronounced stem cell exhaustion, blood AIC reflects inflammaging, muscle SCE is consistent with sarcopenia, and liver TA reflects hepatocyte telomere attrition. The Discussion later qualifies these as simulated patterns, but the Results wording still reads as biological interpretation. For simulated data, the Results should say these profiles recover embedded tissue offsets.

3. **Causal framing persists in figures and wording.**  
   Although the section title has improved, Figure 11 is still named/captioned as "Estimated causal structure" and the caption says "putative causal relationships." The visualization code also labels directed edges as "Directed (causal)." This should be changed throughout to "oriented conditional-dependence graph" or similar.

4. **Gene-set curation remains under-specified.**  
   The manuscript still lacks a reproducible curation protocol, versioned source accessions, inclusion/exclusion criteria, and independent curation procedure. It acknowledges this as a limitation but does not fix it.

5. **Simplified ssGSEA remains unvalidated.**  
   The authors acknowledge the scoring method is simplified, but no comparison to GSVA/ssGSEA, singscore, PLAGE, mean z-score, or pathway PCA is performed. For a framework paper, this is still a substantive gap.

6. **No real-data validation.**  
   The reframing makes this no longer fatal, but the manuscript's claim that the framework is "suitable for application to real cohorts" is premature given weak simulated recovery and no empirical stress test.

7. **No shared-gene sensitivity analysis.**  
   Shared genes are still a plausible source of induced correlation between hallmark scores. This remains only a future-work item.

8. **No proper covariate or mixed-model treatment of tissue effects.**  
   The simulation includes tissue offsets, but the main network analyses appear to use pooled samples. Tissue-only baselines help with prediction, but they do not resolve whether network edges are confounded by tissue composition.

## 4. New Issues Found in the Revision

1. **Manuscript and result files are inconsistent.**  
   This is the most serious new issue. Examples:

   - Manuscript: correlation network density = 0.611 and 22 significant edges.  
     `analysis_summary.json`: density = 0.472 and 17 correlation edges.

   - Manuscript: Graphical Lasso identified 36 direct edges.  
     `analysis_summary.json`: partial-correlation edges = 12.

   - Manuscript: GB biological-age MAE = 4.9 years and Elastic Net R2 = 0.243.  
     `analysis_summary.json`: GB MAE = 15.9 years, GB R2 = 0.103, Elastic Net R2 = 0.188.

   - Manuscript: four age-dependent rewiring pairs are discussed as significant.  
     `results/age_rewiring_fdr.csv`: no comparison has `significant_fdr = True`.

   These mismatches suggest that the manuscript was only partially updated after rerunning the pipeline.

2. **Benchmark table contains placeholders instead of values.**  
   Table 4 reports "Reported" for every benchmark metric rather than the actual numbers. For a simulation-benchmarking paper, this is unacceptable. The benchmark values are the central evidence of the paper and must be in the main text.

3. **Benchmark results are weak but described as strong.**  
   The actual benchmark summary is:

   - Edge AUROC: 0.675
   - Best F1: 0.653
   - Centrality rho: 0.45
   - Top-3 hub overlap: 0.333
   - Module ARI: -0.052
   - Age-effect rho: 0.201
   - Stable edges: 22.2%

   These results support at most partial edge recovery under an idealized simulation. They do not support reliable recovery of modules, hubs, age effects, or robust edges.

4. **Age-rewiring claims are contradicted by the corrected FDR output.**  
   The old nominally significant rewiring results still appear in `interaction_age_changes.csv`, but the corrected `age_rewiring_fdr.csv` shows no FDR-significant tests. The manuscript should not claim detection of four significant rewiring events after FDR correction.

5. **Abstract makes unsupported quantitative claims.**  
   The abstract claims high-fidelity recovery, biological-age MAE around 5 years, significant outperformance of baselines, and robust edge stability. These claims are not supported by the current result summaries. Abstract claims should match the final generated outputs exactly.

6. **NARRATIVE_REPORT.md remains stale and misleading.**  
   The narrative report still presents the work largely as biological findings, reports GB MAE = 4.9 years, says the network density is 0.611, says there are four significant rewiring edges, and describes a PC algorithm identifying causal structure. This file conflicts with the revised framing and current result summaries.

7. **Permutation importance is still not fully cross-validated.**  
   The code fits final models on all data for feature importance and computes permutation importance on the same full data. This is acceptable if explicitly labelled as descriptive model interpretation, but not as cross-validated evidence of feature importance. The manuscript currently does not make this distinction clearly.

8. **Biological-age terminology remains somewhat problematic.**  
   Because the data are simulated and chronological age effects are embedded directly into latent hallmark factors, "biological age acceleration" is not a biological phenotype. It is a residual from a simulated age-prediction task. The manuscript should use more cautious language such as "age-prediction residual" unless applied to real data.

9. **PC algorithm description may still overstate orientation validity.**  
   Even with caveats, the text says the orientation pattern is consistent with hierarchy and discusses information flow. Given the approximate implementation and simulated cross-sectional data, this should be moved to a secondary exploratory analysis and not treated as a framework pillar.

10. **The benchmark design itself needs more detail.**  
    The manuscript does not clearly define what counts as a "true edge," how absent edges are defined from the embedded correlation matrix, what threshold is used, how module ground truth is specified, or how robustness is measured. Without these details, the reported AUROC/F1/ARI/stability values are hard to interpret.

## 5. Specific Remaining Suggestions, Prioritized

1. **Synchronize the manuscript with the current pipeline outputs.**  
   Regenerate all manuscript numbers from the current result files. Fix density, edge counts, partial-correlation edge count, biological-age metrics, centrality values, rewiring significance, and benchmark values. This is mandatory before any further review.

2. **Report the actual benchmark values in the main text and interpret them honestly.**  
   Replace Table 4 placeholders with numeric values. The current results should be described as mixed: modest edge recovery, weak centrality and age-effect recovery, failed module recovery, and low edge stability.

3. **Revise the abstract and conclusions to match the weak benchmark results.**  
   Remove "high fidelity," "reliably recovers," "faithfully capture," and "robustness confirmed" unless new benchmark results actually support those claims. A more accurate claim is that the pipeline partially recovers embedded edge structure but does not yet reliably recover modules, hubs, age effects, or stable topology.

4. **Correct the age-rewiring section.**  
   State that nominal young-old differences were observed in some pairs but none survived BH-FDR correction in the current output. Do not claim detection of significant rewiring after FDR correction.

5. **Correct the biological-age section.**  
   Use the out-of-fold MAE and R2 from `analysis_summary.json` or the generated source table. If MAE is 15.9 years, explain that age prediction is modest. Compare directly against mean-age and tissue-only baselines with confidence intervals.

6. **Fix the partial-correlation section.**  
   The text should report 12 edges, not 36, if `analysis_summary.json` is current. Also soften "direct relationships" and "removes indirect associations"; Graphical Lasso estimates conditional associations under model assumptions.

7. **Remove remaining causal language from figure captions and plot labels.**  
   Rename Figure 11 and related code labels from "causal" to "conditional-dependence orientation." Avoid "putative causal relationships" in captions.

8. **Update or remove NARRATIVE_REPORT.md.**  
   This file still reflects the rejected Round 1 framing and stale metrics. It should be revised to match the simulation-benchmarking paper or excluded from the manuscript package.

9. **Add a benchmark methods subsection with exact definitions.**  
   Define true edges, absent edges, thresholds, AUROC labels, F1 threshold selection, ground-truth modules, centrality ground truth, age-effect ground truth, and stability criteria. Include uncertainty across seeds, not only point estimates.

10. **Add sensitivity analyses that explain poor benchmark recovery.**  
    Vary sample size, noise level, gene-set size, shared-gene overlap, correlation threshold, Graphical Lasso penalty, and tissue imbalance. The current weak ARI and low stability require diagnosis.

11. **Add shared-gene and scoring-method sensitivity analyses.**  
    Recompute hallmark scores and networks after removing multi-hallmark genes, and compare simplified ssGSEA against at least one standard pathway scoring method.

12. **Reframe tissue-specific and biological-age results as simulation recovery only.**  
    Avoid language implying real brain, blood, muscle, skin, or liver biology. Say the pipeline recovered embedded tissue offsets and age trends, or failed to recover them where appropriate.

13. **Report confidence intervals wherever claimed in the Methods.**  
    If bootstrap CIs are implemented, include them in the main Results for ML metrics and key correlations, or move the claim to supplementary materials with clear file references.

14. **Use baselines in interpretation, not just in code.**  
    The Results should report majority-class, stratified-random, tissue-only, and mean-age baseline performance numerically. "Outperforming baselines" should be supported with values and CIs.

15. **Clarify intended contribution and venue fit.**  
    If this remains a simulation benchmark, the contribution should be a transparent evaluation of a pipeline, including failures. It is not yet a validated biological discovery tool. A convincing revision would foreground the limitations revealed by the benchmark rather than presenting the benchmark as uniformly successful.

## Summary

The revision deserves credit: the authors addressed the largest ethical and scientific problem from Round 1 by making the simulation nature explicit and by adding the right categories of validation. But the current manuscript is internally inconsistent and overstates the actual benchmark performance. The next revision should be driven by the generated result files, not by the narrative the authors hoped the benchmark would support. If the authors honestly report the current values, diagnose the weak recovery, and add sensitivity analyses, this could become a useful methods-demonstration paper. As written, it still requires major revision.
