# Review Round 1

## 1. Overall Assessment

**Recommendation: Reject**  
**Score: 2/10**

The manuscript is not suitable for a top bioinformatics journal in its current form. The central problem is that the "multi-tissue transcriptomic data" are simulated with an explicitly embedded hallmark correlation structure and age/tissue effects, yet the paper presents the recovered network, hierarchy, tissue patterns, and age rewiring as biological discoveries. This creates circular validation: the main conclusions largely recover assumptions used to generate the data.

The work could become a useful methods demonstration or teaching pipeline if reframed transparently as a simulation study with rigorous benchmarking against known ground truth. As an empirical aging bioinformatics paper, however, it lacks independent biological data, external validation, adequate statistical correction, and defensible causal or intervention claims.

## 2. Summary

The manuscript proposes a pipeline for quantifying cross-talk among the nine hallmarks of aging using curated hallmark gene sets, ssGSEA-like scoring, correlation/partial-correlation/PC networks, and machine-learning models for age prediction and cross-hallmark predictability. It reports a two-module network architecture, hub roles for stem cell exhaustion and altered intercellular communication, age-dependent network rewiring, and hallmark-based biological age estimation.

The contribution is best interpreted as a simulated proof-of-concept workflow rather than a discovery study. The reported biological conclusions are not convincingly supported because the data-generation process already encodes the claimed hallmark correlations, age effects, and tissue effects.

## 3. Strengths

1. **Clear conceptual framing.** The manuscript is organized around a coherent biological question: whether hallmark-level scores can reveal systems-level aging architecture.

2. **Multi-method pipeline.** The authors combine gene set scoring, correlation networks, sparse partial correlations, causal-structure learning, cross-hallmark prediction, PCA, and age modeling. This is a reasonable computational scaffold for a future empirical study.

3. **Readable presentation.** The manuscript is easy to follow, with clear sectioning and interpretable tables/figures. The narrative report concisely captures the main claims.

4. **Attempt to connect results to biology.** The discussion links reported edges and modules to known mechanisms such as DNA damage, telomere attrition, SASP/inflammaging, nutrient sensing, and mitochondrial dysfunction.

5. **Reproducible code appears to exist locally.** The repository contains source modules, generated CSVs, figures, and a pipeline script. This is a useful start, although reproducibility is incomplete and the simulation nature must be explicit.

## 4. Weaknesses Ranked by Severity

### 1. Fatal circularity: the data are simulated with the expected biological structure

The Methods state that the authors "constructed the dataset with latent hallmark activity factors correlated according to published interaction evidence." The code confirms this: `generate_simulated_expression_data()` embeds a predefined correlation matrix, age effects, and tissue effects before generating expression values. Therefore, the reported two-module architecture, hub hallmarks, tissue-specific patterns, and age associations are not independent findings.

This invalidates the core biological claims. For example, the manuscript concludes that the analysis provides "quantitative evidence for the hierarchical organization of aging hallmarks," but the latent data generator already imposes a literature-derived hallmark correlation structure and age trends. Recovering these patterns is expected and does not constitute evidence.

### 2. Misleading presentation as empirical multi-tissue transcriptomics

The manuscript repeatedly describes "multi-tissue transcriptomic data comprising 300 samples across five tissue types" without clearly stating in the abstract and results that these are simulated samples. This is a major transparency problem. A reader could easily infer that the analysis used real GTEx, GEO, or other human aging transcriptomic data.

The paper should not claim "biological age estimation," "tissue-specific hallmark signatures," or "age-dependent rewiring" as empirical results when the samples, ages, tissues, hallmark activities, and expression values were generated synthetically.

### 3. Causal claims are unsupported and methodologically weak

The PC algorithm is presented as identifying "putative causal relationships" and a "predominantly feedforward causal flow." This is not justified. The data are cross-sectional, simulated, and include latent shared drivers. PC assumptions are not discussed: causal sufficiency, faithfulness, no hidden confounding, correct conditional independence tests, adequate sample size, and acyclicity.

The implementation is also only an approximation: it starts from marginal independence, conditions only on one variable, uses residualized Spearman correlations, and then orients edges heuristically. With only nine variables and simulated correlations, the result should not be interpreted as causal architecture. At most it is an exploratory directed graph under strong assumptions.

### 4. Statistical rigor is insufficient, especially multiple testing and uncertainty

The manuscript reports four age-dependent rewiring edges at `p < 0.05`, but 36 hallmark pairs are tested. The CSV shows the nominally significant p-values are 0.0303, 0.0476, 0.0161, and 0.00138. After Bonferroni correction over 36 tests, only the 0.00138 result would remain below 0.05; under FDR, the number may also decrease. The text does not report adjusted p-values, confidence intervals, or effect-size uncertainty.

Similarly, correlation-network significance is described as Bonferroni-corrected in the Methods, but the code appears to threshold raw p-values at 0.05. This is a direct inconsistency between manuscript and implementation. Many other analyses report point estimates only, with no confidence intervals, bootstrap stability, permutation nulls, or sensitivity analyses.

### 5. Machine-learning evaluation is overclaimed and partly in-sample

The biological age section reports Gradient Boosting MAE = 4.9 years, but the code computes MAE after fitting the model on all samples and predicting the same samples. The reported `gb_r2` in `analysis_summary.json` is cross-validated R2 = 0.100, which is weak, while the MAE is in-sample and optimistic. Elastic Net R2 = 0.243 is also apparently in-sample.

Age-group classification accuracy of 50.3% is modest and evaluated on simulated labels generated with known age effects. The binomial test against 33.3% chance is not the right primary benchmark for a 5-fold cross-validation setting with possible class imbalance and repeated model selection. A permutation test preserving tissue/age structure, or nested cross-validation with confidence intervals, would be more appropriate.

### 6. Gene set curation is under-specified and potentially biased

The authors state that genes were assigned to a "primary hallmark based on functional annotation," but do not provide a reproducible curation protocol, inclusion/exclusion rules, versioned database accessions, or independent curator agreement. Several genes are multifunctional, and forcing them into a primary hallmark can bias downstream ssGSEA scores and inferred cross-talk.

The gene overlap analysis is also not well integrated. The paper reports 16 cross-talk nodes but does not quantify whether overlap exceeds expectation, how duplicated/multifunctional genes affect correlation between hallmark scores, or whether results persist after removing shared genes.

### 7. The ssGSEA implementation is simplified and not validated

The scoring equation is a rank-average transformation, not standard ssGSEA as implemented in GSVA, ssGSEA2.0, singscore, or related tools. It does not discuss directionality of genes, mixed up/down gene sets, gene-set size effects, or calibration against known pathway activation. Since expression was simulated directly from latent hallmark factors, the scoring method is not meaningfully validated on real biological signal.

### 8. Partial-correlation network result is suspicious

The manuscript states that Graphical Lasso identified 36 direct edges among nine hallmarks, i.e., all possible pairwise edges. A complete "sparse" partial-correlation network contradicts the purpose of sparse inverse covariance estimation and suggests the thresholding/regularization is not informative. The claim that this network "removes indirect associations" is therefore not credible.

### 9. Biological interpretation exceeds the evidence

The paper makes intervention-oriented claims, e.g., upstream hallmark targeting may have cascading beneficial effects while targeting integrative hallmarks may provide "only symptomatic relief." These statements are speculative even for real observational data and are not supportable from simulated cross-sectional correlations.

Tissue-specific interpretations are also overextended. For example, brain aging being "dominated by stem cell exhaustion" is biologically contentious from bulk transcriptomic hallmark scores and impossible to infer from simulated tissue effects. Similar issues apply to liver telomere attrition and blood inflammaging claims.

### 10. Novelty is limited

The conceptual hierarchy of primary, antagonistic, and integrative hallmarks is inherited from López-Otín et al. The specific edges discussed are well-known mechanisms from the aging literature. The paper's novelty would need to come from empirical quantification across real cohorts, robust benchmarking, or a new validated method. In its current form, it mostly re-encodes known relationships in a simulated dataset and recovers them.

## 5. Missing Experiments or Analyses

1. **Independent real-data validation.** Apply the full pipeline to real human aging transcriptomic cohorts, ideally GTEx or well-curated GEO/ArrayExpress datasets with age, tissue, sex, batch, RIN/PMI where applicable, and disease status metadata.

2. **External replication.** Train/derive networks in one cohort and replicate edge signs, modules, centrality ranks, and age associations in another cohort. Report reproducibility metrics, not only visual agreement.

3. **Simulation benchmarking against known ground truth.** If simulation remains part of the paper, explicitly benchmark edge recovery, centrality recovery, module recovery, causal-orientation accuracy, and age-effect recovery against the known latent structure. Include null simulations and sensitivity to sample size, noise, tissue imbalance, and gene-set overlap.

4. **Multiple-testing correction.** Report BH-FDR and/or Bonferroni-adjusted p-values for all pairwise correlations, age-rewiring tests, tissue-specific age correlations, and any feature/edge selection claims.

5. **Uncertainty quantification.** Add bootstrapped confidence intervals for correlations, partial correlations, centrality metrics, module assignments, ML performance, and cross-hallmark R2 estimates.

6. **Confounder adjustment.** Model age-hallmark and hallmark-hallmark associations while adjusting for tissue, sex, batch, and other available covariates. For multi-tissue data, use mixed models or stratified analyses with appropriate interaction terms.

7. **Shared-gene sensitivity analysis.** Recompute hallmark scores and networks after removing all genes assigned to multiple hallmarks. This is essential because shared genes can induce artificial cross-hallmark correlation.

8. **Compare scoring methods.** Compare the simplified rank score against GSVA/ssGSEA, singscore, PLAGE, mean z-score, and pathway-level PCA. Assess whether network topology is robust to scoring method.

9. **Network null models.** Compare observed density, modularity, hub rankings, and edge weights against degree-preserving, gene-label permutation, sample-label permutation, and random gene-set nulls.

10. **Proper biological age evaluation.** Use nested cross-validation or an external test set. Report cross-validated MAE, R2, Pearson/Spearman correlation, calibration, and age-bias correction. Compare against simple baselines such as tissue-only, sex+tissue, and mean expression PCs.

11. **Causal validation or reframing.** Either remove causal claims or validate against perturbation/longitudinal data. At minimum, present PC results as exploratory conditional-dependence structure under restrictive assumptions.

12. **Functional enrichment of modules.** If modules are found in real data, test whether module genes are enriched for expected pathways, cell types, regulatory programs, or disease phenotypes.

## 6. Writing Quality

The manuscript is generally readable, but several passages overstate the evidence or obscure the simulated nature of the work.

1. **Abstract Methods:**  
   "Using curated gene sets ... on multi-tissue transcriptomic data (300 samples, 5 tissues, ages 20--90)."  
   This should explicitly say "simulated multi-tissue transcriptomic data" if that remains true. Otherwise it is misleading.

2. **Methods, Transcriptomic Data and Preprocessing:**  
   "To embed realistic inter-hallmark correlation structure... we constructed the dataset with latent hallmark activity factors..."  
   This sentence is crucial but buried. It should be in the abstract, title/summary, and first paragraph of Results. The manuscript must distinguish simulated inputs from discoveries.

3. **Abstract Results:**  
   "Network analysis revealed a highly interconnected hallmark landscape..."  
   This should be rewritten as "In simulated data with embedded hallmark correlations, the pipeline recovered..." unless real data are used.

4. **Conclusion:**  
   "provides quantitative evidence for the hierarchical organization of aging hallmarks"  
   This is unsupported. A more accurate statement would be that the pipeline can recover a hypothesized hierarchy from data generated under that hypothesis.

5. **Discussion:**  
   "interventions targeting upstream hallmarks may have cascading beneficial effects..."  
   This is too speculative and should be removed or clearly labeled as a hypothesis from prior literature, not a conclusion from this analysis.

6. **Partial correlation caption:**  
   "This network removes indirect associations mediated through third-party hallmarks."  
   Too strong. Graphical Lasso estimates conditional associations under modeling assumptions; it does not guarantee removal of indirect biological mediation, especially when all 36 edges are retained.

7. **Causal Structure section:**  
   "Genomic Instability showed causal influence..."  
   Replace "causal influence" with "an oriented edge was inferred under the PC heuristic" or remove this result from the main claims.

8. **Tissue-specific section:**  
   "suggesting pronounced stem cell exhaustion with brain aging" and similar tissue interpretations.  
   These passages are unsupported without real data and cell-composition controls.

There is also a LaTeX/reference issue: `Table~\ref{tab:cross_talk_evidence}` is referenced in the Methods, but no such table appears in the manuscript. The build log reports this as an undefined reference.

## 7. Technical Issues

1. **Manuscript-code mismatch on Bonferroni correction.** The Methods say correlation edges use `p < 0.05` after Bonferroni correction, but the implementation thresholds raw p-values. This must be corrected.

2. **Age rewiring uses nominal p-values.** The four reported rewiring edges are not corrected for 36 pairwise tests. Most do not survive Bonferroni correction.

3. **In-sample MAE for biological age.** The reported GB MAE is computed on predictions from a model fit to all data, not out-of-fold predictions. This inflates performance.

4. **Scaler fitted before cross-validation.** Standardization is performed before cross-validation in the ML code. With only nine features this may be minor, but preprocessing should be inside a pipeline to avoid leakage.

5. **Elastic Net performance is in-sample.** The reported Elastic Net R2 and coefficients are from a model fit on the full dataset. This should not be presented as predictive performance.

6. **PC implementation is incomplete.** Conditioning only on single variables and heuristic orientation is not a valid full PC algorithm. It does not provide reliable causal direction.

7. **Graphical Lasso returns a complete network.** A partial-correlation network with 36/36 edges should not be described as sparse or as revealing selected direct relationships.

8. **No model selection protocol.** Hyperparameters are chosen without nested CV or held-out validation. Feature importances from tree models are not accompanied by stability or permutation importance.

9. **No correction for tissue composition or covariates.** Tissue is a major source of transcriptomic variation. Pooling tissues without mixed models or tissue-adjusted residuals can induce spurious correlations.

10. **Gene directionality is ignored.** Some hallmark genes increase with age, others decrease, and some are context-specific. A single unsigned enrichment score may not represent pathway activity.

11. **Simulated expression is generated from hallmark scores, then hallmark scores are recomputed from expression.** This is a circular validation of the scoring pipeline unless explicitly framed as recovery from known latent variables.

12. **No availability details.** "The complete pipeline is available as open-source code" is not sufficient for publication. Provide repository URL, commit hash, environment lock file, commands, random seeds, input data provenance, and expected outputs.

## 8. Concrete Suggestions

1. **Decide what the paper is.** Either reframe it as a simulation/benchmarking paper or replace the simulated dataset with real aging transcriptomic data. It cannot be presented as empirical biological discovery in the current form.

2. **If reframed as simulation, change the title, abstract, results, and conclusions.** Use language such as "simulation framework," "known ground truth," and "pipeline recovery." Remove claims about discovered aging biology.

3. **If aiming for a bioinformatics discovery paper, analyze real cohorts.** Use GTEx and/or multiple public aging datasets, document preprocessing rigorously, and reserve at least one cohort for external validation.

4. **Make all data provenance explicit.** State which samples are simulated, which are real, how ages/tissues were generated or obtained, and what assumptions were embedded.

5. **Add a ground-truth benchmark section.** Since the simulation has a known latent correlation matrix, report precision, recall, F1, AUROC/AUPRC for edge recovery; module adjusted Rand index; rank correlation for centrality; and robustness across random seeds.

6. **Correct all statistical testing.** Apply BH-FDR or Bonferroni correction consistently. Report adjusted p-values and confidence intervals for all edge and rewiring claims.

7. **Fix ML evaluation.** Use scikit-learn pipelines, nested CV where appropriate, and out-of-fold predictions for MAE/R2. Compare against simple baselines and report uncertainty.

8. **Remove or heavily qualify causal language.** Rename "causal network" to "conditional-dependence orientation analysis" unless supported by longitudinal or perturbational data.

9. **Add sensitivity analyses.** Vary correlation thresholds, Graphical Lasso penalties, age bins, scoring methods, gene-set sizes, shared-gene removal, and random seeds.

10. **Control for tissue and sex.** Perform tissue-stratified or mixed-effect modeling and show that hallmark associations are not artifacts of tissue composition.

11. **Validate gene sets.** Provide the full gene set table as supplementary material with sources, assignment rationale, multifunctional gene labels, and database versions.

12. **Temper biological interpretations.** Restrict claims to what is actually supported. Do not infer intervention strategy, tissue vulnerability, or causal hierarchy from simulated cross-sectional associations.

13. **Fix manuscript inconsistencies.** Add or remove the missing `tab:cross_talk_evidence` reference, align the Results with `analysis_summary.json`, and explicitly state that the strongest correlation reported in the JSON is negative CS-SCE (`r = -0.572`), which requires biological explanation given that both are described as increasing with age.

14. **Improve reproducibility.** Include exact run commands, environment versions, random seeds, a dependency lock file, and checksums or scripts for regenerating every table and figure.

15. **Reconsider publication venue expectations.** For Bioinformatics/PLOS Computational Biology/Genome Biology, the current manuscript would need either a validated new method with rigorous benchmarks or convincing biological insights from real data. The present version provides neither at sufficient depth.

