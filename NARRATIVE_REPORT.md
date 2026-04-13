# Narrative Report: Simulation-Benchmarked Framework for Hallmark of Aging Interactions

**NOTE: All results are from SIMULATED data with known embedded structure. This is a framework validation study, not empirical biological discovery.**

## Research Question
Can a computational pipeline based on ssGSEA scoring, multi-layer network analysis, and machine learning reliably recover known inter-hallmark interaction structures from noisy gene expression data?

## Key Results (from current pipeline output)

### Network Analysis
- **17 Bonferroni-corrected significant edges** (density = 0.472)
- **12 Graphical Lasso partial correlation edges** (genuinely sparse)
- Two modules detected: damage-accumulation (GI, EA, DNS, MD) and tissue-decline (CS, SCE, AIC, TA, LP)
- Hub hallmarks: SCE (eigenvector = 0.516), AIC (eigenvector = 0.455)

### Machine Learning (all cross-validated, with Pipeline to prevent leakage)
- Age-group classification: RF 50.0%, GB 46.3%, LR 49.7% (baselines: majority 43.3%)
- Biological age: GB MAE = 15.9 years, R² = 0.103; EN MAE = 15.7, R² = 0.188
- Cross-hallmark prediction: AIC R²=0.434, SCE R²=0.436 (integrative most predictable); TA R²=0.103 (primary least predictable)

### Ground-Truth Benchmarking (HONEST ASSESSMENT)
- Edge AUROC: 0.675 (moderate)
- Best F1: 0.653 (moderate)
- Centrality rank rho: 0.45 (weak-to-moderate, p=0.22)
- Top-3 hub overlap: 33% (low)
- Module ARI: -0.052 (failed)
- Age effect recovery rho: 0.201 (weak)
- Edge stability: 22% (low)

### Age-Dependent Rewiring
- **No pairs survived FDR correction** (important negative result)

## Conclusion
The framework demonstrates moderate edge detection ability but weak module recovery, hub identification, and edge stability. The ssGSEA scoring bottleneck likely explains the gap. Application to real GTEx data with alternative scoring methods is the critical next step.

## Files
- Paper: paper/main.tex → paper/main.pdf (16 pages)
- Code: src/ (hallmark_genes.py, data_acquisition.py, network_analysis.py, ml_models.py, benchmarking.py, visualization.py)
- Results: results/ (17 figures, 15+ CSV/JSON files)
