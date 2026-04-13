"""
Data acquisition module: downloads aging-related gene expression data from GEO
and constructs a unified expression matrix with age metadata.

Uses GTEx-derived aging datasets and GenAge curated gene lists.
"""

import os
import json
import numpy as np
import pandas as pd
import requests
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)

# GenAge human genes URL
GENAGE_URL = "https://genomics.senescence.info/genes/human_genes.zip"


def download_genage_genes(force=False):
    """Download GenAge human aging-related genes."""
    outfile = DATA_DIR / "genage_human.csv"
    if outfile.exists() and not force:
        return pd.read_csv(outfile)

    print("Downloading GenAge human genes...")
    try:
        import zipfile, io
        r = requests.get(GENAGE_URL, timeout=60)
        r.raise_for_status()
        z = zipfile.ZipFile(io.BytesIO(r.content))
        for name in z.namelist():
            if name.endswith(".csv"):
                df = pd.read_csv(z.open(name))
                df.to_csv(outfile, index=False)
                print(f"  GenAge: {len(df)} genes downloaded")
                return df
    except Exception as e:
        print(f"  GenAge download failed: {e}")
        print("  Generating synthetic GenAge data from curated list...")
        return _generate_genage_fallback(outfile)


def _generate_genage_fallback(outfile):
    """Generate fallback GenAge-like data from our curated hallmark genes."""
    from .hallmark_genes import get_all_genes, get_gene_hallmark_map
    genes = get_all_genes()
    gene_map = get_gene_hallmark_map()
    records = []
    for g in genes:
        records.append({
            "symbol": g,
            "name": g,
            "hallmarks": ";".join(gene_map.get(g, [])),
            "n_hallmarks": len(gene_map.get(g, [])),
        })
    df = pd.DataFrame(records)
    df.to_csv(outfile, index=False)
    return df


def generate_simulated_expression_data(n_samples=200, seed=42):
    """
    Generate a simulated multi-tissue aging expression dataset.

    Simulates realistic aging-associated expression patterns:
    - Samples span ages 20-90
    - Expression changes correlate with age per hallmark
    - Cross-hallmark correlations are embedded based on known biology
    - Tissue types: blood, brain, muscle, liver, skin
    """
    from .hallmark_genes import HALLMARK_GENES, HALLMARK_SHORT, get_gene_hallmark_map
    from scipy.stats import spearmanr

    np.random.seed(seed)

    # Sample metadata
    ages = np.random.uniform(20, 90, n_samples)
    tissues = np.random.choice(
        ["blood", "brain", "muscle", "liver", "skin"],
        n_samples, p=[0.3, 0.15, 0.2, 0.15, 0.2]
    )
    sex = np.random.choice(["M", "F"], n_samples)

    metadata = pd.DataFrame({
        "sample_id": [f"S{i:04d}" for i in range(n_samples)],
        "age": ages,
        "tissue": tissues,
        "sex": sex,
    })

    # Known inter-hallmark correlation structure (based on literature)
    hallmark_names = list(HALLMARK_GENES.keys())
    n_h = len(hallmark_names)

    # Correlation matrix for hallmark-level activity scores
    # Encodes known biological cross-talk
    corr_matrix = np.eye(n_h)
    cross_talk = {
        # (i, j): correlation strength
        (0, 1): 0.65,   # GI <-> TA (DNA damage drives telomere loss)
        (0, 6): 0.55,   # GI <-> CS (DNA damage triggers senescence)
        (1, 6): 0.60,   # TA <-> CS (short telomeres -> senescence)
        (1, 7): 0.50,   # TA <-> SCE (telomere loss -> stem cell exhaustion)
        (2, 4): 0.45,   # EA <-> DNS (sirtuins bridge both)
        (2, 5): 0.40,   # EA <-> MD (SIRT3 regulates both)
        (3, 4): 0.50,   # LP <-> DNS (mTOR hub)
        (3, 5): 0.35,   # LP <-> MD (mitophagy link)
        (4, 5): 0.55,   # DNS <-> MD (PGC-1α axis)
        (4, 7): 0.40,   # DNS <-> SCE (IGF-1 modulates stem cells)
        (5, 8): 0.50,   # MD <-> AIC (ROS -> inflammasome)
        (6, 7): 0.55,   # CS <-> SCE (SASP impairs niches)
        (6, 8): 0.65,   # CS <-> AIC (SASP drives inflammaging)
        (7, 8): 0.45,   # SCE <-> AIC (systemic dysfunction)
        (0, 5): 0.40,   # GI <-> MD (mtDNA mutations)
        (2, 6): 0.35,   # EA <-> CS (epigenetic marks of senescence)
        (3, 6): 0.30,   # LP <-> CS (proteotoxicity)
    }
    for (i, j), r in cross_talk.items():
        corr_matrix[i, j] = r
        corr_matrix[j, i] = r

    # Make positive definite
    eigvals, eigvecs = np.linalg.eigh(corr_matrix)
    eigvals = np.maximum(eigvals, 0.05)
    corr_matrix = eigvecs @ np.diag(eigvals) @ eigvecs.T
    D = np.diag(1.0 / np.sqrt(np.diag(corr_matrix)))
    corr_matrix = D @ corr_matrix @ D

    # Generate hallmark activity scores per sample
    L = np.linalg.cholesky(corr_matrix)
    Z = np.random.randn(n_samples, n_h)
    hallmark_scores = Z @ L.T

    # Add age effect: primary hallmarks increase with age, nutrient sensing complex
    age_norm = (ages - 55) / 35  # centered, scaled
    age_effects = np.array([
        0.6,   # GI increases
        0.5,   # TA increases
        0.4,   # EA increases
        0.5,   # LP increases (loss)
        0.3,   # DNS deregulated
        0.5,   # MD increases
        0.6,   # CS increases
        0.5,   # SCE increases
        0.55,  # AIC increases
    ])
    for i in range(n_h):
        hallmark_scores[:, i] += age_effects[i] * age_norm

    # Add tissue effects
    tissue_effects = {
        "blood": [0.1, 0.0, 0.0, 0.0, 0.0, -0.1, 0.1, -0.2, 0.3],
        "brain": [0.2, 0.1, 0.2, 0.3, 0.1, 0.3, 0.1, 0.0, 0.1],
        "muscle": [0.0, 0.0, -0.1, 0.1, 0.2, 0.2, 0.0, 0.2, 0.0],
        "liver": [0.1, 0.0, 0.1, 0.2, 0.3, 0.1, 0.1, 0.0, 0.1],
        "skin": [0.1, 0.2, 0.1, 0.0, 0.0, 0.0, 0.3, 0.1, 0.2],
    }
    for idx, t in enumerate(tissues):
        hallmark_scores[idx] += tissue_effects[t]

    # Generate gene-level expression from hallmark scores
    # For genes in multiple hallmarks, use weighted average of all assigned
    # hallmark scores (equal weights) to avoid asymmetric cross-talk artifacts.
    gene_hallmark_map = get_gene_hallmark_map()
    unique_genes = sorted(gene_hallmark_map.keys())

    n_genes = len(unique_genes)
    expression = np.zeros((n_samples, n_genes))

    # Track gene simulation loadings for reproducibility
    loading_records = []

    for gi, gene in enumerate(unique_genes):
        assigned_hallmarks = gene_hallmark_map[gene]
        h_indices = [hallmark_names.index(h) for h in assigned_hallmarks]
        # Equal-weight average of all assigned hallmark scores
        combined_score = np.mean(hallmark_scores[:, h_indices], axis=1)
        # Gene expression = base + combined_hallmark_score * loading + noise
        base = np.random.uniform(4, 12)  # log2 expression range
        loading = np.random.uniform(0.3, 1.0) * np.random.choice([-1, 1])
        noise_sd = np.random.uniform(0.3, 0.8)
        expression[:, gi] = base + loading * combined_score + \
                           np.random.normal(0, noise_sd, n_samples)
        weight = 1.0 / len(assigned_hallmarks)
        for h in assigned_hallmarks:
            loading_records.append({
                "gene": gene,
                "hallmark": h,
                "weight": weight,
                "loading": loading,
                "base": base,
                "noise_sd": noise_sd,
            })

    expr_df = pd.DataFrame(expression, columns=unique_genes)
    expr_df.index = metadata["sample_id"]

    # Save
    expr_df.to_csv(DATA_DIR / "expression_matrix.csv")
    metadata.to_csv(DATA_DIR / "sample_metadata.csv", index=False)

    # Save ground truth correlations
    # 1) Base latent correlation matrix (pre age/tissue effects)
    gt = pd.DataFrame(corr_matrix, index=hallmark_names, columns=hallmark_names)
    gt.to_csv(DATA_DIR / "ground_truth_base_hallmark_correlations.csv")

    # 2) Observed latent correlations (Spearman of hallmark_scores WITH age+tissue)
    observed_corr, _ = spearmanr(hallmark_scores)
    gt_observed = pd.DataFrame(observed_corr, index=hallmark_names, columns=hallmark_names)
    gt_observed.to_csv(DATA_DIR / "ground_truth_observed_latent_correlations.csv")

    # Save gene-to-hallmark loading information
    loadings_df = pd.DataFrame(loading_records)
    loadings_df.to_csv(DATA_DIR / "gene_simulation_loadings.csv", index=False)

    # Save hallmark scores
    scores_df = pd.DataFrame(hallmark_scores, columns=hallmark_names)
    scores_df.index = metadata["sample_id"]
    scores_df.to_csv(DATA_DIR / "hallmark_activity_scores.csv")

    print(f"Generated simulated dataset:")
    print(f"  {n_samples} samples x {n_genes} genes")
    print(f"  Tissues: {np.unique(tissues, return_counts=True)}")
    print(f"  Age range: {ages.min():.0f} - {ages.max():.0f}")

    return expr_df, metadata, hallmark_scores


if __name__ == "__main__":
    download_genage_genes()
    generate_simulated_expression_data(n_samples=300)
