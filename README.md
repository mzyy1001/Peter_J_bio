# Hallmarks of Aging Interaction Analysis Framework

A simulation-benchmarked computational framework for quantifying cross-talk among the nine hallmarks of aging.

## Quick Start

```bash
# Setup
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Run full pipeline
python run_pipeline.py
```

## Output

- `results/` — CSV/JSON data files, benchmark metrics
- `results/figures/` — 20 publication-quality PDF figures
- `paper/main.tex` — LaTeX manuscript (compile with `pdflatex` + `bibtex`)

## Project Structure

```
src/
  hallmark_genes.py      # Curated gene sets for 9 hallmarks (250 genes)
  data_acquisition.py    # Simulated multi-tissue expression data
  network_analysis.py    # Correlation, partial correlation, PC networks
  ml_models.py           # Age classification, biological age, cross-prediction
  benchmarking.py        # Ground-truth recovery benchmarks
  sensitivity.py         # Shared-gene, scoring method, tissue adjustment
  visualization.py       # Publication-quality figures
run_pipeline.py          # Main pipeline (runs all analyses)
paper/
  main.tex               # LaTeX manuscript
  references.bib         # BibTeX references
```

## Key Parameters

- `n_samples=300` — number of simulated samples
- `seed=42` — random seed for reproducibility
- 5 tissues: blood, brain, muscle, liver, skin
- Age range: 20-90 years

## Citation

If you use this framework, please cite the accompanying manuscript.

## License

MIT
