# Emergence Suite — Phase A

Foundational experiments testing spontaneous Lorentzian-like phase separation in constraint networks.

## Scripts

| # | File | Tests | Key result |
|---|------|-------|------------|
| 01 | `01_biphasage.py` | Spontaneous phase separation | Δ = 0.31 ± 0.06 |
| 02 | `02_robustness.py` | 6 cohesion metrics, finite-size scaling | 5/6 metrics, plateau N=50–300 |
| 03 | `03_arrow.py` | Directional coherence of temporal arrow | Coherence ≈ 91% |
| 04 | `04_continuum.py` | Hausdorff dimension on RGG torus | d_H ≈ 2.84 (d=3), 3.38 (d=4) |
| 05 | `05_universality.py` | 5 graph families + permutation null | 5/5 pass, Cohen's d > 8 |
| 06 | `06_locality.py` | Local d_H homogeneity, isotropy | CV = 0.08, isotropy = 0.82 (d=3) |
| 07 | `07_coarsegraining.py` | Robust CG aggregators | 59% retention (d=3) |

## Running

```bash
pip install numpy scipy matplotlib networkx
python 01_biphasage.py
```

Each script is self-contained. All use fixed random seeds (documented in `seeds_used.txt`).

## Naming Convention

Scripts were originally numbered by development order (v4, v5, v7d, v8–v10).
Renamed here to logical order (01–07). Mapping:

| Current | Original | Content |
|---------|----------|---------|
| 01_biphasage | v5_biphasage | Phase separation |
| 02_robustness | v5b_robustness | Metric + scaling |
| 03_arrow | v4_arrow | Directional coherence |
| 04_continuum | v7d_continuum | Hausdorff dimension |
| 05_universality | v8_universality | 5 graph families |
| 06_locality | v9_locality | Local homogeneity |
| 07_coarsegraining | v10_rg | Coarse-graining |

## Detailed results

See `RESULTS.md` and `results.csv` for full tables.
