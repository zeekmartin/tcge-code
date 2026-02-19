# TCGE — Computational Evidence for Lorentzian Phase Separation

**Release v1.0** — February 2026

## Overview

This repository provides computational verification that the TCGE (Theory of 
Emergent Global Constraints) framework produces spontaneous Lorentzian phase 
separation on constraint networks. Specifically:

1. **Biphasage** (v5): Edges in high-clustering regions remain symmetric 
   (proto-spatial), while bottleneck edges polarize (proto-temporal), without 
   injecting any T/S label.

2. **Robustness** (v5b): The phase separation is not specific to triangles — 
   5 out of 6 local cohesion metrics produce equivalent results. It is stable 
   across graph sizes (N=30–300).

3. **Arrow** (v4): Directional asymmetry emerges from a product cost on 
   irregular graphs with intrinsic node properties (degree).

4. **Continuum compatibility** (v7d): On Random Geometric Graphs with periodic 
   boundary conditions (torus), Hausdorff dimension exhibits stable plateaus 
   (d_H ≈ 2.84 for d=3, d_H ≈ 3.38 for d=4). The biphasage mechanism 
   preserves the underlying geometric dimensionality.

## Scripts

| File | Purpose | Key result |
|------|---------|------------|
| `v4_arrow.py` | Directional arrow on irregular graphs | Arrow coherence ~91% |
| `v5_biphasage.py` | Phase separation via cohesion protection | Δ = 0.31 ± 0.06 |
| `v5b_robustness.py` | Robustness, TCGE-native reformulation, scaling | 5/6 metrics pass |
| `v7d_continuum.py` | Continuum diagnostics on RGG torus | d_H plateau stable |

## Reproducing results

Each script is self-contained. Requirements: Python 3.8+, NumPy, Matplotlib, SciPy.
For v7d: also `python-louvain` (`pip install python-louvain`).

```bash
python v5_biphasage.py      # ~15s, produces fig_biphasage.png
python v5b_robustness.py    # ~15s, produces fig_robustness.png
python v4_arrow.py          # ~10s, produces fig_arrow.png
python v7d_continuum.py     # ~50s, produces fig_continuum.png
```

All random seeds are fixed for full reproducibility. See `seeds_used.txt`.

## Figures

- `fig_biphasage.png` — Phase separation Δ vs parameters, convergence, distribution
- `fig_robustness.png` — Robustness across cohesion metrics, scaling, TCGE-native comparison
- `fig_arrow.png` — Arrow mechanism on irregular graph, coherence measurements
- `fig_continuum.png` — d_H(r) plateau, d_s(t), volume-radius, CG retention on torus

## Key formulation

> Phase separation between polarized (proto-temporal) and symmetric (proto-spatial)
> edges emerges from competition between a symmetry-breaking product cost and local
> coherence protection. The protection term measures compatibility destruction cost:
> the number of coherent triplets broken by polarization. The separation is robust
> across cohesion metrics, stable across graph sizes, and compatible with geometric
> substrates where it preserves Hausdorff and spectral dimensions.

## Limitations

- d_H remains ~0.5 below embedding dimension (finite-size effect at N=5000)
- d_s exhibits significant variance (±0.5) due to limited mixing window
- Coarse-graining retention reaches 44–60% (d=3); 50% threshold not yet 
  consistently exceeded
- These are toy models on constraint networks, not full spacetime simulations

## Author

Zeek (TCGE project)  
AI-assisted implementation: Claude (Anthropic)  
All conceptual direction and scientific judgment: human researcher  
