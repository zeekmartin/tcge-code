# TCGE Emergence Suite v1.0

Computational evidence for Lorentzian phase separation in a Boolean constraint framework.

**Repository:** Part of [tcge-code](https://github.com/zeekmartin/tcge-code), subdirectory `gap-emergence/`

## Overview

These scripts test whether a single cost functional on constraint graphs
spontaneously produces Lorentzian-like phase separation — edges splitting into
proto-temporal (high |α|) and proto-spatial (low |α|) classes — without injecting
any causal or geometric labels.

The mechanism uses three competing terms:
- A **symmetry-breaking** product term that favors polarization
- A **local coherence protector** (triangle count) that resists polarization
  on locally cohesive edges
- A **directional** term coupling asymmetry to intrinsic node properties

## Scripts

| # | Script | Tests | Key result |
|---|--------|-------|------------|
| 01 | `01_biphasage.py` | Spontaneous phase separation | Δ = 0.31 ± 0.06 |
| 02 | `02_robustness.py` | 6 cohesion metrics, finite-size scaling | 5/6 metrics, plateau N=50–300 |
| 03 | `03_arrow.py` | Directional coherence of temporal arrow | Coherence ≈ 91% |
| 04 | `04_continuum.py` | Hausdorff dimension on RGG torus | d_H ≈ 2.84 (d=3), 3.38 (d=4) |
| 05 | `05_universality.py` | 5 graph families + permutation null | 5/5 pass, p < 0.005, Cohen's d > 8 |
| 06 | `06_locality.py` | Local d_H homogeneity, isotropy | CV = 0.08, isotropy = 0.82 (d=3) |
| 07 | `07_coarsegraining.py` | Robust CG aggregators | 59% retention (d=3), structural limit |

## Running

```bash
pip install numpy scipy matplotlib python-louvain networkx
python 01_biphasage.py
```

Each script is self-contained. All use fixed random seeds (documented in `seeds_used.txt`)
for full reproducibility.

## Key findings

1. **Phase separation is spontaneous.** The cost functional produces two distinct edge
   populations without any injected labels (Δ = 0.31 ± 0.06, p < 0.005 vs permutation null).

2. **Phase separation is universal.** The mechanism works on 5 graph families (ER, RGG,
   Watts-Strogatz, Barabási-Albert, Configuration Model) with large effect sizes
   (Cohen's d > 8 in all cases). Amplitude is reduced on scale-free graphs, consistent
   with concentrated triangle heterogeneity (r = −0.51).

3. **The mechanism is compatible with continuum structure.** On geometric substrates
   (RGG torus), Hausdorff dimension plateaus are stable (σ < 0.03) and locally homogeneous
   (CV = 0.08 in d=3). Biphasage does not distort the underlying geometry.

4. **Coarse-graining partially preserves the signal.** Approximately 40–60% retention in
   d=3, with attenuation consistent with statistical averaging of fine-scale edge
   properties — a structural limit, not a failure of the mechanism.

## Limitations

- All models are toy models (N ≤ 10,000). The fundamental hypothesis — that physical
  reality functions as a compatibility network — is not tested here.
- Spectral dimension d_s remains noisy (σ ≈ 0.5).
- Coarse-graining does not reach stable >50% retention in d=4.
- The question of *why d=4* is not addressed.

## Citation

Venti, D. M. (2026). Theory of Emergent Global Constraints: Reality as Coherence.
Zenodo. https://doi.org/10.5281/zenodo.18297823

## Methodology note

Mathematical formalisations and computational implementations were developed with
AI system assistance. The author maintains full conceptual direction and responsibility
for all scientific claims and interpretations.

## License

MIT
