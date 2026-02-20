# TCGE — Computational Evidence

Computational evidence for emergent phase separation in constraint networks, developed as part of the Theory of Emergent Global Constraints (TCGE) programme.

**Preprint:** Venti, D. M. (2026). *Theory of Emergent Global Constraints: Reality as Coherence.* Zenodo. [doi:10.5281/zenodo.18297823](https://doi.org/10.5281/zenodo.18297823)

## Structure

The repository is organized chronologically by research phase:

```
tcge-code/
├── emergence_suite/          Phase A — Foundational tests (7 scripts)
├── dimensional_analysis/     Phase B — Exploratory (triangle ≠ dimension)
├── causal_controls/          Phase C — Core causal result ★
├── combinatorial/            Sorkin response (combinatorial reformulation)
└── archive/                  Development versions (reference only)
```

### Phase A: Emergence Suite

Seven self-contained experiments establishing that a cost functional on constraint graphs spontaneously produces Lorentzian-like phase separation.

| # | Test | Key result |
|---|------|------------|
| 01 | Biphasage | Δ = 0.31 ± 0.06, spontaneous separation |
| 02 | Robustness | 5/6 cohesion metrics, stable N=50–300 |
| 03 | Arrow | 91% directional coherence |
| 04 | Continuum | d_H ≈ 2.84 (3D), 3.38 (4D) on RGG torus |
| 05 | Universality | 5/5 graph families, Cohen's d > 8 |
| 06 | Locality | CV = 0.08, isotropy = 0.82 |
| 07 | Coarse-graining | 59% retention (3D) |

### Phase B: Dimensional Analysis

Exploratory work revealing that biphasage is controlled by local triangle abundance (⟨tri⟩), not embedding dimension. Key finding: Δ ~ −⟨tri⟩ (r = −0.98, 15 conditions). This motivated Phase C.

### Phase C: Causal Controls ★

**The core result.** Double dissociation establishing that signature-like phase separation requires both local cohesive structure (triangles) and a coherence-protection term. Includes bifurcation theorem proving structural inevitability and dynamic stability analysis.

| Test | Result |
|------|--------|
| Double dissociation | F(1,96) = 4497, p < 10⁻⁶⁰ |
| Bifurcation theorem | Universal for any convex g with g″(0) > 0 |
| Dynamic stability | Unique basin, corr > 0.94 at ε = 0.5 |
| Scale | Effect strengthens N=50→180 |

→ Start here: [`causal_controls/README.md`](causal_controls/README.md)

### Combinatorial Reformulation

Response to the objection that weighted DAGs are less simple than posets. Research log documenting the progression from bare degree (v1, fails) to directed boolean constraints (v5b, succeeds). Key conclusion: Lorentzian-like separation emerges from anisotropic constraint density alone — no real numbers required in the foundations.

## Quick Start

```bash
pip install numpy scipy matplotlib networkx
cd causal_controls/scripts
python 01_double_dissociation.py    # 5s — the central result
python 05_bifurcation_theorem.py    # 2s — the theorem
python 06_dynamic_stability.py      # 2s — basin uniqueness
```

All scripts are self-contained. Total runtime for Phase C: ~1 minute.

## Requirements

```
numpy >= 1.20
scipy >= 1.7
matplotlib >= 3.4
networkx >= 2.6
```

## Methodology

Mathematical formalization and computational implementation were developed in collaboration with AI systems (Claude, Anthropic). The author maintains full conceptual direction and responsibility for all scientific claims and interpretations.

## License

MIT
