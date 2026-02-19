# TCGE Phase B — Dimensional Analysis (Exploratory)

**Status:** Exploratory. Not part of v1.0 claims.

## Scripts

| Script | What it does |
|--------|-------------|
| `dim_scan.py` | Scan d=2→7 at fixed N, plus scaled-N continuum diagnostics for d=2,3,4 |
| `dim_robustness.py` | Vary ⟨k⟩ at fixed d to test whether the dimensional sweet spot is real |

## Key finding

The triangle protector does not control the **amount** of polarization —
it controls its **spatial distribution**.

The linear law Δ ≈ −0.034 × ⟨tri⟩ + 0.51 (r = −0.98, 15 conditions) holds
when varying ⟨k⟩, but the causal control (degree-preserving rewiring) reveals
the deeper mechanism: without triangles, ALL edges polarize maximally
(⟨|α|⟩ → 0.98). Triangles selectively protect cohesive edges, creating the
heterogeneity between proto-temporal and proto-spatial populations that
constitutes the Lorentzian-like phase separation.

Embedding dimension has no independent effect (coefficient ≈ 0 in multivariate fit).

The apparent "sweet spot" at d ≈ 4–5 seen in the d-scan at ⟨k⟩ = 8 is not a dimensional
selection. It reflects the triangle density at that degree: dimensions 4–5 happen to have
⟨tri⟩ ≈ 2.5–3 at ⟨k⟩ = 8, which is the range where biphasage is strongest.

## Interpretation

Triangles are the coherence protector in the cost functional. More triangles → more edges
resist polarization → lower biphasage. The relationship is mechanistically direct: it is
the competition term C_protect × tri(e) in the gradient that controls the equilibrium.

## What this means for TCGE

- TCGE does **not** select d = 4
- Dimension is a proxy for triangle density at a given mean degree
- The triangle protector's role is not to limit polarization, but to **structure** it:
  it creates the spatial heterogeneity between proto-temporal and proto-spatial edges
- Without the protector, all edges polarize uniformly — no Lorentzian signature
- The protector is therefore **necessary** for emergent spacetime structure, not just
  a regularizer

## Causal control: rewiring (dim_causal.py)

Degree-preserving rewiring (Maslov-Sneppen) at fixed d=3, ⟨k⟩=8 reveals
the deepest result: destroying triangles does NOT reduce polarization.
It makes ALL edges polarize to |α| ≈ 0.98.

| rewire% | ⟨tri⟩ | ⟨|α|⟩ | % polarized |
|---------|-------|-------|-------------|
| 0%      | 3.80  | 0.43  | 36%         |
| 10%     | 2.08  | 0.65  | 61%         |
| 50%     | 0.21  | 0.95  | 95%         |
| 90%     | 0.04  | 0.98  | 98%         |

**The triangle protector does not control the amount of polarization —
it controls its spatial distribution.** Without triangles, polarization
is uniform and the T/S distinction vanishes. With triangles, cohesive
edges resist polarization (proto-spatial) while others polarize
(proto-temporal), creating the heterogeneity required for Lorentzian
signature.

## Trajectory of understanding

1. **Initial hypothesis** (dim_scan.py): Δ ~ tri/diam, r = 0.986 on 6 points
2. **Falsified** (dim_robustness.py): varying ⟨k⟩ reverses the sign
3. **Refined hypothesis**: Δ ~ −⟨tri⟩, r = −0.98 on 15 points
4. **Deepened** (dim_causal.py): tri doesn't control polarization magnitude
   but its spatial heterogeneity — the protector creates the T/S structure

Each step killed the previous interpretation and replaced it with
something more fundamental. This is documented honestly as a model
of scientific self-correction.
