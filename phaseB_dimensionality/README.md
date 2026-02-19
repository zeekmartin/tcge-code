# TCGE Phase B — Dimensional Analysis (Exploratory)

**Status:** Exploratory. Not part of v1.0 claims.

## Scripts

| Script | What it does |
|--------|-------------|
| `dim_scan.py` | Scan d=2→7 at fixed N, plus scaled-N continuum diagnostics for d=2,3,4 |
| `dim_robustness.py` | Vary ⟨k⟩ at fixed d to test whether the dimensional sweet spot is real |

## Key finding

The biphasage amplitude is controlled by a single quantity: **mean triangle count per edge**.

```
Δ ≈ −0.034 × ⟨tri⟩ + 0.51
r = −0.98, R² = 0.97, n = 15 conditions (d = 3–5, ⟨k⟩ = 5–20)
```

Embedding dimension has no independent effect (coefficient ≈ 0).

The apparent "sweet spot" at d ≈ 4–5 seen in the d-scan at ⟨k⟩ = 8 is not a dimensional
selection. It reflects the triangle density at that degree: dimensions 4–5 happen to have
⟨tri⟩ ≈ 2.5–3 at ⟨k⟩ = 8, which is the range where biphasage is strongest.

## Interpretation

Triangles are the coherence protector in the cost functional. More triangles → more edges
resist polarization → lower biphasage. The relationship is mechanistically direct: it is
the competition term C_protect × tri(e) in the gradient that controls the equilibrium.

## What this means for TCGE

- TCGE does **not** select d = 4
- TCGE biphasage depends on **local clustering**, not dimension
- Dimension is a proxy for triangle density at a given mean degree
- This is a structural law, not a cosmological selection mechanism

## Caveat

The initial 6-point correlation (Δ vs tri/diam, r = 0.986) was a spurious fit
on too few points. The robustness test (15 points, varying both d and ⟨k⟩)
revealed the true driver and killed the tri/diam hypothesis.
