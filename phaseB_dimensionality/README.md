# TCGE Phase B — Dimensional Analysis (Exploratory)

**Status:** Exploratory. Not part of v1.0 claims.

## Scripts

| Script | What it does |
|--------|-------------|
| `dim_scan.py` | Scan d=2→7 at fixed N, plus scaled-N continuum diagnostics for d=2,3,4 |
| `dim_robustness.py` | Vary ⟨k⟩ at fixed d to test whether the dimensional sweet spot is real |

| `dim_rewire_inverse.py` | Inverse rewiring: ER → add triangles → watch T/S emerge |
| `metric_S.py` | Mixture-model separation score S (independent of tri classifier) |

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
- The triangle protector's role is not to limit polarization, but to **create two distinct
  populations**: without it, all edges polarize to |α| ≈ 0.98 (one population), and the
  proto-spatial class vanishes entirely (π₁ → 0 in mixture model)
- The protector is therefore **constitutive** of emergent spacetime structure: it is the
  mechanism that produces the T/S distinction, not a regularizer of a pre-existing one

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

## Tri-independent separation metric S (metric_S.py)

The classical biphasage Δ uses tri(e) as both protector AND classifier. To measure
separation independently, we fit a 2-component Beta mixture model on |α|:

```
S = |μ₂ - μ₁| × 2·min(π₁, 1-π₁)
```

where μ₁, μ₂ are component means and π₁ is the weight of the low-|α| component.

**Validation:** S correlates with Δ at r = 0.947 on standard graphs (ER + RGG).

**Rewiring result:** S drops FASTER than Δ — and this is correct. The mixture model
reveals that the proto-spatial population (low |α|) physically vanishes (π₁: 22% → 0%)
rather than just becoming unmeasurable. The protector doesn't just enable measurement
of two populations — it creates them.

| rewire% | π₁ (proto-spatial) | S | Δ |
|---------|-------------------|-------|-------|
| 0% | 22% | 0.112 | 0.369 |
| 10% | 18% | 0.157 | 0.502 |
| 30% | 7% | 0.038 | 0.299 |
| 90% | 0% | 0.001 | 0.139 |

This corrects the earlier histogram-based bimodality test that suggested
separation persists to 50% rewiring — it was an artifact of valley-detection
in a nearly unimodal distribution.

## Trajectory of understanding

1. **Initial hypothesis** (dim_scan.py): Δ ~ tri/diam, r = 0.986 on 6 points
2. **Falsified** (dim_robustness.py): varying ⟨k⟩ reverses the sign
3. **Refined hypothesis**: Δ ~ −⟨tri⟩, r = −0.98 on 15 points
4. **Deepened** (dim_causal.py): tri doesn't control polarization magnitude
   but its spatial heterogeneity — the protector creates the T/S structure
5. **Confirmed** (metric_S.py): tri-independent mixture model shows the
   proto-spatial population physically vanishes without triangles (π₁→0),
   not just that the measurement fails

Each step killed the previous interpretation and replaced it with
something more fundamental. This is documented honestly as a model
of scientific self-correction.

## Inverse rewiring: mirror proof (dim_rewire_inverse.py)

Starting from ER (⟨tri⟩ ≈ 0), triangle-closing rewiring progressively
adds triangles while preserving degrees exactly. The separation score S
(mixture-model based, independent of tri classifier) rises from 0.003 to 0.173:

| ⟨tri⟩ | S | ⟨\|α\|⟩ | μ₁ | μ₂ | Interpretation |
|-------|------|---------|------|------|----------------|
| 0.02 | 0.003 | 0.987 | 0.83 | 0.98 | one population, uniform polarization |
| 0.30 | 0.013 | 0.955 | 0.02 | 0.68 | first proto-spatial nuclei |
| 1.04 | 0.148 | 0.857 | 0.07 | 0.90 | two distinct populations |
| 1.27 | 0.173 | 0.821 | 0.06 | 0.90 | stable T/S separation |

The proto-spatial population (μ₁ ≈ 0.06) did not exist in the ER graph.
It was **created** by the triangles.

**Symmetry with destruction experiment:**
- Destroying triangles (RGG → rewire): S drops, 2 populations → 1
- Creating triangles (ER → tri-close): S rises, 1 population → 2

## Separation metric S (metric_S.py)

S = |μ₂ − μ₁| × 2·min(π₁, 1−π₁) from a 2-component Beta mixture on |α|.
Validated against Δ on standard graphs: corr(S, Δ) = 0.947.
S remains interpretable when tri → 0 (where Δ collapses as a measurement artifact).

## Causal controls (dim_controls.py)

Double-dissociation experiment at identical graph size (N=3000, ⟨k⟩=8):

| Condition | ⟨tri⟩ | S | μ₁ | Interpretation |
|-----------|-------|------|------|----------------|
| ER baseline | 0.02 | 0.003 | 0.84 | no separation |
| C1: random rewire (no tri created) | 0.02 | 0.003 | 0.83 | no separation |
| C2: tri-close + protector OFF | 1.03 | 0.002 | 0.87 | triangles present but invisible to cost |
| REF: tri-close + protector ON | 1.03 | **0.144** | **0.07** | two populations emerge |

**Neither triangles alone nor rewiring alone produce separation.**
Only triangles + active protector in the cost functional creates the T/S structure.
This is a double dissociation: the causal chain is
`triangles × protector → heterogeneous polarization → Lorentzian signature`.
