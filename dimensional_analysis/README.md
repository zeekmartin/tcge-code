# TCGE Phase B — Dimensional Analysis (Exploratory)

**Status:** Exploratory. Not part of v1.0 claims.  
**Sessions:** Conducted in previous sessions (see conversation history).

## Purpose

Test whether the biphasage mechanism preferentially selects d=4 dimensions.

## Scripts

| Script | What it does |
|--------|-------------|
| `dim_scan.py` | Scan d=2→7 at fixed N, continuum diagnostics for d=2,3,4 |
| `dim_robustness.py` | Vary ⟨k⟩ at fixed d to test dimensional sweet spot |
| `dim_causal.py` | Degree-preserving rewiring to isolate triangle effect |

## Key Findings (chronological — each killed the previous)

### 1. Initial hypothesis (dim_scan.py)
Δ ~ tri/diam, r = 0.986 on 6 points → apparent "sweet spot" at d ≈ 4–5.

### 2. Falsified (dim_robustness.py)
Varying ⟨k⟩ at fixed d reverses the correlation. The fit was spurious on too few points.

### 3. Refined hypothesis
Δ ~ −⟨tri⟩ alone, r = −0.98 on 15 conditions (d=3–5, ⟨k⟩=5–20).
Dimension has **no independent effect** (multivariate coefficient ≈ 0).

```
Δ ≈ −0.034 × ⟨tri⟩ + 0.51
r = −0.98, R² = 0.97, n = 15
```

### 4. Deepened (dim_causal.py) — THE key result
Degree-preserving rewiring revealed the deepest mechanism:

| rewire% | ⟨tri⟩ | ⟨|α|⟩ | Δ     |
|---------|-------|-------|-------|
| 0%      | 3.80  | 0.43  | 0.369 |
| 10%     | 2.08  | 0.65  | 0.502 |
| 50%     | 0.21  | 0.95  | 0.211 |
| 90%     | 0.04  | 0.98  | 0.139 |

**The triangle protector does not control the amount of polarization —
it controls its spatial distribution.** Destroying triangles drives ALL
edges to maximal polarization (⟨|α|⟩ → 0.98). The Lorentzian-like phase
separation requires the protector to selectively shield cohesive edges.

## Interpretation

- The "sweet spot d=4–5" at ⟨k⟩=8 was an artefact: those dimensions
  just happen to have ⟨tri⟩ ≈ 2.5–3 at that degree.
- Biphasage is controlled by **local 3-clique abundance**, which measures
  how costly it is to break compatibility. Dimension has no independent role.
- This insight directly motivates the causal controls work (see
  `../causal_controls/`).

## Note on scripts

The original scripts were created in previous Claude sessions. If they
are not present here, they can be reconstructed from the conversation
history or from the GitHub repo `tcge-code/`. The key results and data
are documented above and in the causal_controls package which supersedes
this exploratory analysis.
