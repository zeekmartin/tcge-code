# TCGE — Combinatorial Reformulation: Research Log

**Author:** David Martin Venti  
**Date:** 19 February 2026  
**Context:** Response to Rafael Sorkin's objection that weighted DAGs are less simple than posets, importing "all the apparatus of the real numbers."

---

## 1. Research Question

**Can TCGE's key results (metric axioms, Lorentzian-like separation, causal foliation, gravity) be derived from a purely combinatorial substrate — no real numbers in the foundations?**

If yes, the real numbers in TCGE are emergent (continuum limit), not primitive. This directly addresses Sorkin's simplicity objection.

**Key framing (developed in §6):** In causal set theory, the physical content is the *causal order* — not any arbitrary partial order. In TCGE, the physical content is the *anisotropic constraint structure* — not any arbitrary hypergraph. Both frameworks require a physically meaningful datum; neither is inherently simpler than the other.

---

## 2. Methodology

We tested four progressively refined approaches on toy models (N = 81–108 atoms, 4D lattices). All computations use only integers and booleans as primitive data types.

### Substrate (common to all versions)
- **Layer 1:** Finite set Ω + hard constraints H_hard ⊆ P(Ω) — a hypergraph. No real numbers.
- **Layer 2:** All quantitative structure derived from graph topology.
- **Layer 3:** Real numbers only appear as averages/limits.

---

## 3. Results by Version

### v1: Bare degree (w_sym = 1, Δw = N(i) − N(j))

| Test | Result | Status |
|------|--------|--------|
| Metric axioms | 100% (graph distance) | ✅ |
| Lorentzian-like separation | 0/20 seeds | ❌ |
| Causal foliation | 100% (16 N-classes) | ✅ |

**Diagnosis:** With w_sym = 1 for all pairs, s² = 1 − Δw² is negative for almost all pairs (spatial AND temporal) whenever |ΔN| > 1. No directional separation.

### v2: Enriched symmetric weight (w_sym = |common_neighbors| + 1)

Three variants tested:

| w_sym definition | Temporal s² < 0 | Spatial s² > 0 | Lorentzian? |
|-----------------|-----------------|----------------|-------------|
| common_neighbors + 1 | 5.6% | 84.0% | ❌ |
| min(N(i), N(j)) | 0.0% | 99.7% | ❌ |
| (N(i) + N(j)) // 2 | 0.0% | 100.0% | ❌ |

**Diagnosis:** Spatial classification works well. Temporal fails because spatial and temporal neighbors share similar structural properties — common neighbor count doesn't distinguish directions.

### v3: Gradient of N(i) as temporal direction

| Test | Result | Status |
|------|--------|--------|
| Temporal direction identification | 1/20 (5%) | ❌ |
| Lorentzian-like separation | 0/20 | ❌ |
| Causal foliation | 20/20 | ✅ |

**Diagnosis:** N(i) gradient does not reliably identify the temporal direction in isotropic lattices. The gradient is noisy and depends on random constraint placement, not on systematic directional structure.

### v4: Min-cut (Menger) weights on anisotropic hypergraph

| Test | Result | Status |
|------|--------|--------|
| Min-cut anisotropy | 0/10 seeds | ❌ |
| Lorentzian-like separation | 0/10 | ❌ |
| Causal foliation | 10/10 | ✅ |

Anisotropy phase diagram (spatial rate fixed at 5%):

| Temporal constraint rate | mc_temporal | mc_spatial | Ratio |
|--------------------------|-------------|------------|-------|
| 0% | 4.96 | 4.90 | 1.01 |
| 20% | 4.65 | 4.60 | 1.01 |
| 40% | 4.30 | 4.19 | 1.02 |
| 60% | 3.98 | 3.92 | 1.01 |

**Diagnosis:** Min-cut is a GLOBAL connectivity measure. Even with 44% of temporal edges removed, nodes remain well-connected through spatial detour paths. The min-cut ratio stays near 1.0 regardless of anisotropy. Min-cut does not capture directional tension.

---

## 4. What Works Without Real Numbers (robust across all versions)

| Result | Mechanism | Robustness |
|--------|-----------|------------|
| **Metric axioms** | Graph distance = integer, satisfies axioms by construction | 100% (all versions, all seeds) |
| **Causal foliation** | N(i) = integer degree → temporal classes with inter-class ordering | 100% (all versions, all seeds) |
| **Acyclicity** | Inter-class order is DAG; cycles confined to N-degenerate classes | 100% |

These three results are **trivially combinatorial** and require no real numbers whatsoever.

## 5. What Does NOT Work Without Real Numbers

| Result | Difficulty | Root cause |
|--------|-----------|------------|
| **Lorentzian-like separation** | 0% across all versions (v1-v4) | See analysis below |
| **Gravity (geodesic deviation)** | Insufficient data | Requires Lorentzian first |
| **Einstein equations** | Not tested | Requires continuous approximation |

### Root Cause Analysis: Why Lorentzian-like Separation Fails (v1–v4)

The Lorentzian signature (−,+,+,+) requires that the signed interval s² be **systematically negative** along one direction and **systematically positive** along the other three. This demands a clean **directional asymmetry** in some quantitative measure.

**The fundamental problem:** In a compatibility graph derived from hard constraints alone, the combinatorial quantities (degree, common neighbors, min-cut, graph distance) are all **isotropic by construction** on a regular lattice, even when the constraint density is anisotropic. This is because:

1. **Degree N(i)** aggregates all directions — it doesn't know which direction the constraints come from.
2. **Common neighbors** are shared across directions.
3. **Min-cut** uses all available paths (including cross-directional ones).
4. **Graph distance** is shortest-path, which routes around removed edges.

The constraint asymmetry gets **washed out** by the global connectivity of the graph. Removing 40% of temporal edges is invisible to global measures because spatial paths provide abundant alternatives.

### Comparison with Standard TCGE (real-valued weights)

In standard TCGE, the Lorentzian signature emerges because **directed weights** w(a→b) ≠ w(b→a) carry directional information at the level of each individual edge. This per-edge asymmetry is the key ingredient. It is fundamentally a **local, directional** datum — not a global graph property.

A purely combinatorial replacement would need a **local, directional, integer-valued** quantity. Candidates not yet tested:

- **Directional edge density** in a local neighborhood (count of temporal vs spatial edges around node i)
- **Directed min-cut** on a DAG derived from the N(i) ordering
- **Path asymmetry:** number of paths from i to j through temporal edges only vs spatial edges only
- **Local constraint count** per direction (how many H_hard hyperedges involve edges in direction μ near node i)

---

## 6. Theoretical Interpretation

### What Sorkin's Objection Really Asks

Sorkin's objection is not merely about real numbers. It is about whether the **directed, weighted structure** of TCGE can be reduced to something simpler. Our results show:

- **3 of 4 key results** (metric, foliation, acyclicity) survive reduction to pure combinatorics.
- **The 4th result** (Lorentzian-like separation) requires per-edge directional information that global combinatorial measures cannot provide.

### The Structural Analogy

This is actually symmetric with Sorkin's own situation:

| | Causal Sets (Sorkin) | TCGE Combinatorial |
|---|---|---|
| **Substrate** | Partial order (poset) | Constraint hypergraph |
| **Physical content** | The poset is *causal* (not any partial order) | The hypergraph is *anisotropic* (not any hypergraph) |
| **What's "extra"** | Causality structure | Directional constraint density |
| **What emerges** | Geometry (Order + Number) | Geometry + Dynamics (Coherence minimization) |

In both cases, the **physical content** — what distinguishes our universe from a generic mathematical structure — must be supplied. Sorkin supplies it as "causal order." TCGE supplies it as "anisotropic constraint structure." Neither is more or less fundamental than the other.

### The Honest Answer to Sorkin

> The real numbers in TCGE are not ontologically primitive. Three of TCGE's four key results (metric geometry, causal foliation, temporal acyclicity) can be derived from a purely combinatorial substrate using integer-valued quantities. The fourth result (Lorentzian-like separation) requires anisotropic constraint density — specifically, that constraints are distributed differently along one direction vs the others. This is the physical content of the theory, analogous to the requirement that a causal set be "causal" rather than an arbitrary partial order. The weights are the simplest encoding of this anisotropic information. Whether they can be further reduced to a purely order-theoretic datum remains an open question — and a well-defined research programme.
>
> **Update (post v5b–v5d):** Subsequent experiments confirmed that Lorentzian-like separation emerges from boolean edge markings (constrained / unconstrained) with anisotropic density (Cohen's d = 7.5 vs permutation null, p < 0.001). Constraint *orientation* does not contribute to this separation (d = 0.00 in knife-edge test). The arrow of time remains an open problem requiring additional structure.

---

## 7. Open Questions (Future Work)

1. **Local directional quantities:** Can per-direction edge counts or constraint densities in local neighborhoods provide the missing directional information? (Promising: this is still integer-valued.)

2. **Continuum limit:** Do the combinatorial results (metric, foliation) converge to smooth structures as N → ∞? (Requires spectral dimension analysis, Gromov-Hausdorff convergence tests.)

3. **Dynamics from combinatorial min-cut:** Can Poisson/Einstein equations be recovered from optimization on a graph with min-cut weights? (Requires larger lattices.)

4. **Hypergraph-native formulation:** Can n-ary constraints (not just pairwise) provide directional information that pairwise measures miss?

---

## 8. Files

| File | Description |
|------|-------------|
| `tcge_combinatorial_reformulation.py` | v1: bare degree |
| `tcge_combinatorial_v2.py` | v2: enriched symmetric weights |
| `tcge_combinatorial_v3.py` | v3: N(i) gradient approach |
| `tcge_mincut_v4.py` | v4: min-cut (Menger) weights |
| `tcge_directed_v5.py` | v5: directed boolean constraints |
| `tcge_combinatorial_results.json` | v1 output |
| `tcge_mincut_results.json` | v4 output |
| `tcge_directed_v5_results.json` | v5 output |
| `RESEARCH_LOG.md` | This document |

---

## ADDENDUM: v5 — Directed Boolean Constraints (Breakthrough)

### Approach

Replace real-valued weights with **boolean directed constraints**: for some compatible pairs (i,j), a relation directed[i,j] = True means "i constrains j." This is the physical content — analogous to causal order in causal sets. Data types: booleans only.

The directed constraints are anisotropic: ~70% of temporal edges carry a directed constraint, vs ~5-10% of spatial edges. This is the physical datum, not a parameter.

### Key Innovation: Per-edge Asymmetry Δ_v3

Three versions of Δ(i,j) were tested:

| Δ version | Definition | Temporal/Spatial ratio | Directional? |
|-----------|-----------|:---:|:---:|
| v1: global out-degree | out(i) - out(j) | 0.98 | ❌ |
| v2: local neighborhood | directed influence on neighborhoods | 0.00 (degenerate) | ❌ |
| **v3: per-edge** | **directed[i,j] − directed[j,i] ∈ {−1,0,+1}** | **8–48×** | **✅** |

The per-edge measure is massively directional because it is the **most local** measure possible. It doesn't aggregate across directions. This is exactly what v1-v4 were missing.

### Results

| Test | Result | Status |
|------|--------|--------|
| Δ directional | 20/20 (100%) | ✅ |
| Lorentzian-like separation | 6/20 (30%) | ⚠️ Partial |
| Causal foliation | 20/20 (100%) | ✅ |

Signature detail at best threshold (|Δ| ≤ 1, λ=2):
- Temporal edges: 88.9% timelike ✅
- Spatial edges: only 30.8% spacelike ⚠️ (most also negative)

### Comparison Across All Versions

| Version | Primitive data | Δ directional | Lorentzian | Foliation |
|---------|---------------|:---:|:---:|:---:|
| v1 (degree) | integers | — | 0/20 (0%) | 20/20 |
| v2 (common neighbors) | integers | — | 0/20 (0%) | 20/20 |
| v3 (gradient) | integers | 1/20 (5%) | 0/20 (0%) | 20/20 |
| v4 (min-cut) | integers | 0/10 (0%) | 0/10 (0%) | 10/10 |
| **v5 (directed boolean)** | **booleans** | **20/20 (100%)** | **6/20 (30%)** | **20/20 (100%)** |

### Diagnosis of Remaining Issue

The spatial classification fails because the temporal distance t(i,j) = |rank(i) − rank(j)| is nonzero for some spatial neighbors (they differ in directed-graph rank even though they are spatially adjacent). This contaminates the spatial s² values. Possible fixes:

1. Use a more local temporal measure (e.g. per-edge Δ instead of global rank)
2. Better G_sym filtering (exclude edges with ANY directed constraint)
3. Two-step procedure: first identify temporal direction from Δ, then compute spatial distance only within temporal slices

### Theoretical Significance

This is the first version where:
- **Directional asymmetry is robust** (100% detection rate)
- **Lorentzian-like separation partially emerges** (30% vs 0% in all previous attempts)
- **The primitive data type is boolean** (simpler than integers, much simpler than reals)

The directed boolean constraint is structurally comparable to a poset relation (it's a directed graph). But unlike a bare poset, the TCGE compatibility + directed constraint structure carries enough information for metric geometry, causal foliation, AND partial Lorentzian separation.

### Updated Honest Answer to Sorkin

> TCGE's quantitative structure can be reduced to boolean directed relations on a compatibility graph — no real numbers, no integers even. The directed constraints are the physical content, directly analogous to the causal order in causal set theory. From this boolean substrate, we derive: (i) metric geometry from graph distance (100% robust), (ii) causal foliation from directed-graph rank (100% robust), and (iii) partial Lorentzian-like separation where temporal edges are correctly classified as timelike 89% of the time. The remaining challenge — clean spatial classification — is a well-defined technical problem, not a fundamental obstacle. The real numbers in TCGE's published formulation are the continuum shadow of this boolean directed structure.

---

*This research was conducted in collaboration with Claude (Anthropic), used as a computational and analytical tool. All conceptual direction, evaluation, and conclusions are the author's responsibility.*

---

## ADDENDUM 2: v5c — Statistical Validation (Permutation Tests)

### Results

| Comparison | Cohen's d | p-value |
|-----------|:---------:|:-------:|
| Real vs Global permutation (score) | **+7.51** | **< 0.001** |
| Real vs Global permutation (D statistic) | **−1.18** | — |
| Real vs Local permutation (score) | **−0.09** | 0.70 |

At λ = 10: **14/15 real** pass Lorentzian vs **0/150 permuted** (global). Distributions completely non-overlapping.

### Critical Discovery

The **local permutation** (keep which pairs have constraints, randomize direction) gives **identical** scores to the real signal (d = −0.09). This means the signal comes from **constraint density anisotropy**, not from **constraint orientation**.

**Note on λ calibration:** λ = 10 was selected by grid search over λ ∈ {1, ..., 15}, choosing the value that maximizes Youden's J (TPR − FPR) against the global permutation null. The effect is not sensitive to the exact value: all λ ≥ 10 produce identical results (score plateaus), and the permutation null gives score ≈ 0 across the entire range. This rules out parameter tuning as an explanation for the signal.

---

## ADDENDUM 3: v5d — Knife-Edge Test (Definitive Discrimination)

### Protocol

Three conditions with EQUAL constraint density (p = 0.3, 0.5, 0.7) for temporal and spatial edges:
- **A**: Structured orientation (aligned with lattice time)
- **B**: Random orientation (same support map, direction randomized)
- **C**: Shuffled support (constraint assignments randomized globally)

### Results (all densities, all λ values)

| Comparison | Cohen's d | Interpretation |
|-----------|:---------:|:-------------|
| A vs B (orientation effect) | **0.00** | Orientation contributes **nothing** |
| A vs C (support effect) | −0.32 to −0.78 | Support map matters (when unequal) |
| Anisotropic A vs Anisotropic B | **0.00** | Even with 70%/5% density asymmetry, direction irrelevant |

Arrow alignment: Structured = 1.00, Random = 0.49 — yet identical Lorentzian scores.

### Definitive Conclusion

**The Lorentzian-like separation in TCGE toy models emerges from ONE primitive:**

> **Anisotropic constraint density** — which edges carry constraints, not which direction they point.

This is even simpler than "boolean directed constraints." The minimal primitive is a **binary edge marking** (constrained / unconstrained) with anisotropic distribution.

### Two Separable GAPs

| GAP | Question | Status |
|-----|----------|--------|
| **GAP-Signature** | Can Lorentzian-like separation emerge without ℝ? | **CLOSED ✅** (d = 7.5, p < 0.001) |
| **GAP-Arrow** | Can past→future emerge without postulate? | **OPEN ⚠️** (d = 0.00, orientation irrelevant in current models) |

**Scope:** GAP-Signature refers to combinatorial timelike/spacelike classification on toy lattices, not to full Lorentzian metric structure (local Minkowski invariance, smooth limit). Orientation is irrelevant for this separation in current models; the arrow of time requires additional structure.

### Comparison Table (Complete Series)

| Version | Primitive | Lorentzian | Directional | Foliation |
|---------|-----------|:---:|:---:|:---:|
| v1 (degree) | integer | 0% | — | 100% |
| v2 (common neighbors) | integer | 0% | — | 100% |
| v3 (gradient) | integer | 0% | 5% | 100% |
| v4 (min-cut) | integer | 0% | 0% | 100% |
| v5 (directed boolean) | boolean directed | 30% | 100% | 100% |
| v5b (slice separation) | boolean directed | **95%** | 100% | 100% |
| v5c (validated) | boolean directed | **93%** vs **0%** null | — | 100% |
| **v5d (knife-edge)** | **binary edge marking** | **Signal from density only** | — | — |

### Files (Complete)

| File | Description |
|------|-------------|
| `tcge_combinatorial_reformulation.py` | v1: bare degree |
| `tcge_combinatorial_v2.py` | v2: enriched symmetric weights |
| `tcge_combinatorial_v3.py` | v3: N(i) gradient approach |
| `tcge_mincut_v4.py` | v4: min-cut (Menger) weights |
| `tcge_directed_v5.py` | v5: directed boolean constraints |
| `tcge_directed_v5b.py` | v5b: slice-based separation |
| `tcge_v5c_validation.py` | v5c: permutation tests + effect sizes |
| `tcge_v5d_knife_edge.py` | v5d: knife-edge density vs orientation |
| `RESEARCH_LOG.md` | This document |

---

*This research was conducted in collaboration with Claude (Anthropic), used as a computational and analytical tool. All conceptual direction, evaluation, and conclusions are the author's responsibility.*
