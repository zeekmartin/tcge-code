# Why Triangle Protection? A Derivation from Minimal Coherence Principles

## The Question

The computational results establish that the Lorentzian-like phase separation requires a cost term of the form P(e) = C · tri(e) · w². A reviewer trained in foundations will ask: *Why this term? Is it ad hoc or derivable?*

We show it follows from three requirements that any coherence-based cost functional must satisfy.

---

## Axiom 1: Locality

The cost functional decomposes as a sum over edges:

    E = Σ_e  f(w_e, local_structure(e))

This is the standard assumption in lattice field theory, spin models, and causal set dynamics. Non-local terms (depending on distant edges) would violate the principle that constraints propagate through adjacency.

## Axiom 2: Symmetry-breaking drive

There exists a term in f that favors |w_e| → W (maximal polarization). In our implementation this is −A·w². In TCGE, this arises from the product structure of the cost: the incompatibility between adjacent nodes is minimized when edge weights are extremal. This is not controversial — it is the analogue of the double-well potential in Landau theory, or the kinetic term in a lattice gauge theory that favors alignment.

## Axiom 3: Coherence penalty

The cost should penalize the breaking of *local coherence*. The question is: what constitutes local coherence, and what is its minimal expression?

### Definition: Local coherence at edge e

An edge e = (i,j) is *locally coherent* if its polarization is compatible with the constraints seen by its immediate neighborhood. The minimal unit of mutual constraint beyond the edge itself is a **triangle**: if i−j−k form a triangle, then the constraints i→j, j→k, k→i form a closed cycle. Polarizing edge (i,j) (setting |w_{ij}| → W) while the cycle constraints demand consistency creates a tension.

### Claim: The triangle is the minimal coherence unit

**Proof sketch.** Consider what local structures an edge participates in:

1. **The edge itself** — provides no context for coherence (it's the object being polarized).

2. **A path of length 2 through e** (i−j−k, no closing edge) — the constraints i→j and j→k are independent. Polarizing (i,j) creates no tension with (j,k) because there is no closure condition. No coherence to break.

3. **A triangle through e** (i−j−k−i) — now the three constraints form a cycle. If (i,j) polarizes but (j,k) and (k,i) remain symmetric, the cycle is frustrated: the polarization of one edge is inconsistent with the compatibility structure of the other two. This IS coherence.

4. **Larger cycles** (4-cliques, 5-cliques, ...) — these contain triangles as substructures. Their coherence contribution decomposes into triangular units (by inclusion-exclusion or by the flag complex structure).

Therefore: **the triangle is the minimal non-trivial unit of local coherence.** Any cost functional that penalizes coherence-breaking must, at leading order, depend on the triangle count.

### From coherence to cost term

The penalty for breaking coherence at edge e should be:

1. **Proportional to the amount of coherence present** — measured by tri(e), the number of triangles containing e. More triangles = more mutual constraints to violate.

2. **Proportional to the degree of polarization** — measured by some increasing function of |w_e|. An unpolarized edge (|w_e| ≈ 0) breaks no coherence regardless of its neighborhood.

3. **Positive** — it is a penalty, not a reward.

The simplest term satisfying (1)-(3) is:

    P(e) = C · tri(e) · g(w_e)

where g is an even, increasing function of |w| with g(0) = 0. The lowest-order choice is g(w) = w². Higher orders (w⁴, etc.) work but are less canonical.

### Robustness to functional form

The computational tests confirm this: replacing w² by w⁴, or tri(e) by √tri(e), changes the amplitude of phase separation but preserves its existence and the interaction structure. Only the binary threshold (tri(e) > 2 → 1) fails, because it discards the graded coupling between coherence and protection strength.

This is the expected pattern: **the phenomenon requires graded coherence coupling, not a specific functional form.** The triangle count is the essential structural variable; the coupling function modulates amplitude.

---

## The Protection Threshold

The effective potential per edge is:

    V(w_e) = (−A + β + C·tri(e)) · w_e²

Edge e is protected (stays near w = 0, proto-spatial) when:

    C · tri(e) > A − β

and polarizes (|w| → W, proto-temporal) when:

    C · tri(e) < A − β

This gives a local, analytically predictable threshold. The Lorentzian-like separation is the *set of edges above vs. below threshold* — determined entirely by the local triangle count.

### Why this is not ad hoc

The threshold emerges from the competition between two terms (drive and protection) that are independently motivated:

- The drive term (−A·w²) comes from the TCGE minimization principle: incompatibility is reduced by extremal weights.
- The protection term (C·tri(e)·w²) comes from Axiom 3: breaking local coherence has a cost proportional to the coherence present.

Neither term was designed to produce a threshold. The threshold is an *emergent property* of their competition. The 2×2 factorial experiment confirms this: removing either term eliminates the phase separation.

---

## Addressing Potential Objections

**"Why not penalize 4-cliques, 5-cliques, etc.?"**

These are higher-order corrections. The triangle is the leading term because:
(a) It is the minimal coherence unit (shown above).
(b) Higher cliques are rarer and their contribution is absorbed by the triangle term at leading order.
(c) Computationally, adding k-clique terms (k > 3) does not change the qualitative result — only the amplitude.

**"In causal set theory, the causal structure doesn't use triangles."**

Correct. Causal sets use the partial order directly. TCGE is not causal set theory. In TCGE, the constraint network is pre-geometric: there is no partial order yet. The triangle structure is the *source* from which causal structure will emerge, not a consequence of it. The analogy: in CDT, the triangulation rules are pre-geometric; the metric emerges from the ensemble. In TCGE, the triangle density is pre-geometric; the Lorentzian-like separation emerges from the cost functional.

**"Is C a free parameter?"**

C controls the strength of coherence protection relative to the polarization drive. It is not a fundamental constant but a parameter of the cost functional, analogous to a coupling constant in lattice field theory. The dose-response curve shows that the qualitative phenomenon exists for all C above a threshold (C_crit ≈ 0.3 for typical graphs), and saturates above C ≈ 0.5. The phase separation is structurally stable over a wide range of C.

**"What about non-triangular coherence (e.g., algebraic constraints)?"**

In a more general TCGE framework, coherence could be measured by higher-order compatibility structures. The triangle is the graph-theoretic proxy. If the constraint network has algebraic structure (e.g., gauge symmetry), the coherence measure would be the gauge-invariant analogue. This is a direction for future work — the current results establish the mechanism at the graph-theoretic level.

---

## Summary

The triangle protector is not a design choice — it is the minimal local coherence penalty consistent with three axioms (locality, symmetry-breaking drive, coherence penalty). The triangle is the smallest substructure where coherence is defined. The cost term P(e) = C · tri(e) · w² is the leading-order expression of this penalty. The Lorentzian-like phase separation emerges from the competition between this term and the polarization drive, producing a local threshold that separates proto-spatial (protected) from proto-temporal (polarized) edges without any injected label.
