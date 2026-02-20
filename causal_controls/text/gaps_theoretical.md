# GAP 2: From "Lorentzian-like" to Lorentzian — What Remains

## What we have established

The computational results demonstrate a binary separation of edge populations:

| Population | |α| | Role | Properties |
|---|---|---|---|
| Proto-spatial | ≈ 0 | Symmetric | High local coherence (tri↑), resist polarization |
| Proto-temporal | ≈ 1 | Polarized | Low coherence (tri↓), direction-broken |

This is "Lorentzian-like" in a precise morphological sense: it mirrors the distinction between spacelike (symmetric under reversal) and timelike (direction-privileged) in a Lorentzian manifold.

## What a Lorentzian signature actually requires

A true Lorentzian metric signature (−,+,+,...) implies four things:

1. **Two distinct classes of directions** — ✅ We have this. The bifurcation theorem guarantees two populations on any graph with heterogeneous tri(e).

2. **Dimensional ratio** — The number of "spatial" directions should exceed the number of "temporal" directions. In 3+1D: 3 spatial for 1 temporal. 
   
   **Status: Partially addressed.** In our WS models, f_S ≈ 0.85, giving roughly 85% spatial / 15% temporal. The exact ratio depends on the triangle distribution of the substrate, not on a postulated dimension. Whether specific graph ensembles reproduce the 3:1 ratio is an open question.

3. **Causal cone structure** — Timelike-connected regions should form cones (future/past), not arbitrary subsets.

   **Status: Open, but groundwork exists.** The arrow analysis (v4, 91% directional coherence) shows that polarized edges carry a consistent direction. If proto-temporal edges form connected subgraphs with coherent directionality, these are proto-cones. The key test would be: does the subgraph induced by high-|α| edges have small diameter (cone-like) or is it scattered?

4. **Indefinite inner product** — The metric should have signature (−,+,...,+), not just "two types of things."

   **Status: Open.** Translating the discrete edge variable |α| into a continuous metric tensor g_μν requires a coarse-graining procedure (averaging over neighborhoods). The natural candidate is:

       g_eff(x) ~ Σ_{e near x}  sign(e) · w_e ⊗ w_e

   where sign(e) = −1 for proto-temporal edges and +1 for proto-spatial. If the coarse-grained tensor has eigenvalues of mixed sign, we have a Lorentzian signature. This is computable but requires a continuum embedding — it belongs to the next phase of work.

## The honest statement

We can say: *"The phase separation has the structure of a Lorentzian-like signature: two populations with qualitatively different symmetry properties (symmetric vs. direction-broken), separated by a local threshold determined by the graph's coherence structure. Whether this maps to a Lorentzian metric tensor in a coarse-grained continuum limit is an open question that requires (a) embedding into a geometric substrate and (b) computing the effective metric."*

This is defensible, precise, and doesn't overclaim.

## Path to closing this gap

Three concrete steps, in order of difficulty:

**Step 1 (easy):** On geometric substrates (RGG in d=4), measure the subgraph of proto-temporal edges. Is it connected? What is its dimension? If it forms 1D filaments, those are proto-worldlines.

**Step 2 (medium):** Compute the coarse-grained metric. Partition the graph into regions. For each region, compute the average α-weighted tensor. Check if eigenvalues have mixed sign.

**Step 3 (hard):** Compute geodesics on the effective metric. Do null geodesics form cones?

---

# GAP 4: From Discrete Threshold to Continuum Signature — The Bridge

## The discrete result

On a graph, each edge has a protection potential P(e) = C·tri(e). The bifurcation occurs at P(e) = A − β. Edges above threshold: protected (spatial). Below: polarized (temporal).

This is a **local, edge-by-edge** classification. There is no global structure yet.

## The continuum question

In a Lorentzian manifold, the metric signature is a *field* — it varies smoothly (or at least continuously) over the manifold. How does a discrete edge-by-edge threshold produce a smooth signature field?

## The Landau analogy

The closest physical analogy is the **Ising model / Landau theory of phase transitions**:

| Ising/Landau | TCGE |
|---|---|
| Spin s_i ∈ {−1, +1} | Edge weight w_e ∈ [−W, W] |
| Double-well potential V(s) = −as² + bs⁴ | Drive term −A·w² |
| Coupling J·s_i·s_j | Protection C·tri(e)·g(w) |
| Magnetization m = ⟨s⟩ | Polarization α = ⟨|w|/W⟩ |
| Ordered (T < T_c) | Protected (P(e) > A−β) |
| Disordered (T > T_c) | Polarized (P(e) < A−β) |

In the Ising model, the discrete spin-by-spin configuration produces a *smooth magnetization field* m(x) in the continuum limit, via coarse-graining over mesoscopic blocks. The Landau free energy F[m] = ∫ (−a·m² + b·m⁴ + c·(∇m)²) dx is the continuum theory.

## The TCGE analogue

By analogy, we expect the continuum limit of the edge polarization to be a **signature field** σ(x):

    σ(x) ~ 1/|B_r(x)| · Σ_{e ∈ B_r(x)}  (classification of e)

where B_r(x) is a ball of coarse-graining radius r. The key requirements for this to produce a smooth field:

1. **Spatial coherence** — Nearby edges should have similar |α|. This is guaranteed if the triangle distribution is spatially correlated (which it is on geometric substrates like RGG, where nearby vertices share neighbors).

2. **Scale separation** — The coarse-graining radius r should be much larger than the lattice spacing but much smaller than the system size. This is the standard assumption in any lattice-to-continuum passage.

3. **Stability under coarse-graining** — The classification (spatial vs. temporal) should not flip when averaging over neighborhoods. Our previous tests showed 40–60% retention of the biphasage signal under coarse-graining in d=3, which is partial but non-zero.

## The CDT analogy

In Causal Dynamical Triangulations (CDT), the Lorentzian signature is *built into* the triangulation rules: time-like and space-like edges are labeled from the start. The signature is an input, not an output.

In TCGE, the signature is **emergent**: there are no labels. The analogue of CDT's triangulation rules is the cost functional, and the signature emerges from the competition between drive and protection.

This is a stronger claim (emergence vs. postulation) but comes with a higher burden of proof: we must show that the emergent classification is *consistent* across scales, which CDT gets for free.

## The renormalization intuition (qualitative)

A block-spin renormalization group (RG) step in TCGE would:

1. Partition vertices into blocks of size b
2. Define block-level weights as some aggregate of edge weights
3. Define block-level triangles from the block graph
4. Repeat the optimization

If the phase separation survives, the RG fixed point has a meaningful signature field. Our coarse-graining tests suggest partial survival (~50% retention), which is promising but not conclusive.

A more rigorous approach would be to define a **real-space RG** on the cost functional:

    E_eff[{w_block}] = RG(E[{w_edge}])

and check whether the bifurcation structure is preserved. This is feasible but goes beyond the current work.

## The honest statement

*"The discrete phase separation defines a local classification (spatial vs. temporal) for each edge, driven by the bifurcation threshold P(e) = A−β. On geometric substrates, this classification is spatially coherent because triangle density is spatially correlated. The passage to a continuum signature field requires coarse-graining, which our preliminary tests show preserves 40–60% of the signal. A full continuum limit, analogous to the Landau theory of magnetization, remains an open problem. The key open question is whether a real-space renormalization group preserves the bifurcation structure."*

## Summary of all four gaps

| Gap | Status | What we showed | What remains |
|---|---|---|---|
| **1. Bifurcation universality** | ✅ Closed | Pitchfork bifurcation for any convex g with g''(0)>0. Verified on 6 families. | — |
| **2. Lorentzian nature** | 🟡 Partial | Two populations with correct symmetry properties. Arrow coherence 91%. | Cone structure, effective metric tensor, dimensional ratio. |
| **3. Dynamic stability** | ✅ Closed | Unique basin, corr > 0.94 under ε=0.5 perturbation, 3 initial conditions converge. | — |
| **4. Continuum bridge** | 🟡 Partial | Landau analogy precise. Spatial coherence on RGG. CG retention 40–60%. | Real-space RG on cost functional. Full coarse-grained metric. |
