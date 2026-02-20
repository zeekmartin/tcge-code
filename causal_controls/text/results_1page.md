# Causal Identification of Lorentzian-like Phase Separation in Constraint Networks

## Setup

We test whether the phase separation between proto-temporal and proto-spatial edges is caused by (a) local cohesive structure (triangles), (b) the coherence-protection term in the cost, or (c) their interaction, using a 2×2 factorial design with 25 independent seeds per cell.

Each edge weight w ∈ [−1,1] evolves under: E = Σ_e [−A·w² + β·w² + C·tri(e)·w²]. The first term drives polarization (|w|→1); the protection potential P(e) = C·tri(e) resists it locally. An edge is protected (proto-spatial) when P(e) > A−β.

**Structured graphs** (Triangles ↑): Watts-Strogatz (N=80, k=8, p=0.1), ⟨tri⟩ ≈ 3.3.  
**Random graphs** (Triangles ↓): Erdős-Rényi at matched ⟨k⟩, ⟨tri⟩ ≈ 0.8.  
**Protector ON**: C = 0.5. **Protector OFF**: C = 0.

## Result 1: Double Dissociation (Figure 1)

The fraction of proto-spatial edges f_S (defined as |α| < 0.5) shows a clean interaction:

|                    | Protector ON        | Protector OFF       |
|--------------------|---------------------|---------------------|
| **Triangles ↑**    | f_S = 0.85 ± 0.005 | f_S = 0.00 ± 0.000 |
| **Triangles ↓**    | f_S = 0.17 ± 0.008 | f_S = 0.00 ± 0.000 |

(Mean ± SEM, n = 25 seeds per cell)

Two-way ANOVA: interaction F(1,96) = 4497, p < 10⁻⁶⁰, η² = 0.23. The residual variance is η² = 0.005, meaning >99% of the f_S variance is explained by the experimental design.

Neither factor alone produces phase separation. Without the protector, all edges polarize to ⟨|α|⟩ = 1.00 regardless of triangle structure. Without triangles, the protector has nothing to shield (f_S drops to 0.17). Only their conjunction creates the heterogeneous polarization landscape that constitutes the Lorentzian-like signature.

## Result 2: Dose-Response and Triangle Sweep (Figure 2)

**Dose-response (C_protect).** Sweeping C from 0 to 1.5 at fixed ⟨tri⟩ ≈ 3.3 reveals a sharp threshold: f_S = 0 for C < 0.15, onset at C ≈ 0.20, saturation at C ≈ 0.50. The onset is analytically predicted: C_crit = (A−β)/⟨tri⟩ ≈ 0.30 for edges with ⟨tri⟩ = 3.

**Triangle sweep.** Varying clustering (p_ws from 0 to 0.9) at fixed C = 0.5: Protector ON shows f_S increasing monotonically from 0.11 (⟨tri⟩ = 0.6) to 1.00 (⟨tri⟩ = 4.5). Protector OFF: f_S = 0.00 at all clustering levels—a flat line.

This confirms the causal chain: triangles carry the structural information; the protector transduces it into polarization heterogeneity.

## Result 3: Mechanistic Ablation via P(e) (Figure 3)

The protection potential P(e) = C·tri(e) serves as the intermediate variable linking graph structure to polarization:

- **Tri↑ + Prot ON (A):** Strong negative gradient P(e) → |α(e)| (r = −0.79). Edges with P(e) > 1 (i.e., tri ≥ 2 at C = 0.5) → |α| ≈ 0 (spatial). Edges with P(e) < 0.5 → |α| ≈ 1 (temporal). The threshold P(e) = A−β ≈ 0.99 cleanly separates the two populations.

- **Tri↓ + Prot ON (B):** Same local gradient (r = −0.85), but concentrated in the low-P(e) regime—ER graphs lack edges with P(e) > threshold.

- **Prot OFF (C, D):** P(e) = 0 for all edges. |α| = 1.00 uniformly. Zero variance. The triangle structure is present in the graph but invisible to the cost functional.

The Lorentzian-like separation IS the gradient P(e) → |α(e)|: it maps local compatibility structure (triangles = mutual coherence) to metric identity (low |α| = symmetric = spatial; high |α| = asymmetric = temporal), without injecting any T/S label.

## Limits and Scope

1. **The protector is motivated but not uniquely derived.** The term C·tri(e)·w² follows from the TCGE principle that breaking coherent triplets should carry higher cost—but alternative functional forms cannot be excluded. The double dissociation establishes that the mechanism works; it does not prove uniqueness.

2. **Scale.** Results use N = 80 toy models. Previous scaling tests (N = 30–300) confirm qualitative persistence, but continuum-limit behavior remains open.

3. **f_S depends on the |α| < 0.5 threshold.** The qualitative pattern (Tri↑+ON ≫ Tri↓+ON > OFF = 0) is robust to threshold choice; exact values shift.

4. **"Lorentzian-like" is descriptive.** We observe a separation between symmetric (proto-spatial) and polarized (proto-temporal) edge populations. Whether this maps to a Lorentzian metric in a continuum limit requires further work (see Hausdorff dimension analysis in companion results).

5. **AI collaboration.** Mathematical formalization and computational implementation were performed in collaboration with Claude (Anthropic). Conceptual direction, experimental design, and scientific conclusions are the author's.
