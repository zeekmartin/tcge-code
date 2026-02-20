# TCGE Emergence Suite — Results

All results below are reproducible using the scripts in this directory
with the seeds documented in `seeds_used.txt`.

## 1. Phase separation (01_biphasage.py)

| Metric | Value |
|--------|-------|
| Biphasage Δ | 0.31 ± 0.06 |
| |α|_low (proto-temporal) | 0.55 ± 0.04 |
| |α|_high (proto-spatial) | 0.24 ± 0.03 |
| Trials with Δ > 0.2 | 97% (29/30) |
| Graph | Erdős–Rényi, N=100, ⟨k⟩≈10 |

## 2. Robustness (02_robustness.py)

| Cohesion metric | Δ | Pass |
|-----------------|---|------|
| Triangles | 0.31 ± 0.06 | ✅ |
| Jaccard | 0.29 ± 0.05 | ✅ |
| Edge clustering | 0.28 ± 0.06 | ✅ |
| Quadrangles | 0.25 ± 0.05 | ✅ |
| Truss number | 0.27 ± 0.06 | ✅ |
| Betweenness | 0.04 ± 0.03 | ❌ |

Finite-size scaling: Δ stable across N = 50–300 (plateau).

## 3. Directional arrow (03_arrow.py)

| Metric | Value |
|--------|-------|
| Arrow coherence | ~91% |
| Sign correlation with degree gradient | > 0.90 |
| Mechanism | α direction couples to ΔN(e) |

## 4. Continuum diagnostics (04_continuum.py)

RGG on torus, N = 5000, ⟨k⟩ ≈ 8, 5 trials.

| dim | d_H (plateau) | σ(d_H) | d_H (global) | Biphasage Δ |
|-----|---------------|--------|--------------|-------------|
| 3 | 2.84 | 0.03 | 2.63 | 0.376 |
| 4 | 3.38 | 0.02 | 3.24 | 0.400 |

## 5. Universality (05_universality.py)

N = 150, ⟨k⟩ ≈ 12, 10 trials per family, 200 permutations for null.

| Family | Δ | p (perm) | Cohen's d | Arrow | Verdict |
|--------|---|----------|-----------|-------|---------|
| ER | 0.246 | < 0.005 | 11.9 | 0.98 | STRONG |
| RGG 3D | 0.283 | < 0.005 | 9.7 | 1.00 | STRONG |
| WS | 0.192 | < 0.005 | 11.6 | 0.93 | STRONG |
| BA | 0.110 | < 0.005 | 12.2 | 0.97 | PASS |
| CM | 0.250 | < 0.005 | 8.3 | 0.98 | STRONG |

Pass criterion: Δ > 0.05 AND p < 0.05.
Strong criterion: Δ > 0.15 AND p < 0.01 AND Cohen's d > 0.8.
Result: 5/5 PASS, 4/5 STRONG.

Correlation Δ vs var(tri): Pearson r = −0.51 (n = 50).

## 6. Local homogeneity (06_locality.py)

RGG torus, N = 5000.

| Test | d=3 | d=4 |
|------|-----|-----|
| d_H local CV | 0.082 (✅ < 0.1) | — (diameter too short) |
| Distance distribution | Unimodal, mode=11 | Unimodal, mode=8 |
| Isotropy score | 0.82 (✅ > 0.8) | 0.72 (⚠️) |

## 7. Coarse-graining (07_coarsegraining.py)

RGG torus, N = 5000, Louvain k ≈ 200.

| Aggregator × Classifier | d=3 retention | d=4 retention |
|--------------------------|---------------|---------------|
| mean × density (baseline) | 51% | 36% |
| mean × n_edges (best) | **59%** | **42%** |
| max × density | −7% | −11% |
| q90 × density | −1% | −7% |
| trimmed × n_edges | 41% | 30% |

Robust aggregators (max, q90) invert the signal.
The improvement comes from the classifier (n_edges > density), not the aggregator.
Structural limit: CG averaging compresses the fine-scale signal.

---

## CSV format

```csv
experiment,family,N,k,metric,value,std,p_value,cohens_d
biphasage,ER,100,10,delta,0.31,0.06,,
biphasage,ER,100,10,alpha_low,0.55,0.04,,
biphasage,ER,100,10,alpha_high,0.24,0.03,,
universality,ER,150,12,delta,0.246,0.028,0.005,11.86
universality,RGG_3D,150,12,delta,0.283,0.034,0.005,9.74
universality,WS,150,12,delta,0.192,0.012,0.005,11.55
universality,BA,150,12,delta,0.110,0.022,0.005,12.15
universality,CM,150,12,delta,0.250,0.039,0.005,8.26
continuum,RGG_3D,5000,8,dH_plateau,2.84,0.03,,
continuum,RGG_4D,5000,8,dH_plateau,3.38,0.02,,
locality,RGG_3D,5000,8,dH_local_CV,0.082,,,
locality,RGG_3D,5000,8,isotropy,0.82,0.087,,
locality,RGG_4D,5000,8,isotropy,0.72,0.096,,
coarsegraining,RGG_3D,5000,8,retention_best,0.59,0.05,,
coarsegraining,RGG_4D,5000,8,retention_best,0.42,0.03,,
```
