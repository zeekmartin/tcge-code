# TCGE — Quantitative Results

## 1. Biphasage T/S (v5)

**Setup:** Erdős–Rényi graph, N=60, p=0.25, 30 trials.  
**Mechanism:** E = Σ_e [A·w_fwd·w_bwd + C·tri(e)·α² + γ·Δw·ΔN]  
**Optimal parameters:** A=1.0, C_protect=0.5, γ=0.5

| Metric | Value |
|--------|-------|
| |α|_lowTri (proto-temporal) | 0.697 ± 0.057 |
| |α|_highTri (proto-spatial) | 0.382 ± 0.051 |
| **Biphasage Δ** | **0.315 ± 0.055** |
| Directional coherence (low) | 91.8% |
| Directional coherence (high) | 91.0% |
| Trials with Δ > 0.2 | 97% |
| Trials with Δ > 0.1 | 100% |

### Detailed |α| vs tri(e) (median trial)

| tri(e) | N edges | |α| mean | Class |
|--------|---------|---------|-------|
| 0 | 8 | 0.999 | TEMPORAL |
| 1 | 34 | 0.963 | TEMPORAL |
| 2 | 75 | 0.737 | TEMPORAL |
| 3 | 75 | 0.509 | TEMPORAL |
| 4 | 108 | 0.397 | TEMPORAL |
| 5 | 83 | 0.333 | SPATIAL |
| 6 | 47 | 0.237 | SPATIAL |
| 7 | 26 | 0.207 | SPATIAL |
| 8 | 10 | 0.226 | SPATIAL |
| 9 | 6 | 0.117 | SPATIAL |

### Parameter scan (C_protect)

| C_protect | |α|_low | |α|_high | Δ | Status |
|-----------|--------|---------|---|--------|
| 0.0 | 0.995 | 0.997 | -0.002 | ❌ |
| 0.1 | 0.926 | 0.875 | 0.051 | ❌ |
| 0.2 | 0.889 | 0.707 | 0.182 | ⚠️ |
| **0.3** | **0.822** | **0.563** | **0.259** | **✅** |
| **0.5** | **0.694** | **0.373** | **0.321** | **✅** |
| **0.8** | **0.533** | **0.232** | **0.301** | **✅** |
| **1.0** | **0.458** | **0.184** | **0.274** | **✅** |
| 1.5 | 0.336 | 0.120 | 0.216 | ✅ |
| 2.0 | 0.266 | 0.090 | 0.176 | ⚠️ |
| 3.0 | 0.190 | 0.059 | 0.131 | ⚠️ |


## 2. Robustness (v5b)

### Test A — Multiple cohesion metrics (N=60, 25 trials)

| Metric | Δ | Status |
|--------|---|--------|
| triangles | 0.321 | ✅ |
| jaccard | 0.318 | ✅ |
| edge_clustering | 0.204 | ✅ |
| quadrangles | 0.213 | ✅ |
| truss | 0.208 | ✅ |
| compat_destroy | 0.139 | ⚠️ |

### Test B — TCGE-native reformulation

| Measure | Δ |
|---------|---|
| Triangles | 0.321 ± 0.057 |
| Compatibility destruction | 0.139 ± 0.069 |
| Trial-by-trial correlation | r = 0.822 |

### Test C — Finite-size scaling

| N | Δ | % > 0.2 |
|---|---|---------|
| 30 | 0.150 | — |
| 50 | 0.304 | — |
| 70 | 0.313 | 95% |
| 100 | 0.324 | 85% |
| 150 | 0.259 | — |
| 200 | 0.232 | — |
| 300 | 0.225 | — |

Log-log slope: **+0.063 ≈ 0** (plateau, not finite-size artifact)


## 3. Arrow (v4)

**Setup:** Erdős–Rényi graph, N=50, p=0.25, 30 trials per config.

| Config | |α| mean | α_T | α_S | Aniso(T-S) | Coherence |
|--------|---------|-----|-----|------------|-----------|
| Product alone | 0.714 | 0.715 | 0.714 | 0.024 | 0.50 |
| Product+thermo (γ=0.5) | 0.973 | 0.999 | 0.954 | 0.045 | 1.00 |
| Strong product (A=5) | 0.999 | 0.999 | 0.999 | 0.000 | 0.50 |
| Strong product+thermo | 0.999 | 0.999 | 0.999 | 0.000 | 1.00 |


## 4. Continuum diagnostics (v7d)

**Setup:** RGG on torus, N=5000, ⟨k⟩≈8, 5 trials per dimension.

### Hausdorff dimension

| d (embedding) | d_H (plateau) | d_H (global) | d_H std |
|:---:|:---:|:---:|:---:|
| 3 | **2.84** | 2.63 | 0.03 |
| 4 | **3.38** | 3.24 | 0.02 |

### Spectral dimension

| d | d_s (plateau) | d_s (global) | d_s std |
|:---:|:---:|:---:|:---:|
| 3 | 2.37 | 2.53 | 0.61 |
| 4 | 2.71 | 3.17 | 0.52 |

### Biphasage on geometric substrates

| d | Δ | |α|_low | |α|_high |
|:---:|:---:|:---:|:---:|
| 3 | **0.376** | 0.601 | 0.220 |
| 4 | **0.400** | 0.677 | 0.274 |

### Coarse-graining retention (k~100, density-based cohesion)

| d | Best retention (density) | Best retention (Jaccard) |
|:---:|:---:|:---:|
| 3 | **60%** (trial 5) | 18% |
| 4 | **39%** | 15% |

Mean retention (density, k~100): d=3: 44%, d=4: 34%

### Improvement from torus (vs v7c boundary)

| d | d_H (boundary) | d_H (torus) | Δd_H |
|:---:|:---:|:---:|:---:|
| 3 | 2.50 | 2.84 | +0.34 |
| 4 | 3.01 | 3.38 | +0.37 |


## 5. GAP status summary

| GAP | Status | Key evidence |
|-----|--------|-------------|
| Signature (T/S) | ✅ Closed | v5: Δ=0.31, v5b: 5/6 metrics, plateau scaling |
| Arrow (direction) | ✅ Closed | v4: 91% coherence, thermo+degree |
| Continuum (d_H) | 🟢 Quasi-closed | v7d: stable plateau, torus correction |
| Continuum (d_s) | 🟡 Partial | v7d: signal but σ≈0.5 |
| Continuum (CG) | 🟡 Substantial | v7d: up to 60%, mean 44% |
