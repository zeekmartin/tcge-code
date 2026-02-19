#!/usr/bin/env python3
"""
TCGE GAP-Dimensionality — Test de robustesse de l'invariant tri/diam

Question : Δ ≈ a × (tri/diam) + b tient-il aussi quand on varie ⟨k⟩ 
à dimension fixe ?

Si oui → c'est un invariant structurel, pas un artefact du scan en d.
Si non → la relation n'est qu'une corrélation accidentelle sur 6 points.

Protocole :
  - d fixe ∈ {3, 4, 5}
  - ⟨k⟩ ∈ {5, 8, 12, 16, 20}  
  - N = 3000, 6 trials par (d, k)
  - Mesurer Δ, <tri>, diam, ratio tri/diam

Author: Zeek (TCGE project) — AI-assisted implementation
Date: 2026-02-19
"""

import numpy as np
from collections import defaultdict, deque
from math import gamma as gamma_func
from scipy.spatial import cKDTree
import time

class RGG_Torus:
    def __init__(self, n_nodes, dim, target_degree=8, seed=None):
        if seed is not None:
            np.random.seed(seed)
        self.n = n_nodes
        self.dim = dim
        self.positions = np.random.random((n_nodes, dim))
        vol_d = np.pi**(dim/2) / gamma_func(dim/2 + 1)
        r = (target_degree / ((n_nodes - 1) * vol_d))**(1.0/dim)
        self.r_connect = r
        tree = cKDTree(self.positions, boxsize=np.ones(dim))
        pairs = tree.query_pairs(r)
        self.edges = [(min(i,j), max(i,j)) for i,j in pairs]
        self.adj = defaultdict(set)
        for i, j in self.edges:
            self.adj[i].add(j)
            self.adj[j].add(i)
        self.n_edges = len(self.edges)
        self.degree = np.array([len(self.adj[i]) for i in range(n_nodes)])
        self.N_prop = self.degree.astype(float)
        self.tri = np.zeros(self.n_edges, dtype=int)
        for idx, (i, j) in enumerate(self.edges):
            self.tri[idx] = len(self.adj[i] & self.adj[j])
        diam = 0
        for s in np.random.choice(n_nodes, size=min(5, n_nodes), replace=False):
            dist = {s: 0}
            queue = deque([s])
            while queue:
                v = queue.popleft()
                for u in self.adj[v]:
                    if u not in dist:
                        dist[u] = dist[v] + 1
                        queue.append(u)
            if dist:
                diam = max(diam, max(dist.values()))
        self.diameter_est = diam

    def biphasage_optimize(self, n_steps=2000, seed=None):
        if seed is not None:
            np.random.seed(seed)
        n_e = self.n_edges
        if n_e == 0:
            self.biphasage = 0
            return 0
        alpha_e = np.random.randn(n_e) * 0.01
        delta_N = np.array([self.N_prop[i] - self.N_prop[j] for i, j in self.edges])
        tri_f = self.tri.astype(float)
        for _ in range(n_steps):
            grad = (2 * alpha_e * (-0.25 + 0.5*tri_f + 0.01) + 0.5 * delta_N)
            alpha_e -= 0.005 * grad
            alpha_e = np.clip(alpha_e, -0.999, 0.999)
        abs_alpha = np.abs(alpha_e)
        med = np.median(self.tri)
        if med == self.tri.min():
            med = np.mean(self.tri)
        low = [i for i in range(n_e) if self.tri[i] <= med]
        high = [i for i in range(n_e) if self.tri[i] > med]
        a_l = np.mean(abs_alpha[low]) if low else 0
        a_h = np.mean(abs_alpha[high]) if high else 0
        self.biphasage = a_l - a_h
        return self.biphasage


def run():
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  Test de robustesse : Δ ~ tri/diam en variant ⟨k⟩         ║")
    print("║  d ∈ {3,4,5}, ⟨k⟩ ∈ {5,8,12,16,20}, N=3000              ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    dims = [3, 4, 5]
    ks = [5, 8, 12, 16, 20]
    N = 3000
    n_trials = 6
    
    all_points = []  # (tri/diam, Δ, d, k) for global fit
    
    for dim in dims:
        print(f"\n{'═'*65}")
        print(f"  d = {dim}")
        print(f"{'═'*65}")
        print(f"  {'⟨k⟩':<6} {'Δ':<12} {'<tri>':<8} {'diam':<6} "
              f"{'tri/diam':<10} {'<k>_real':<8}")
        print(f"  {'-'*55}")
        
        for k in ks:
            deltas = []
            tri_means = []
            diams = []
            k_reals = []
            
            for trial in range(n_trials):
                g = RGG_Torus(N, dim, target_degree=k, seed=1000*trial + dim*100 + k)
                bi = g.biphasage_optimize(seed=2000*trial + 7)
                deltas.append(bi)
                tri_means.append(np.mean(g.tri))
                diams.append(g.diameter_est)
                k_reals.append(g.degree.mean())
            
            md = np.mean(deltas)
            sd = np.std(deltas)
            mt = np.mean(tri_means)
            mdiam = np.mean(diams)
            ratio = mt / mdiam if mdiam > 0 else 0
            mk = np.mean(k_reals)
            
            all_points.append((ratio, md, dim, k, mt, mdiam))
            
            print(f"  {k:<6} {md:.3f}±{sd:.3f} {mt:<8.1f} {mdiam:<6.0f} "
                  f"{ratio:<10.3f} {mk:<8.1f}")
    
    # ── GLOBAL FIT ──
    print(f"\n\n{'═'*65}")
    print("  FIT GLOBAL : Δ = a × (tri/diam) + b")
    print(f"{'═'*65}\n")
    
    ratios = np.array([p[0] for p in all_points])
    deltas = np.array([p[1] for p in all_points])
    dims_arr = np.array([p[2] for p in all_points])
    ks_arr = np.array([p[3] for p in all_points])
    
    # Global correlation
    r_global = np.corrcoef(ratios, deltas)[0, 1]
    coeffs = np.polyfit(ratios, deltas, 1)
    predicted = np.polyval(coeffs, ratios)
    ss_res = np.sum((deltas - predicted)**2)
    ss_tot = np.sum((deltas - np.mean(deltas))**2)
    r2_global = 1 - ss_res / ss_tot
    
    print(f"  Points: {len(all_points)} ({len(dims)}×{len(ks)} conditions)")
    print(f"  Pearson r = {r_global:.3f}")
    print(f"  R² = {r2_global:.3f}")
    print(f"  Fit: Δ = {coeffs[0]:.3f} × (tri/diam) + {coeffs[1]:.3f}")
    print(f"  Residual std = {np.std(deltas - predicted):.4f}")
    
    # Per-dimension correlation
    print(f"\n  Per-dimension breakdown:")
    for dim in dims:
        mask = dims_arr == dim
        r_d = np.corrcoef(ratios[mask], deltas[mask])[0, 1]
        print(f"    d={dim}: r = {r_d:.3f} (n={sum(mask)})")
    
    # ── VERDICT ──
    print(f"\n{'═'*65}")
    print("  VERDICT")
    print(f"{'═'*65}\n")
    
    if r_global > 0.9:
        print(f"  ✅ INVARIANT CONFIRMÉ (r={r_global:.3f}, {len(all_points)} points)")
        print(f"     Δ ~ tri/diam tient en variant BOTH d et ⟨k⟩.")
        print(f"     Ce n'est pas un artefact du scan dimensionnel.")
        print()
        print(f"  Formulation :")
        print(f"  'The biphasage amplitude is governed by a single structural")
        print(f"   ratio — mean triangle count per unit graph diameter —")
        print(f"   across both embedding dimensions (d=3,4,5) and mean")
        print(f"   degrees (k=5–20), with R²={r2_global:.2f} on {len(all_points)} conditions.'")
    elif r_global > 0.7:
        print(f"  ⚠️ RELATION PARTIELLE (r={r_global:.3f})")
        print(f"     tri/diam capture la tendance mais pas toute la variance.")
    else:
        print(f"  ❌ PAS D'INVARIANT (r={r_global:.3f})")
        print(f"     La relation ne tient pas quand on varie ⟨k⟩.")
    
    return all_points


if __name__ == "__main__":
    t0 = time.time()
    pts = run()
    print(f"\n  Temps total: {time.time()-t0:.0f}s")
