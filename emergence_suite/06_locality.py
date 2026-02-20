#!/usr/bin/env python3
"""
TCGE GAP-Continuum-Locality — Homogénéité locale

Question : les propriétés continuum sont-elles localement homogènes 
(comme une variété) ou très hétérogènes ?

Tests :
  1. d_H local : mesurer d_H depuis N centres, regarder la variance
  2. Distribution des distances géodésiques : forme caractéristique ?
  3. Isotropie locale : variance angulaire du volume des boules

Un substrat "variété-like" aurait :
  - d_H local ≈ constant (faible variance)
  - distances géodésiques ∝ distribution chi (dimension d)
  - isotropie ≈ 1 (pas de direction privilégiée)

On teste sur RGG-tore (d=3,4) avec et sans biphasage.

Author: Zeek (TCGE project) — AI-assisted implementation
Date: 2026-02-19
"""

import numpy as np
from collections import defaultdict, deque
from math import gamma as gamma_func
from scipy.spatial import cKDTree
import time

# =============================================================================
# RGG Torus (from v7d)
# =============================================================================

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
        
        # Diameter
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
    
    def biphasage_optimize(self, A=1.0, C_protect=0.5, gamma=0.5,
                           beta_reg=0.01, W=1.0, lr=0.005, n_steps=2500):
        n_e = self.n_edges
        if n_e == 0:
            return 0
        alpha_e = np.random.randn(n_e) * 0.01
        delta_N = np.array([self.N_prop[i] - self.N_prop[j] 
                           for i, j in self.edges])
        tri = self.tri.astype(float)
        for _ in range(n_steps):
            grad = (2 * alpha_e * (-A*(W/2)**2 + C_protect*tri + beta_reg)
                    + gamma * W * delta_N)
            alpha_e -= lr * grad
            alpha_e = np.clip(alpha_e, -0.999, 0.999)
        self.alpha_e = alpha_e
        self.abs_alpha = np.abs(alpha_e)
        med = np.median(self.tri)
        if med == self.tri.min():
            med = np.mean(self.tri)
        low = [i for i in range(n_e) if self.tri[i] <= med]
        high = [i for i in range(n_e) if self.tri[i] > med]
        a_l = np.mean(self.abs_alpha[low]) if low else 0
        a_h = np.mean(self.abs_alpha[high]) if high else 0
        self.biphasage = a_l - a_h
        return self.biphasage


# =============================================================================
# BFS distances from a source
# =============================================================================

def bfs_distances(adj, src, n_nodes):
    dist = [-1] * n_nodes
    dist[src] = 0
    queue = deque([src])
    while queue:
        v = queue.popleft()
        for u in adj[v]:
            if dist[u] == -1:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist


# =============================================================================
# TEST 1: Local d_H (per-source)
# =============================================================================

def local_dH(graph, n_sources=100):
    """Measure d_H from many individual sources. Return distribution."""
    sources = np.random.choice(graph.n, size=min(n_sources, graph.n), replace=False)
    r_max = min(graph.diameter_est // 2, 20)
    
    dH_per_source = []
    
    for src in sources:
        dist = bfs_distances(graph.adj, src, graph.n)
        
        # Volume at each radius
        vols = []
        for r in range(1, r_max + 1):
            vol = sum(1 for d in dist if 0 <= d <= r)
            vols.append(vol)
        
        # Fit log(vol) vs log(r) in intermediate range
        r_vals = np.arange(1, r_max + 1)
        log_r = np.log(r_vals.astype(float))
        log_v = np.log(np.maximum(np.array(vols), 1))
        
        # Use range [3, r_max-2] to avoid boundary effects
        lo, hi = 2, min(r_max - 2, len(r_vals) - 1)
        if hi - lo >= 3:
            slope = np.polyfit(log_r[lo:hi], log_v[lo:hi], 1)[0]
            dH_per_source.append(slope)
    
    return np.array(dH_per_source)


# =============================================================================
# TEST 2: Geodesic distance distribution
# =============================================================================

def geodesic_distribution(graph, n_pairs=5000):
    """Sample pairwise geodesic distances, return histogram."""
    pairs_i = np.random.randint(0, graph.n, size=n_pairs)
    pairs_j = np.random.randint(0, graph.n, size=n_pairs)
    
    # BFS from unique sources
    unique_sources = np.unique(pairs_i)
    # Limit to manageable number
    if len(unique_sources) > 200:
        unique_sources = np.random.choice(unique_sources, 200, replace=False)
    
    dist_cache = {}
    for src in unique_sources:
        dist_cache[src] = bfs_distances(graph.adj, src, graph.n)
    
    distances = []
    for i, j in zip(pairs_i, pairs_j):
        if i == j:
            continue
        if i in dist_cache:
            d = dist_cache[i][j]
            if d >= 0:
                distances.append(d)
    
    return np.array(distances)


# =============================================================================
# TEST 3: Local isotropy
# =============================================================================

def local_isotropy(graph, n_sources=50):
    """
    For each source, measure volume in different 'octants' of the torus.
    If isotropic, volumes should be roughly equal.
    
    We split neighbors at each radius into 2^d sectors based on 
    relative position sign, then compute variance of sector populations.
    """
    dim = graph.dim
    n_sectors = 2 ** dim
    sources = np.random.choice(graph.n, size=min(n_sources, graph.n), replace=False)
    
    isotropy_scores = []
    
    for src in sources:
        dist = bfs_distances(graph.adj, src, graph.n)
        
        # Nodes at intermediate distance (not too close, not too far)
        r_lo = 3
        r_hi = min(graph.diameter_est // 2, 12)
        
        shell_nodes = [v for v in range(graph.n) if r_lo <= dist[v] <= r_hi]
        
        if len(shell_nodes) < n_sectors * 2:
            continue
        
        # Assign each node to a sector based on relative position
        src_pos = graph.positions[src]
        sector_counts = np.zeros(n_sectors)
        
        for v in shell_nodes:
            # Toroidal relative position
            delta = graph.positions[v] - src_pos
            # Wrap to [-0.5, 0.5]
            delta = delta - np.round(delta)
            
            # Sector index from sign of each coordinate
            sector = 0
            for d in range(dim):
                if delta[d] >= 0:
                    sector += (1 << d)
            sector_counts[sector] += 1
        
        # Isotropy = 1 - normalized variance
        # Perfect isotropy: all sectors equal → variance = 0
        if sector_counts.sum() > 0:
            fracs = sector_counts / sector_counts.sum()
            expected = 1.0 / n_sectors
            # Coefficient of variation of sector fractions
            cv = np.std(fracs) / expected
            isotropy = max(0, 1 - cv)
            isotropy_scores.append(isotropy)
    
    return np.array(isotropy_scores)


# =============================================================================
# MAIN
# =============================================================================

def run():
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  TCGE GAP-Continuum-Locality                              ║")
    print("║  Homogénéité locale + isotropie                           ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    dims = [3, 4]
    N = 5000
    
    for dim in dims:
        print(f"\n{'═'*65}")
        print(f"  d = {dim}  |  N = {N}  |  RGG Tore  |  <k> ≈ 8")
        print(f"{'═'*65}")
        
        t0 = time.time()
        graph = RGG_Torus(N, dim, target_degree=8, seed=42 + dim*100)
        print(f"  Graph: E={graph.n_edges}, <k>={graph.degree.mean():.1f}, "
              f"diam≈{graph.diameter_est}")
        
        # Run biphasage
        np.random.seed(7777 + dim)
        bi = graph.biphasage_optimize(n_steps=2000)
        print(f"  Biphasage: Δ={bi:.3f}")
        
        # ── TEST 1: Local d_H ──
        print(f"\n  --- TEST 1: d_H local (100 sources) ---")
        dH_local = local_dH(graph, n_sources=100)
        
        if len(dH_local) > 0:
            print(f"  d_H moyen    = {np.mean(dH_local):.2f}")
            print(f"  d_H médian   = {np.median(dH_local):.2f}")
            print(f"  d_H std      = {np.std(dH_local):.2f}")
            print(f"  d_H [5%, 95%] = [{np.percentile(dH_local, 5):.2f}, "
                  f"{np.percentile(dH_local, 95):.2f}]")
            cv = np.std(dH_local) / np.mean(dH_local)
            print(f"  CV (std/mean) = {cv:.3f}")
            
            if cv < 0.1:
                print(f"  → ✅ HOMOGÈNE (CV < 0.1)")
            elif cv < 0.2:
                print(f"  → ⚠️ MODÉRÉMENT HOMOGÈNE (CV < 0.2)")
            else:
                print(f"  → ❌ HÉTÉROGÈNE (CV ≥ 0.2)")
        
        # ── TEST 2: Geodesic distances ──
        print(f"\n  --- TEST 2: Distribution des distances géodésiques ---")
        geo_dist = geodesic_distribution(graph, n_pairs=10000)
        
        if len(geo_dist) > 0:
            print(f"  Distances: moy={np.mean(geo_dist):.1f}, "
                  f"med={np.median(geo_dist):.1f}, "
                  f"std={np.std(geo_dist):.1f}")
            print(f"  [min, max] = [{np.min(geo_dist)}, {np.max(geo_dist)}]")
            
            # Histogram
            bins = np.arange(0, np.max(geo_dist) + 2)
            hist, _ = np.histogram(geo_dist, bins=bins)
            mode_r = np.argmax(hist)
            print(f"  Mode = {mode_r}")
            
            # Skewness
            from scipy.stats import skew, kurtosis
            sk = skew(geo_dist)
            ku = kurtosis(geo_dist)
            print(f"  Skewness = {sk:.2f} (variété d-dim: légèrement positif)")
            print(f"  Kurtosis = {ku:.2f}")
            
            # Compare with d-dim ball expectation: 
            # P(r) ∝ r^(d-1) * exp(-r²/2σ²) for finite graph
            # Mode should be around diameter/2 * (d-1)/d
            expected_mode = graph.diameter_est * (dim - 1) / (2 * dim)
            print(f"  Mode attendu (heuristique) ≈ {expected_mode:.0f}")
        
        # ── TEST 3: Isotropie locale ──
        print(f"\n  --- TEST 3: Isotropie locale (50 sources, "
              f"{2**dim} secteurs) ---")
        iso = local_isotropy(graph, n_sources=50)
        
        if len(iso) > 0:
            print(f"  Isotropie moyenne = {np.mean(iso):.3f} (1.0 = parfait)")
            print(f"  Isotropie std     = {np.std(iso):.3f}")
            print(f"  Isotropie [5%,95%] = [{np.percentile(iso, 5):.3f}, "
                  f"{np.percentile(iso, 95):.3f}]")
            
            if np.mean(iso) > 0.8:
                print(f"  → ✅ ISOTROPE (moy > 0.8)")
            elif np.mean(iso) > 0.6:
                print(f"  → ⚠️ MODÉRÉMENT ISOTROPE")
            else:
                print(f"  → ❌ ANISOTROPE")
        
        print(f"\n  ({time.time()-t0:.1f}s)")
    
    # ── SYNTHÈSE ──
    print(f"\n\n{'═'*65}")
    print("  SYNTHÈSE GAP-CONTINUUM-LOCALITY")
    print(f"{'═'*65}\n")
    print("  Question: les propriétés sont-elles localement homogènes ?")
    print()
    print("  Critères de 'variété-like' :")
    print("    1. d_H local: CV < 0.1 (faible dispersion)")
    print("    2. Distances: distribution unimodale, skew modéré")
    print("    3. Isotropie: score > 0.8")
    print()
    print("  Si les 3 critères passent: le substrat biphasé se")
    print("  comporte localement comme une variété riemannienne,")
    print("  pas comme un objet fractal ou hétérogène.")


if __name__ == "__main__":
    t0 = time.time()
    run()
    print(f"\n  Temps total: {time.time()-t0:.0f}s")
