#!/usr/bin/env python3
"""
TCGE — Contrôle causal : rewiring à degré fixe

Question : Δ suit-il ⟨tri⟩ quand on casse les triangles SANS changer ⟨k⟩ ?

Protocole :
  1. Générer un RGG-tore (d=3, N=3000, ⟨k⟩=8)
  2. Appliquer des niveaux croissants de rewiring (0%, 10%, 20%, ..., 90%)
     - Le rewiring échange les extrémités de deux arêtes aléatoires
     - Préserve exactement la séquence de degrés
     - Détruit progressivement les triangles
  3. Mesurer (⟨tri⟩, Δ) à chaque niveau
  4. Vérifier que la relation Δ ~ -⟨tri⟩ tient

Si oui → le driver est CAUSALEMENT le clustering, pas le degré.

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
        tree = cKDTree(self.positions, boxsize=np.ones(dim))
        pairs = tree.query_pairs(r)
        self.edges = [(min(i,j), max(i,j)) for i,j in pairs]
        self.adj = defaultdict(set)
        for i, j in self.edges:
            self.adj[i].add(j)
            self.adj[j].add(i)
        self.n_edges = len(self.edges)
        self.degree = np.array([len(self.adj[i]) for i in range(n_nodes)])


def degree_preserving_rewire(edges, adj, n_nodes, n_swaps, seed=None):
    """
    Degree-preserving rewiring (Maslov-Sneppen).
    Pick two random edges (a-b, c-d), swap to (a-d, c-b).
    Reject if self-loop or multi-edge.
    Preserves degree sequence exactly.
    """
    if seed is not None:
        np.random.seed(seed)
    
    edges = list(edges)
    edge_set = set(edges)
    adj = defaultdict(set)
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    
    n_e = len(edges)
    swaps_done = 0
    attempts = 0
    max_attempts = n_swaps * 10
    
    while swaps_done < n_swaps and attempts < max_attempts:
        attempts += 1
        
        # Pick two random edges
        idx1 = np.random.randint(n_e)
        idx2 = np.random.randint(n_e)
        if idx1 == idx2:
            continue
        
        a, b = edges[idx1]
        c, d = edges[idx2]
        
        # Choose swap direction randomly
        if np.random.random() < 0.5:
            new1 = (min(a, d), max(a, d))
            new2 = (min(c, b), max(c, b))
        else:
            new1 = (min(a, c), max(a, c))
            new2 = (min(b, d), max(b, d))
        
        # Reject self-loops
        if new1[0] == new1[1] or new2[0] == new2[1]:
            continue
        
        # Reject multi-edges
        if new1 in edge_set or new2 in edge_set:
            continue
        
        # Reject if new1 == new2
        if new1 == new2:
            continue
        
        # Accept swap
        edge_set.discard(edges[idx1])
        edge_set.discard(edges[idx2])
        
        adj[a].discard(b); adj[b].discard(a)
        adj[c].discard(d); adj[d].discard(c)
        
        edges[idx1] = new1
        edges[idx2] = new2
        edge_set.add(new1)
        edge_set.add(new2)
        
        adj[new1[0]].add(new1[1]); adj[new1[1]].add(new1[0])
        adj[new2[0]].add(new2[1]); adj[new2[1]].add(new2[0])
        
        swaps_done += 1
    
    return edges, adj


def compute_biphasage(n_nodes, edges, adj, seed=None):
    if seed is not None:
        np.random.seed(seed)
    
    n_e = len(edges)
    if n_e < 10:
        return 0, 0, 0
    
    degree = np.zeros(n_nodes)
    for i, j in edges:
        degree[i] += 1
        degree[j] += 1
    
    tri = np.zeros(n_e, dtype=int)
    for idx, (i, j) in enumerate(edges):
        tri[idx] = len(adj[i] & adj[j])
    
    mean_tri = np.mean(tri)
    
    alpha_e = np.random.randn(n_e) * 0.01
    delta_N = np.array([degree[i] - degree[j] for i, j in edges])
    tri_f = tri.astype(float)
    
    for _ in range(2500):
        grad = (2 * alpha_e * (-0.25 + 0.5*tri_f + 0.01) + 0.5 * delta_N)
        alpha_e -= 0.005 * grad
        alpha_e = np.clip(alpha_e, -0.999, 0.999)
    
    abs_alpha = np.abs(alpha_e)
    med = np.median(tri)
    if med == tri.min():
        med = np.mean(tri)
    low = [i for i in range(n_e) if tri[i] <= med]
    high = [i for i in range(n_e) if tri[i] > med]
    
    if not low or not high:
        return 0, mean_tri, np.mean(degree)
    
    biphasage = np.mean(abs_alpha[low]) - np.mean(abs_alpha[high])
    return biphasage, mean_tri, np.mean(degree)


def run():
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  Contrôle causal : rewiring à degré fixe                  ║")
    print("║  RGG d=3, N=3000, ⟨k⟩=8                                  ║")
    print("║  Rewire 0% → 90% des arêtes                              ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    N = 3000
    dim = 3
    n_trials = 6
    rewire_fracs = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
    
    all_points = []
    
    print(f"\n  {'rewire%':<10} {'Δ':<12} {'⟨tri⟩':<10} {'⟨k⟩':<8}")
    print(f"  {'-'*45}")
    
    for frac in rewire_fracs:
        deltas = []
        tris = []
        ks = []
        
        for trial in range(n_trials):
            # Generate base graph
            graph = RGG_Torus(N, dim, target_degree=8, seed=1000*trial + 42)
            edges = list(graph.edges)
            adj = defaultdict(set)
            for i, j in edges:
                adj[i].add(j)
                adj[j].add(i)
            
            # Rewire
            n_swaps = int(frac * len(edges))
            if n_swaps > 0:
                edges, adj = degree_preserving_rewire(
                    edges, adj, N, n_swaps, seed=3000*trial + int(frac*100))
            
            # Measure
            bi, tri_mean, k_mean = compute_biphasage(N, edges, adj, seed=2000*trial + 7)
            deltas.append(bi)
            tris.append(tri_mean)
            ks.append(k_mean)
        
        md = np.mean(deltas)
        sd = np.std(deltas)
        mt = np.mean(tris)
        mk = np.mean(ks)
        
        all_points.append((frac, md, sd, mt, mk))
        
        print(f"  {frac*100:>5.0f}%    {md:.3f}±{sd:.3f} {mt:<10.2f} {mk:<8.1f}")
    
    # ── FIT ──
    print(f"\n\n{'═'*55}")
    print("  FIT : Δ vs ⟨tri⟩ (rewiring control)")
    print(f"{'═'*55}\n")
    
    tris_arr = np.array([p[3] for p in all_points])
    deltas_arr = np.array([p[1] for p in all_points])
    ks_arr = np.array([p[4] for p in all_points])
    
    r = np.corrcoef(tris_arr, deltas_arr)[0, 1]
    coeffs = np.polyfit(tris_arr, deltas_arr, 1)
    pred = np.polyval(coeffs, tris_arr)
    r2 = 1 - np.sum((deltas_arr - pred)**2) / np.sum((deltas_arr - np.mean(deltas_arr))**2)
    
    print(f"  Δ = {coeffs[0]:.4f} × ⟨tri⟩ + {coeffs[1]:.3f}")
    print(f"  r = {r:.3f}")
    print(f"  R² = {r2:.3f}")
    print(f"  n = {len(all_points)} conditions")
    
    # Degree control
    k_range = ks_arr.max() - ks_arr.min()
    print(f"\n  Contrôle degré : ⟨k⟩ varie de {ks_arr.min():.1f} à {ks_arr.max():.1f} "
          f"(range = {k_range:.2f})")
    print(f"  → degré {'CONSTANT' if k_range < 0.5 else 'variable'} "
          f"({'✅ contrôle valide' if k_range < 0.5 else '⚠️ léger drift'})")
    
    # ── COMPARISON WITH K-SCAN LAW ──
    print(f"\n{'═'*55}")
    print("  COMPARAISON avec la loi du k-scan")
    print(f"{'═'*55}\n")
    print(f"  k-scan (d=3-5, k=5-20):  Δ = -0.034 × ⟨tri⟩ + 0.507")
    print(f"  Rewire (d=3, k=8 fixe):  Δ = {coeffs[0]:.4f} × ⟨tri⟩ + {coeffs[1]:.3f}")
    print(f"  Pentes : {-0.034:.4f} vs {coeffs[0]:.4f}")
    
    if abs(coeffs[0] - (-0.034)) < 0.015:
        print(f"  → ✅ PENTES CONCORDANTES (Δ < 0.015)")
    else:
        print(f"  → ⚠️ Pentes différentes (écart = {abs(coeffs[0]-(-0.034)):.4f})")
    
    # ── VERDICT ──
    print(f"\n{'═'*55}")
    print("  VERDICT")
    print(f"{'═'*55}\n")
    
    if abs(r) > 0.95 and k_range < 0.5:
        print(f"  ✅ CAUSALITÉ CONFIRMÉE")
        print(f"     Δ suit ⟨tri⟩ quand on détruit les triangles")
        print(f"     par rewiring, à séquence de degrés EXACTE.")
        print(f"     Le degré n'est pas un confondeur.")
        print()
        print(f"  Formulation :")
        print(f"  'Degree-preserving rewiring confirms the causal role of")
        print(f"   triangle abundance: destroying triangles while holding")
        print(f"   the degree sequence exactly fixed increases biphasage")
        print(f"   (r = {r:.2f}, n = {len(all_points)}), ruling out degree")
        print(f"   as a confounding variable.'")
    elif abs(r) > 0.8:
        print(f"  ⚠️ RELATION FORTE mais pas parfaite (r={r:.3f})")
    else:
        print(f"  ❌ RELATION FAIBLE (r={r:.3f})")
    
    return all_points


if __name__ == "__main__":
    t0 = time.time()
    pts = run()
    print(f"\n  Temps total: {time.time()-t0:.0f}s")
