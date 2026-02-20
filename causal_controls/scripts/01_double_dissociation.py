#!/usr/bin/env python3
"""
TCGE Double Dissociation — 4-condition causal control
=====================================================

Conditions (2×2 factorial):
  ER  : natural graph + protector ON   (triangles ✓, protect ✓)
  C1  : rewired graph + protector ON   (triangles ✗, protect ✓)
  C2  : natural graph + protector OFF  (triangles ✓, protect ✗)
  REF : rewired graph + protector OFF  (triangles ✗, protect ✗)

Metrics:
  S   = biphasage separation (⟨|α|⟩_lowTri − ⟨|α|⟩_highTri)
  μ₁  = spectral gap proxy (variance of |α| distribution)
  Also report: ⟨tri⟩, ⟨|α|⟩

N=60, p=0.25, 25 seeds. Should run in <90s.
"""

import numpy as np
import networkx as nx
import time
from collections import defaultdict

# ──────────────────────────────────────────────────────────────
# GRAPH GENERATION
# ──────────────────────────────────────────────────────────────

def make_graph(N, k=8, p_ws=0.1, seed=42):
    """Watts-Strogatz graph: high clustering that rewiring can destroy."""
    G = nx.watts_strogatz_graph(N, k, p_ws, seed=seed)
    return G


def degree_preserving_rewire(G, fraction=0.9, seed=42):
    """Destroy triangles via degree-preserving rewiring using NetworkX."""
    rng = np.random.RandomState(seed)
    G2 = G.copy()
    M = G2.number_of_edges()
    n_swaps = int(M * fraction * 5)
    
    try:
        nx.double_edge_swap(G2, nswap=n_swaps, max_tries=n_swaps * 10, seed=seed)
    except nx.NetworkXAlgorithmError:
        pass  # Hit max tries, that's OK
    
    return G2


# ──────────────────────────────────────────────────────────────
# TRIANGLE COUNT PER EDGE
# ──────────────────────────────────────────────────────────────

def edge_triangles(G):
    """Count triangles per edge."""
    tri = {}
    for u, v in G.edges():
        common = len(set(G.neighbors(u)) & set(G.neighbors(v)))
        tri[(u, v)] = common
        tri[(v, u)] = common
    return tri


# ──────────────────────────────────────────────────────────────
# OPTIMIZATION
# ──────────────────────────────────────────────────────────────

def optimize(G, C_protect=0.5, A=1.0, beta_reg=0.01, n_iter=3000, seed=42):
    """
    Optimize edge weights w ∈ [-W, W] under:
      cost = Σ_edges [ -A*w² + β*w²  +  C_protect * tri(e) * (w/W)² ]
    
    -A*w² drives |w| → W (polarization, symmetry breaking)
    C_protect*tri(e)*(w/W)² resists polarization on high-tri edges
    β*w² is a small regularization
    
    Net effect per edge: if A > C_protect*tri(e) + β, edge polarizes.
    Uses gradient ascent on -cost (= gradient descent on cost) + noise.
    """
    rng = np.random.RandomState(seed)
    edges = list(G.edges())
    M = len(edges)
    W = 1.0
    
    # Get triangle counts
    tri_dict = edge_triangles(G)
    tri_arr = np.array([tri_dict.get((u, v), 0) for u, v in edges], dtype=float)
    
    # Initialize with small random perturbation (breaks symmetry)
    w = rng.uniform(-0.15, 0.15, M)
    
    lr = 0.02
    noise_scale = 0.02
    
    for it in range(n_iter):
        # Effective coefficient per edge: (-A + β + C_protect*tri/W²)
        # If negative → edge wants to polarize
        # If positive → edge wants to stay at zero
        eff = -A + beta_reg + C_protect * tri_arr / (W**2)
        
        # Gradient of cost w.r.t. w: 2 * eff * w
        grad = 2 * eff * w
        
        # Gradient descent on cost
        w -= lr * grad
        
        # Add noise (annealing)
        if it < n_iter * 0.7:
            noise = rng.normal(0, noise_scale * (1 - it / n_iter), M)
            w += noise
        
        # Clip to [-W, W]
        w = np.clip(w, -W, W)
        
        # Decrease lr at midpoint
        if it == n_iter // 2:
            lr *= 0.3
    
    # Compute alpha = |w| / W (polarization)
    alpha = np.abs(w) / W
    
    return w, alpha, tri_arr, edges


# ──────────────────────────────────────────────────────────────
# METRICS
# ──────────────────────────────────────────────────────────────

def compute_metrics(alpha, tri_arr):
    """
    S  = separation = ⟨|α|⟩_lowTri − ⟨|α|⟩_highTri   (>0 means biphasage)
    μ₁ = std(|α|)  (proxy for spectral gap / bimodality)
    """
    median_tri = np.median(tri_arr)
    
    # Handle edge case where all tri are equal
    if median_tri == 0:
        # Use tri > 0 vs tri == 0
        low_mask = tri_arr == 0
        high_mask = tri_arr > 0
    else:
        low_mask = tri_arr <= median_tri
        high_mask = tri_arr > median_tri
    
    if np.sum(low_mask) == 0 or np.sum(high_mask) == 0:
        S = 0.0
    else:
        S = np.mean(alpha[low_mask]) - np.mean(alpha[high_mask])
    
    mu1 = np.std(alpha)
    mean_alpha = np.mean(alpha)
    mean_tri = np.mean(tri_arr)
    
    return S, mu1, mean_alpha, mean_tri


# ──────────────────────────────────────────────────────────────
# RUN ONE CONDITION
# ──────────────────────────────────────────────────────────────

def run_condition(condition, N=80, n_seeds=25):
    """Run one condition across multiple seeds.
    
    For tri✓ conditions: Watts-Strogatz (high clustering)
    For tri✗ conditions: Erdős-Rényi at same ⟨k⟩ (random baseline triangles)
    """
    results = []
    k_ws = 8  # WS parameter
    p_er = k_ws / (N - 1)  # ER parameter to match mean degree
    
    for seed in range(n_seeds):
        if condition in ('C1', 'REF'):
            # Random graph — triangles at random baseline
            G = nx.erdos_renyi_graph(N, p_er, seed=1000 + seed)
            # Ensure connected
            if not nx.is_connected(G):
                comps = list(nx.connected_components(G))
                for c in comps[1:]:
                    u = list(c)[0]
                    v = list(comps[0])[0]
                    G.add_edge(u, v)
        else:
            # Structured graph — high clustering
            G = make_graph(N, k=k_ws, p_ws=0.1, seed=1000 + seed)
        
        C_prot = 0.0 if condition in ('C2', 'REF') else 0.5
        
        w, alpha, tri_arr, edges = optimize(G, C_protect=C_prot, seed=3000 + seed)
        S, mu1, mean_alpha, mean_tri = compute_metrics(alpha, tri_arr)
        
        # Extra metrics
        frac_spatial = np.mean(alpha < 0.5)  # fraction proto-spatial
        
        results.append({
            'S': S, 'mu1': mu1, 
            'mean_alpha': mean_alpha, 'mean_tri': mean_tri,
            'frac_spatial': frac_spatial
        })
    
    return results


# ──────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    t0 = time.time()
    N = 80
    n_seeds = 25
    
    conditions = ['ER', 'C1', 'C2', 'REF']
    labels = {
        'ER':  'WS  (tri↑, prot ON)',
        'C1':  'C1  (tri↓, prot ON)',
        'C2':  'C2  (tri↑, prot OFF)',
        'REF': 'REF (tri↓, prot OFF)'
    }
    
    all_results = {}
    
    for cond in conditions:
        print(f"Running {cond}... ", end='', flush=True)
        t1 = time.time()
        all_results[cond] = run_condition(cond, N=N, n_seeds=n_seeds)
        dt = time.time() - t1
        
        S_vals = [r['S'] for r in all_results[cond]]
        mu1_vals = [r['mu1'] for r in all_results[cond]]
        tri_vals = [r['mean_tri'] for r in all_results[cond]]
        alpha_vals = [r['mean_alpha'] for r in all_results[cond]]
        fs_vals = [r['frac_spatial'] for r in all_results[cond]]
        
        print(f"{dt:.1f}s  |  S={np.mean(S_vals):.3f}±{np.std(S_vals):.3f}  "
              f"μ₁={np.mean(mu1_vals):.3f}  "
              f"⟨tri⟩={np.mean(tri_vals):.2f}  ⟨|α|⟩={np.mean(alpha_vals):.3f}  "
              f"f_S={np.mean(fs_vals):.2f}")
    
    elapsed = time.time() - t0
    print(f"\nTotal: {elapsed:.0f}s")
    
    # ── Summary table ──
    print(f"\n{'='*85}")
    print(f"{'Condition':<25} {'S':>8} {'σ_S':>6} {'μ₁':>8} {'σ_μ₁':>6} {'⟨tri⟩':>7} {'⟨|α|⟩':>7} {'f_S':>6}")
    print(f"{'-'*85}")
    for cond in conditions:
        S_vals = [r['S'] for r in all_results[cond]]
        mu1_vals = [r['mu1'] for r in all_results[cond]]
        tri_vals = [r['mean_tri'] for r in all_results[cond]]
        alpha_vals = [r['mean_alpha'] for r in all_results[cond]]
        fs_vals = [r['frac_spatial'] for r in all_results[cond]]
        print(f"{labels[cond]:<25} {np.mean(S_vals):>8.3f} {np.std(S_vals):>6.3f} "
              f"{np.mean(mu1_vals):>8.3f} {np.std(mu1_vals):>6.3f} "
              f"{np.mean(tri_vals):>7.2f} {np.mean(alpha_vals):>7.3f} "
              f"{np.mean(fs_vals):>6.2f}")
    print(f"{'='*85}")
    
    # ── Save data for figure ──
    import json
    data = {}
    for cond in conditions:
        data[cond] = {
            'S': [r['S'] for r in all_results[cond]],
            'mu1': [r['mu1'] for r in all_results[cond]],
            'mean_alpha': [r['mean_alpha'] for r in all_results[cond]],
            'mean_tri': [r['mean_tri'] for r in all_results[cond]],
            'frac_spatial': [r['frac_spatial'] for r in all_results[cond]]
        }
    
    with open('/home/claude/dd_results.json', 'w') as f:
        json.dump(data, f)
    
    print("\nData saved → dd_results.json")
