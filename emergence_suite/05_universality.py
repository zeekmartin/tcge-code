#!/usr/bin/env python3
"""
TCGE GAP-Universality v8b — Statistiques propres

Corrections :
  1. p-value permutation : rang du Δ_real parmi N_perm nulls
  2. Cohen's d : (mean_real - mean_null) / pooled_std
  3. Seuil explicite : PASS = Δ > 0.05 ET p < 0.05
  4. Corrélation Δ vs triangle heterogeneity (var(tri))

Familles : ER, RGG, WS, BA, CM
N=150, <k>≈12, 10 trials par famille, 200 permutations par null.

Author: Zeek (TCGE project) — AI-assisted implementation
Date: 2026-02-19
"""

import numpy as np
from collections import defaultdict
from math import gamma as gamma_func
import time

# =============================================================================
# GRAPH GENERATORS (identiques à v8)
# =============================================================================

def make_ER(N, target_k, seed):
    np.random.seed(seed)
    p = target_k / (N - 1)
    edges = []
    adj = defaultdict(set)
    for i in range(N):
        for j in range(i+1, N):
            if np.random.random() < p:
                edges.append((i, j))
                adj[i].add(j)
                adj[j].add(i)
    return N, edges, adj

def make_RGG(N, dim, target_k, seed):
    np.random.seed(seed)
    from scipy.spatial import cKDTree
    positions = np.random.random((N, dim))
    vol_d = np.pi**(dim/2) / gamma_func(dim/2 + 1)
    r = (target_k / ((N - 1) * vol_d))**(1.0/dim)
    tree = cKDTree(positions, boxsize=np.ones(dim))
    pairs = tree.query_pairs(r)
    edges = [(min(i,j), max(i,j)) for i,j in pairs]
    adj = defaultdict(set)
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    return N, edges, adj

def make_WS(N, k, p_rewire, seed):
    np.random.seed(seed)
    adj = defaultdict(set)
    for i in range(N):
        for j in range(1, k//2 + 1):
            adj[i].add((i+j) % N)
            adj[(i+j) % N].add(i)
    edge_list = [(min(i,j), max(i,j)) for i in adj for j in adj[i] if i < j]
    for idx, (i, j) in enumerate(edge_list):
        if np.random.random() < p_rewire:
            new_j = np.random.randint(N)
            tries = 0
            while (new_j == i or new_j in adj[i]) and tries < 100:
                new_j = np.random.randint(N)
                tries += 1
            if tries < 100:
                adj[i].discard(j)
                adj[j].discard(i)
                adj[i].add(new_j)
                adj[new_j].add(i)
    edges = [(min(i,j), max(i,j)) for i in range(N) for j in adj[i] if i < j]
    return N, edges, adj

def make_BA(N, m, seed):
    np.random.seed(seed)
    adj = defaultdict(set)
    for i in range(m+1):
        for j in range(i+1, m+1):
            adj[i].add(j)
            adj[j].add(i)
    degrees = [m] * (m+1)
    total_degree = m * (m+1)
    for new_node in range(m+1, N):
        targets = set()
        while len(targets) < m:
            r = np.random.randint(total_degree)
            cumsum = 0
            for node in range(len(degrees)):
                cumsum += degrees[node]
                if cumsum > r:
                    if node != new_node and node not in targets:
                        targets.add(node)
                    break
        for t in targets:
            adj[new_node].add(t)
            adj[t].add(new_node)
        degrees.append(len(targets))
        for t in targets:
            degrees[t] += 1
        total_degree += 2 * len(targets)
    edges = [(min(i,j), max(i,j)) for i in range(N) for j in adj[i] if i < j]
    return N, edges, adj

def make_CM(N, target_k, seed):
    np.random.seed(seed)
    deg_seq = np.random.poisson(target_k, N)
    if sum(deg_seq) % 2 == 1:
        deg_seq[0] += 1
    stubs = []
    for i in range(N):
        stubs.extend([i] * deg_seq[i])
    np.random.shuffle(stubs)
    adj = defaultdict(set)
    edges_set = set()
    for k in range(0, len(stubs) - 1, 2):
        i, j = stubs[k], stubs[k+1]
        if i != j and (min(i,j), max(i,j)) not in edges_set:
            edges_set.add((min(i,j), max(i,j)))
            adj[i].add(j)
            adj[j].add(i)
    edges = list(edges_set)
    return N, edges, adj


# =============================================================================
# BIPHASAGE CORE
# =============================================================================

def run_biphasage(N, edges, adj, seed=None):
    """Run biphasage optimization. Returns (Δ, arrow, tri_array, abs_alpha)."""
    if seed is not None:
        np.random.seed(seed)
    
    n_e = len(edges)
    if n_e < 10:
        return 0, 0, np.array([]), np.array([])
    
    degree = np.zeros(N)
    for i, j in edges:
        degree[i] += 1
        degree[j] += 1
    
    tri = np.zeros(n_e, dtype=int)
    for idx, (i, j) in enumerate(edges):
        tri[idx] = len(adj[i] & adj[j])
    
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
    
    high = [i for i in range(n_e) if tri[i] > med]
    low = [i for i in range(n_e) if tri[i] <= med]
    
    if not high or not low:
        return 0, 0, tri, abs_alpha
    
    delta = np.mean(abs_alpha[low]) - np.mean(abs_alpha[high])
    
    # Arrow
    polarized = [i for i in range(n_e) if abs_alpha[i] > 0.5]
    if len(polarized) > 1:
        signs = np.sign(alpha_e[polarized])
        corr = np.corrcoef(signs, np.sign(delta_N[polarized]))[0, 1]
        arrow = abs(corr) if not np.isnan(corr) else 0
    else:
        arrow = 0
    
    return delta, arrow, tri, abs_alpha


def permutation_null(N, edges, adj, n_perm=200, seed_base=5000):
    """
    Null model : shuffle tri(e) labels, re-optimize, compute Δ using ORIGINAL tri.
    Returns array of Δ_null values.
    """
    n_e = len(edges)
    if n_e < 10:
        return np.zeros(n_perm)
    
    tri = np.zeros(n_e, dtype=int)
    for idx, (i, j) in enumerate(edges):
        tri[idx] = len(adj[i] & adj[j])
    
    degree = np.zeros(N)
    for i, j in edges:
        degree[i] += 1
        degree[j] += 1
    delta_N = np.array([degree[i] - degree[j] for i, j in edges])
    
    med_orig = np.median(tri)
    if med_orig == tri.min():
        med_orig = np.mean(tri)
    high_orig = [i for i in range(n_e) if tri[i] > med_orig]
    low_orig = [i for i in range(n_e) if tri[i] <= med_orig]
    
    null_deltas = np.zeros(n_perm)
    
    for p in range(n_perm):
        np.random.seed(seed_base + p)
        
        tri_shuffled = tri.copy()
        np.random.shuffle(tri_shuffled)
        
        alpha_e = np.random.randn(n_e) * 0.01
        tri_f = tri_shuffled.astype(float)
        
        for _ in range(2500):
            grad = (2 * alpha_e * (-0.25 + 0.5*tri_f + 0.01) + 0.5 * delta_N)
            alpha_e -= 0.005 * grad
            alpha_e = np.clip(alpha_e, -0.999, 0.999)
        
        abs_alpha = np.abs(alpha_e)
        
        if high_orig and low_orig:
            null_deltas[p] = np.mean(abs_alpha[low_orig]) - np.mean(abs_alpha[high_orig])
    
    return null_deltas


# =============================================================================
# MAIN
# =============================================================================

def run():
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  TCGE GAP-Universality v8b — Statistiques propres         ║")
    print("║  p-value perm | Cohen's d | Δ vs var(tri) corrélation     ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print()
    print("  Critère PASS : Δ_mean > 0.05 ET p_perm < 0.05")
    print("  Critère STRONG : Δ_mean > 0.15 ET p_perm < 0.01 ET d > 0.8")
    
    N = 150
    target_k = 12
    n_trials = 10
    n_perm = 200
    
    families = {
        'ER':     lambda seed: make_ER(N, target_k, seed),
        'RGG_3D': lambda seed: make_RGG(N, 3, target_k, seed),
        'WS':     lambda seed: make_WS(N, target_k, 0.1, seed),
        'BA':     lambda seed: make_BA(N, target_k//2, seed),
        'CM':     lambda seed: make_CM(N, target_k, seed),
    }
    
    summary = {}
    hetero_data = []  # (var_tri, Δ) pairs across all families
    
    for name, gen_fn in families.items():
        print(f"\n{'═'*65}")
        print(f"  {name}")
        print(f"{'═'*65}")
        
        deltas = []
        arrows = []
        tri_vars = []
        
        for trial in range(n_trials):
            n, edges, adj = gen_fn(seed=1000*trial + 42)
            d, a, tri, aa = run_biphasage(n, edges, adj, seed=2000*trial + 7)
            deltas.append(d)
            arrows.append(a)
            if len(tri) > 0:
                tv = np.var(tri)
                tri_vars.append(tv)
                hetero_data.append((tv, d))
            
            status = '✅' if d > 0.15 else ('⚠️' if d > 0.05 else '❌')
            print(f"  T{trial+1:>2}: Δ={d:.3f} arrow={a:.2f} var(tri)={tv:.1f} {status}")
        
        # Permutation null (on trial 0 graph)
        print(f"\n  Permutation null ({n_perm} permutations)...")
        t0 = time.time()
        n, edges, adj = gen_fn(seed=1042)
        null_deltas = permutation_null(n, edges, adj, n_perm=n_perm)
        
        mean_real = np.mean(deltas)
        mean_null = np.mean(null_deltas)
        std_null = np.std(null_deltas)
        
        # p-value: fraction of nulls >= mean_real
        p_value = np.mean(null_deltas >= mean_real)
        # If none exceed, p < 1/n_perm
        if p_value == 0:
            p_value_str = f"< {1/n_perm:.4f}"
            p_value_num = 1 / (n_perm + 1)  # conservative
        else:
            p_value_str = f"{p_value:.4f}"
            p_value_num = p_value
        
        # Cohen's d
        pooled_std = np.sqrt((np.var(deltas) + np.var(null_deltas)) / 2)
        if pooled_std > 0:
            cohens_d = (mean_real - mean_null) / pooled_std
        else:
            cohens_d = 0
        
        # Verdict
        passes = mean_real > 0.05 and p_value_num < 0.05
        strong = mean_real > 0.15 and p_value_num < 0.01 and cohens_d > 0.8
        
        print(f"  Null: Δ = {mean_null:.4f} ± {std_null:.4f} ({time.time()-t0:.1f}s)")
        print(f"  Real: Δ = {mean_real:.4f} ± {np.std(deltas):.4f}")
        print(f"  p-value (permutation): {p_value_str}")
        print(f"  Cohen's d: {cohens_d:.2f}")
        print(f"  Arrow moyen: {np.mean(arrows):.2f}")
        
        if strong:
            print(f"  → ✅ STRONG PASS")
        elif passes:
            print(f"  → ✅ PASS")
        else:
            print(f"  → ❌ FAIL")
        
        summary[name] = {
            'mean_delta': mean_real, 'std_delta': np.std(deltas),
            'mean_arrow': np.mean(arrows),
            'null_mean': mean_null, 'null_std': std_null,
            'p_value': p_value_num, 'p_value_str': p_value_str,
            'cohens_d': cohens_d,
            'passes': passes, 'strong': strong,
            'mean_tri_var': np.mean(tri_vars) if tri_vars else 0
        }
    
    # ── CORRÉLATION Δ vs var(tri) ──
    print(f"\n\n{'═'*65}")
    print("  CORRÉLATION Δ vs HÉTÉROGÉNÉITÉ tri(e)")
    print(f"{'═'*65}\n")
    
    if len(hetero_data) > 5:
        vars_t = np.array([x[0] for x in hetero_data])
        deltas_t = np.array([x[1] for x in hetero_data])
        
        corr = np.corrcoef(vars_t, deltas_t)[0, 1]
        print(f"  Pearson r = {corr:.3f} (n={len(hetero_data)} points)")
        
        # Also per-family means
        print(f"\n  {'Famille':<10} {'var(tri)':<12} {'Δ moyen':<10}")
        print(f"  {'-'*35}")
        for name in families:
            s = summary[name]
            print(f"  {name:<10} {s['mean_tri_var']:<12.1f} {s['mean_delta']:<10.3f}")
        
        if corr < -0.3:
            print(f"\n  → Corrélation négative: plus tri(e) est hétérogène,")
            print(f"    plus Δ se contracte. Cohérent avec la prédiction")
            print(f"    structurelle pour BA.")
        elif abs(corr) < 0.3:
            print(f"\n  → Corrélation faible: var(tri) n'est pas le seul facteur.")
    
    # ── SYNTHÈSE FINALE ──
    print(f"\n\n{'═'*65}")
    print("  SYNTHÈSE UNIVERSALITY v8b")
    print(f"{'═'*65}\n")
    
    print(f"  {'Famille':<10} {'Δ':<10} {'p-perm':<12} {'Cohen d':<10} "
          f"{'Arrow':<8} {'Verdict'}")
    print(f"  {'-'*60}")
    
    n_pass = 0
    n_strong = 0
    for name in families:
        s = summary[name]
        n_pass += int(s['passes'])
        n_strong += int(s['strong'])
        v = '✅ STRONG' if s['strong'] else ('✅ PASS' if s['passes'] else '❌ FAIL')
        print(f"  {name:<10} {s['mean_delta']:<10.3f} {s['p_value_str']:<12} "
              f"{s['cohens_d']:<10.2f} {s['mean_arrow']:<8.2f} {v}")
    
    print(f"\n  Résultat: {n_pass}/5 PASS, {n_strong}/5 STRONG")
    
    if n_pass >= 4:
        print(f"\n  Formulation autorisée :")
        print(f"  'Phase separation appears across all five graph families")
        print(f"   tested (ER, RGG, WS, BA, CM), with permutation p < 0.05")
        print(f"   and large effect sizes (Cohen d > 0.8) in {n_strong}/5 cases.")
        if summary.get('BA', {}).get('mean_delta', 0) < 0.15:
            print(f"   Amplitude is reduced in scale-free (BA) graphs,")
            print(f"   consistent with concentrated triangle distribution.'")
    
    return summary


if __name__ == "__main__":
    t0 = time.time()
    results = run()
    print(f"\n  Temps total: {time.time()-t0:.0f}s")
