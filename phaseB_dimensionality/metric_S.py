#!/usr/bin/env python3
"""
TCGE — Métrique de séparation S indépendante de tri(e)

Problème : Δ utilise tri(e) comme classifieur, donc s'effondre quand tri→0.
Solution : fit mixture à 2 composantes sur |α|, sans aucune info de tri.

Méthode :
  1. Fit Beta mixture (2 composantes) sur |α| via EM
  2. Extraire : μ₁, μ₂ (moyennes), σ₁, σ₂ (std), p (poids composante 1)
  3. Score de séparation :
     S = |μ₁ - μ₂| × 2×min(p, 1-p)
     → S élevé ssi deux pics éloignés ET populations équilibrées
  4. Ashman D = |μ₁ - μ₂| / sqrt(σ₁² + σ₂²)
     → D > 2 = vraie bimodalité

Validation :
  Part 1 — ER graphs, comparer S vs Δ (doivent corréler)
  Part 2 — Rewiring series (le test décisif : S quand tri→0)

Author: Zeek (TCGE project) — AI-assisted implementation
Date: 2026-02-20
"""

import numpy as np
from collections import defaultdict
from math import gamma as gamma_func
from scipy.spatial import cKDTree
from scipy.special import beta as beta_fn
from scipy.optimize import minimize
import time


# =============================================================================
# BETA MIXTURE EM
# =============================================================================

def beta_pdf(x, a, b):
    """Beta distribution PDF, numerically stable."""
    x = np.clip(x, 1e-10, 1 - 1e-10)
    log_pdf = (a - 1) * np.log(x) + (b - 1) * np.log(1 - x)
    from scipy.special import betaln
    log_pdf -= betaln(a, b)
    return np.exp(log_pdf)


def fit_beta_mixture(data, max_iter=100, tol=1e-6):
    """
    EM algorithm for 2-component Beta mixture on [0, 1] data.
    Returns: (mu1, mu2, std1, std2, pi1, converged)
    """
    data = np.clip(data, 1e-6, 1 - 1e-6)
    n = len(data)
    
    # Initialize: split at median
    med = np.median(data)
    mask_lo = data <= med
    mask_hi = data > med
    
    # Method of moments for initial Beta params
    def mom_beta(x):
        m, v = np.mean(x), np.var(x)
        if v >= m * (1 - m):
            v = m * (1 - m) * 0.9
        common = m * (1 - m) / v - 1
        a = m * common
        b = (1 - m) * common
        return max(a, 0.5), max(b, 0.5)
    
    if sum(mask_lo) > 2 and sum(mask_hi) > 2:
        a1, b1 = mom_beta(data[mask_lo])
        a2, b2 = mom_beta(data[mask_hi])
    else:
        a1, b1 = 2, 5
        a2, b2 = 5, 2
    
    pi1 = 0.5
    
    for iteration in range(max_iter):
        # E-step
        p1 = pi1 * beta_pdf(data, a1, b1)
        p2 = (1 - pi1) * beta_pdf(data, a2, b2)
        total = p1 + p2 + 1e-300
        gamma1 = p1 / total
        gamma2 = p2 / total
        
        # M-step: pi
        n1 = np.sum(gamma1)
        n2 = np.sum(gamma2)
        pi1_new = n1 / n
        
        # M-step: Beta params via weighted method of moments
        def weighted_beta(gamma_k):
            w = gamma_k / (np.sum(gamma_k) + 1e-300)
            m = np.sum(w * data)
            v = np.sum(w * (data - m)**2)
            if v < 1e-8 or v >= m * (1 - m):
                v = m * (1 - m) * 0.5
            common = m * (1 - m) / v - 1
            common = max(common, 1.0)
            a = max(m * common, 0.5)
            b = max((1 - m) * common, 0.5)
            return a, b
        
        a1_new, b1_new = weighted_beta(gamma1)
        a2_new, b2_new = weighted_beta(gamma2)
        
        # Check convergence
        change = (abs(pi1_new - pi1) + abs(a1_new - a1) + abs(b1_new - b1)
                  + abs(a2_new - a2) + abs(b2_new - b2))
        
        pi1, a1, b1, a2, b2 = pi1_new, a1_new, b1_new, a2_new, b2_new
        
        if change < tol:
            break
    
    # Compute means and stds
    mu1 = a1 / (a1 + b1)
    mu2 = a2 / (a2 + b2)
    std1 = np.sqrt(a1 * b1 / ((a1 + b1)**2 * (a1 + b1 + 1)))
    std2 = np.sqrt(a2 * b2 / ((a2 + b2)**2 * (a2 + b2 + 1)))
    
    # Ensure mu1 < mu2 (convention: component 1 = low, 2 = high)
    if mu1 > mu2:
        mu1, mu2 = mu2, mu1
        std1, std2 = std2, std1
        pi1 = 1 - pi1
    
    return mu1, mu2, std1, std2, pi1, iteration < max_iter - 1


def separation_score(data):
    """
    Compute S and Ashman D from |α| distribution.
    Returns dict with all metrics.
    """
    mu1, mu2, std1, std2, pi1, converged = fit_beta_mixture(data)
    
    distance = abs(mu2 - mu1)
    balance = 2 * min(pi1, 1 - pi1)  # 1 if 50/50, 0 if one-sided
    S = distance * balance
    
    # Ashman D
    denom = np.sqrt(std1**2 + std2**2)
    D = distance / denom if denom > 0 else 0
    
    return {
        'S': S, 'D': D,
        'mu1': mu1, 'mu2': mu2,
        'std1': std1, 'std2': std2,
        'pi1': pi1, 'distance': distance,
        'balance': balance, 'converged': converged
    }


# =============================================================================
# GRAPH + BIPHASAGE (compact)
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
                adj[i].add(j); adj[j].add(i)
    return N, edges, adj


def make_RGG(N, dim, target_k, seed):
    np.random.seed(seed)
    positions = np.random.random((N, dim))
    vol_d = np.pi**(dim/2) / gamma_func(dim/2 + 1)
    r = (target_k / ((N - 1) * vol_d))**(1.0/dim)
    tree = cKDTree(positions, boxsize=np.ones(dim))
    pairs = tree.query_pairs(r)
    edges = [(min(i,j), max(i,j)) for i,j in pairs]
    adj = defaultdict(set)
    for i, j in edges:
        adj[i].add(j); adj[j].add(i)
    return N, edges, adj


def run_biphasage(N, edges, adj, seed=None):
    if seed is not None:
        np.random.seed(seed)
    n_e = len(edges)
    if n_e < 10:
        return np.array([]), np.array([]), 0
    
    degree = np.zeros(N)
    for i, j in edges:
        degree[i] += 1; degree[j] += 1
    
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
    
    # Classic Δ
    med = np.median(tri)
    if med == tri.min(): med = np.mean(tri)
    low = [i for i in range(n_e) if tri[i] <= med]
    high = [i for i in range(n_e) if tri[i] > med]
    delta = (np.mean(abs_alpha[low]) - np.mean(abs_alpha[high])) if low and high else 0
    
    return abs_alpha, tri, delta


def degree_preserving_rewire(edges, adj, n_nodes, n_swaps, seed=None):
    if seed is not None:
        np.random.seed(seed)
    edges = list(edges)
    edge_set = set(edges)
    adj = defaultdict(set)
    for i, j in edges:
        adj[i].add(j); adj[j].add(i)
    n_e = len(edges)
    done = 0; attempts = 0
    while done < n_swaps and attempts < n_swaps * 10:
        attempts += 1
        i1, i2 = np.random.randint(n_e), np.random.randint(n_e)
        if i1 == i2: continue
        a,b = edges[i1]; c,d = edges[i2]
        if np.random.random() < 0.5:
            n1,n2 = (min(a,d),max(a,d)), (min(c,b),max(c,b))
        else:
            n1,n2 = (min(a,c),max(a,c)), (min(b,d),max(b,d))
        if n1[0]==n1[1] or n2[0]==n2[1]: continue
        if n1 in edge_set or n2 in edge_set: continue
        if n1==n2: continue
        edge_set.discard(edges[i1]); edge_set.discard(edges[i2])
        adj[a].discard(b); adj[b].discard(a)
        adj[c].discard(d); adj[d].discard(c)
        edges[i1]=n1; edges[i2]=n2
        edge_set.add(n1); edge_set.add(n2)
        adj[n1[0]].add(n1[1]); adj[n1[1]].add(n1[0])
        adj[n2[0]].add(n2[1]); adj[n2[1]].add(n2[0])
        done += 1
    return edges, adj


# =============================================================================
# PART 1: VALIDATION — S vs Δ on standard graphs
# =============================================================================

def part1_validation():
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  PART 1 — Validation : S concorde avec Δ ?                ║")
    print("║  ER (N=150, k=12) + RGG d=3 (N=3000, k=8), 10 trials    ║")
    print("╚══════════════════════════════════════════════════════════════╝\n")
    
    configs = [
        ('ER', lambda seed: make_ER(150, 12, seed)),
        ('RGG_3D', lambda seed: make_RGG(3000, 3, 8, seed)),
    ]
    
    all_deltas = []
    all_S = []
    all_D = []
    
    for name, gen_fn in configs:
        print(f"  {name}:")
        print(f"  {'trial':<7} {'Δ':<8} {'S':<8} {'D':<8} "
              f"{'μ₁':<7} {'μ₂':<7} {'π₁':<7} {'bal':<6}")
        print(f"  {'-'*60}")
        
        for trial in range(10):
            n, edges, adj = gen_fn(1000*trial + 42)
            aa, tri, delta = run_biphasage(n, edges, adj, seed=2000*trial + 7)
            
            if len(aa) < 10:
                continue
            
            metrics = separation_score(aa)
            
            all_deltas.append(delta)
            all_S.append(metrics['S'])
            all_D.append(metrics['D'])
            
            print(f"  {trial+1:<7} {delta:<8.3f} {metrics['S']:<8.3f} "
                  f"{metrics['D']:<8.2f} {metrics['mu1']:<7.2f} "
                  f"{metrics['mu2']:<7.2f} {metrics['pi1']:<7.2f} "
                  f"{metrics['balance']:<6.2f}")
        print()
    
    # Correlation
    r_S = np.corrcoef(all_deltas, all_S)[0, 1]
    r_D = np.corrcoef(all_deltas, all_D)[0, 1]
    
    print(f"  Corrélations (n={len(all_deltas)}):")
    print(f"    Δ vs S : r = {r_S:.3f}")
    print(f"    Δ vs D : r = {r_D:.3f}")
    
    if r_S > 0.8:
        print(f"    → ✅ S concorde avec Δ (validation réussie)")
    elif r_S > 0.5:
        print(f"    → ⚠️ Corrélation modérée")
    else:
        print(f"    → ❌ S ne concorde pas avec Δ")
    
    return r_S, r_D


# =============================================================================
# PART 2: REWIRING — S quand Δ échoue
# =============================================================================

def part2_rewiring():
    print(f"\n\n╔══════════════════════════════════════════════════════════════╗")
    print(f"║  PART 2 — Rewiring : S quand tri→0                        ║")
    print(f"║  RGG d=3, N=3000, k=8, 6 trials par niveau              ║")
    print(f"╚══════════════════════════════════════════════════════════════╝\n")
    
    N = 3000
    n_trials = 6
    rewire_fracs = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.70, 0.90]
    
    results = []
    
    for frac in rewire_fracs:
        deltas, Ss, Ds, tris = [], [], [], []
        mu1s, mu2s, pi1s = [], [], []
        
        for trial in range(n_trials):
            n, edges, adj = make_RGG(N, 3, 8, seed=1000*trial + 42)
            edges = list(edges)
            adj_rw = defaultdict(set)
            for i, j in edges: adj_rw[i].add(j); adj_rw[j].add(i)
            
            n_swaps = int(frac * len(edges))
            if n_swaps > 0:
                edges, adj_rw = degree_preserving_rewire(
                    edges, adj_rw, N, n_swaps, seed=3000*trial + int(frac*1000))
            
            aa, tri, delta = run_biphasage(N, edges, adj_rw, seed=2000*trial + 7)
            
            if len(aa) < 10:
                continue
            
            metrics = separation_score(aa)
            
            deltas.append(delta)
            Ss.append(metrics['S'])
            Ds.append(metrics['D'])
            tris.append(np.mean(tri))
            mu1s.append(metrics['mu1'])
            mu2s.append(metrics['mu2'])
            pi1s.append(metrics['pi1'])
        
        results.append({
            'frac': frac,
            'delta': np.mean(deltas), 'S': np.mean(Ss), 'D': np.mean(Ds),
            'tri': np.mean(tris),
            'mu1': np.mean(mu1s), 'mu2': np.mean(mu2s), 'pi1': np.mean(pi1s),
        })
    
    print(f"  {'rew%':<6} {'⟨tri⟩':<7} {'Δ':<8} {'S':<8} {'D':<8} "
          f"{'μ₁':<7} {'μ₂':<7} {'π₁':<7}")
    print(f"  {'-'*62}")
    
    for r in results:
        print(f"  {r['frac']*100:>4.0f}%  {r['tri']:<7.2f} {r['delta']:<8.3f} "
              f"{r['S']:<8.3f} {r['D']:<8.2f} {r['mu1']:<7.2f} "
              f"{r['mu2']:<7.2f} {r['pi1']:<7.2f}")
    
    # ── ANALYSIS ──
    print(f"\n  {'═'*55}")
    print(f"  Comparaison Δ vs S sur la série de rewiring")
    print(f"  {'═'*55}\n")
    
    deltas_arr = np.array([r['delta'] for r in results])
    S_arr = np.array([r['S'] for r in results])
    D_arr = np.array([r['D'] for r in results])
    
    # Where do they diverge?
    print(f"  {'rew%':<6} {'Δ':<8} {'S':<8} {'Δ/Δ_max':<10} {'S/S_max':<10} {'diverge?'}")
    print(f"  {'-'*55}")
    
    d_max = max(deltas_arr)
    s_max = max(S_arr)
    
    for r in results:
        d_norm = r['delta'] / d_max
        s_norm = r['S'] / s_max
        diverge = abs(d_norm - s_norm) > 0.15
        print(f"  {r['frac']*100:>4.0f}%  {r['delta']:<8.3f} {r['S']:<8.3f} "
              f"{d_norm:<10.2f} {s_norm:<10.2f} {'← DIVERGE' if diverge else ''}")
    
    return results


# =============================================================================
# MAIN
# =============================================================================

def run():
    t0 = time.time()
    
    r_S, r_D = part1_validation()
    results = part2_rewiring()
    
    # ── VERDICT ──
    print(f"\n\n{'═'*65}")
    print("  VERDICT")
    print(f"{'═'*65}\n")
    
    S_arr = [r['S'] for r in results]
    D_arr = [r['D'] for r in results]
    delta_arr = [r['delta'] for r in results]
    
    # Does S stay high where Δ drops?
    # Compare at 30% rewiring (where Δ drops sharply)
    r30 = [r for r in results if r['frac'] == 0.3][0]
    r0 = [r for r in results if r['frac'] == 0.0][0]
    
    delta_retention_30 = r30['delta'] / r0['delta']
    S_retention_30 = r30['S'] / r0['S'] if r0['S'] > 0 else 0
    
    print(f"  À 30% rewiring (point de divergence Δ) :")
    print(f"    Δ : {r0['delta']:.3f} → {r30['delta']:.3f} (rétention {delta_retention_30:.0%})")
    print(f"    S : {r0['S']:.3f} → {r30['S']:.3f} (rétention {S_retention_30:.0%})")
    print()
    
    if S_retention_30 > delta_retention_30 + 0.1:
        print(f"  ✅ S résiste mieux que Δ au rewiring.")
        print(f"     S capture la séparation réelle,")
        print(f"     Δ sous-estime quand tri perd son pouvoir discriminant.")
    elif abs(S_retention_30 - delta_retention_30) < 0.1:
        print(f"  ⚠️ S et Δ se comportent de façon similaire.")
        print(f"     Les deux métriques mesurent le même phénomène.")
    else:
        print(f"  ❌ S chute plus vite que Δ (inattendu).")
    
    # Summary of what S reveals at high rewiring
    r90 = [r for r in results if r['frac'] == 0.9][0]
    print(f"\n  À 90% rewiring (tri ≈ 0) :")
    print(f"    S = {r90['S']:.3f}, D = {r90['D']:.2f}")
    print(f"    μ₁ = {r90['mu1']:.2f}, μ₂ = {r90['mu2']:.2f}, π₁ = {r90['pi1']:.2f}")
    
    if r90['D'] > 2:
        print(f"    → Bimodalité FORTE (D > 2) même sans triangles")
    elif r90['D'] > 1:
        print(f"    → Bimodalité modérée (D > 1)")
    else:
        print(f"    → Pas de bimodalité (D < 1) — distribution unimodale")
    
    print(f"\n  Temps total: {time.time()-t0:.0f}s")
    
    return results


if __name__ == "__main__":
    results = run()
