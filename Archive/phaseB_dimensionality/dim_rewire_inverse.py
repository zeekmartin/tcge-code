#!/usr/bin/env python3
"""
TCGE — Rewiring inverse : créer des triangles → créer la séparation T/S

Preuve miroir du rôle constitutif du protecteur.

Protocole :
  1. Point de départ : ER (N=3000, ⟨k⟩=8) — quasi pas de triangles
  2. Opération : triangle-closing degree-preserving rewiring
     - Choisir deux arêtes (a-b, c-d)
     - Rewire → (a-c, b-d) si ça AUGMENTE le nombre de triangles
     - Degrés strictement conservés
  3. Mesures à chaque palier :
     - Mixture model β : π₁, μ₁, μ₂, S
     - Δ (standard tri-based)  
     - Ashman D
     - ⟨tri⟩, ⟨|α|⟩

Attendu :
  - début : 1 population, S ≈ 0 (polarisation uniforme)
  - puis : π₁ > 0, S monte → proto-spatial apparaît
  - optimum possible à ⟨tri⟩ intermédiaire

Author: Zeek (TCGE project) — AI-assisted implementation
Date: 2026-02-20
"""

import numpy as np
from collections import defaultdict
from scipy.special import betaln
import time


# =============================================================================
# BETA MIXTURE EM (from tcge_metric_S.py)
# =============================================================================

def beta_pdf(x, a, b):
    x = np.clip(x, 1e-10, 1 - 1e-10)
    log_pdf = (a - 1) * np.log(x) + (b - 1) * np.log(1 - x) - betaln(a, b)
    return np.exp(log_pdf)

def fit_beta_mixture(data, max_iter=100, tol=1e-6):
    data = np.clip(data, 1e-6, 1 - 1e-6)
    n = len(data)
    
    med = np.median(data)
    mask_lo = data <= med
    mask_hi = data > med
    
    def mom_beta(x):
        m, v = np.mean(x), np.var(x)
        if v >= m * (1 - m):
            v = m * (1 - m) * 0.9
        common = m * (1 - m) / v - 1
        return max(m * common, 0.5), max((1 - m) * common, 0.5)
    
    if sum(mask_lo) > 2 and sum(mask_hi) > 2:
        a1, b1 = mom_beta(data[mask_lo])
        a2, b2 = mom_beta(data[mask_hi])
    else:
        a1, b1 = 2, 5
        a2, b2 = 5, 2
    
    pi1 = 0.5
    
    for iteration in range(max_iter):
        pdf1 = beta_pdf(data, a1, b1)
        pdf2 = beta_pdf(data, a2, b2)
        denom = pi1 * pdf1 + (1 - pi1) * pdf2 + 1e-300
        gamma = pi1 * pdf1 / denom
        
        pi1_new = np.mean(gamma)
        pi1_new = np.clip(pi1_new, 0.01, 0.99)
        
        def weighted_beta(g):
            w = g / (g.sum() + 1e-300)
            m = np.sum(w * data)
            v = np.sum(w * (data - m)**2)
            m = np.clip(m, 0.01, 0.99)
            if v >= m * (1 - m):
                v = m * (1 - m) * 0.9
            if v < 1e-10:
                v = 1e-4
            common = m * (1 - m) / v - 1
            return max(m * common, 0.1), max((1 - m) * common, 0.1)
        
        a1_new, b1_new = weighted_beta(gamma)
        a2_new, b2_new = weighted_beta(1 - gamma)
        
        change = abs(pi1_new - pi1) + abs(a1_new - a1) + abs(a2_new - a2)
        a1, b1, a2, b2, pi1 = a1_new, b1_new, a2_new, b2_new, pi1_new
        
        if change < tol:
            break
    
    mu1 = a1 / (a1 + b1)
    mu2 = a2 / (a2 + b2)
    var1 = a1 * b1 / ((a1 + b1)**2 * (a1 + b1 + 1))
    var2 = a2 * b2 / ((a2 + b2)**2 * (a2 + b2 + 1))
    
    if mu1 > mu2:
        mu1, mu2 = mu2, mu1
        var1, var2 = var2, var1
        pi1 = 1 - pi1
    
    return mu1, mu2, np.sqrt(var1), np.sqrt(var2), pi1, True


def separation_score(data):
    mu1, mu2, std1, std2, pi1, converged = fit_beta_mixture(data)
    distance = abs(mu2 - mu1)
    balance = 2 * min(pi1, 1 - pi1)
    S = distance * balance
    denom = np.sqrt(std1**2 + std2**2)
    D = distance / denom if denom > 0 else 0
    return {
        'S': S, 'D': D,
        'mu1': mu1, 'mu2': mu2,
        'std1': std1, 'std2': std2,
        'pi1': pi1, 'distance': distance,
        'balance': balance
    }


# =============================================================================
# ER GRAPH
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
    degree = np.array([len(adj[i]) for i in range(N)])
    return N, edges, adj, degree


# =============================================================================
# TRIANGLE-CLOSING REWIRE
# =============================================================================

def triangle_closing_rewire(edges, adj, n_nodes, n_attempts, seed=None):
    """
    Degree-preserving rewiring that ONLY accepts swaps that increase triangles.
    Returns edges, adj, number of successful swaps.
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
    
    for _ in range(n_attempts):
        idx1 = np.random.randint(n_e)
        idx2 = np.random.randint(n_e)
        if idx1 == idx2:
            continue
        
        a, b = edges[idx1]
        c, d = edges[idx2]
        
        # Count current triangles involving these edges
        tri_before = (len(adj[a] & adj[b]) + len(adj[c] & adj[d]))
        
        # Try swap
        if np.random.random() < 0.5:
            new1 = (min(a, d), max(a, d))
            new2 = (min(c, b), max(c, b))
        else:
            new1 = (min(a, c), max(a, c))
            new2 = (min(b, d), max(b, d))
        
        if new1[0] == new1[1] or new2[0] == new2[1]:
            continue
        if new1 in edge_set or new2 in edge_set:
            continue
        if new1 == new2:
            continue
        
        # Temporarily apply swap to count new triangles
        adj[a].discard(b); adj[b].discard(a)
        adj[c].discard(d); adj[d].discard(c)
        adj[new1[0]].add(new1[1]); adj[new1[1]].add(new1[0])
        adj[new2[0]].add(new2[1]); adj[new2[1]].add(new2[0])
        
        tri_after = (len(adj[new1[0]] & adj[new1[1]]) + 
                     len(adj[new2[0]] & adj[new2[1]]))
        
        if tri_after > tri_before:
            # Accept
            edge_set.discard(edges[idx1])
            edge_set.discard(edges[idx2])
            edges[idx1] = new1
            edges[idx2] = new2
            edge_set.add(new1)
            edge_set.add(new2)
            swaps_done += 1
        else:
            # Reject — revert
            adj[new1[0]].discard(new1[1]); adj[new1[1]].discard(new1[0])
            adj[new2[0]].discard(new2[1]); adj[new2[1]].discard(new2[0])
            adj[a].add(b); adj[b].add(a)
            adj[c].add(d); adj[d].add(c)
    
    return edges, adj, swaps_done


# =============================================================================
# BIPHASAGE + S
# =============================================================================

def full_analysis(N, edges, adj, seed=None):
    """Run biphasage, compute Δ and S."""
    if seed is not None:
        np.random.seed(seed)
    
    n_e = len(edges)
    if n_e < 10:
        return {'delta': 0, 'S': 0, 'alpha_mean': 0, 'tri_mean': 0}
    
    degree = np.zeros(N)
    for i, j in edges:
        degree[i] += 1
        degree[j] += 1
    
    tri = np.zeros(n_e, dtype=int)
    for idx, (i, j) in enumerate(edges):
        tri[idx] = len(adj[i] & adj[j])
    
    # Optimize
    alpha_e = np.random.randn(n_e) * 0.01
    delta_N = np.array([degree[i] - degree[j] for i, j in edges])
    tri_f = tri.astype(float)
    
    for _ in range(2500):
        grad = (2 * alpha_e * (-0.25 + 0.5*tri_f + 0.01) + 0.5 * delta_N)
        alpha_e -= 0.005 * grad
        alpha_e = np.clip(alpha_e, -0.999, 0.999)
    
    aa = np.abs(alpha_e)
    
    # Δ standard
    med = np.median(tri)
    if med == tri.min():
        med = np.mean(tri)
    low = [i for i in range(n_e) if tri[i] <= med]
    high = [i for i in range(n_e) if tri[i] > med]
    delta = (np.mean(aa[low]) - np.mean(aa[high])) if low and high else 0
    
    # S metric
    s_result = separation_score(aa)
    
    return {
        'delta': delta,
        'S': s_result['S'],
        'D': s_result['D'],
        'pi1': s_result['pi1'],
        'mu1': s_result['mu1'],
        'mu2': s_result['mu2'],
        'alpha_mean': np.mean(aa),
        'tri_mean': np.mean(tri),
        'tri_var': np.var(tri),
        'k_mean': np.mean(degree),
    }


# =============================================================================
# MAIN
# =============================================================================

def run():
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  Rewiring inverse : ER → ajout de triangles               ║")
    print("║  N=3000, ⟨k⟩=8, triangle-closing rewiring                ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    N = 3000
    target_k = 8
    n_trials = 5
    
    # Progressive triangle-closing: each step attempts more swaps
    # cumulative_attempts controls how many triangle-closing attempts
    steps = [0, 5000, 15000, 40000, 100000, 250000, 500000, 1000000]
    
    print(f"\n  {n_trials} trials, {len(steps)} paliers d'ajout de triangles")
    
    all_results = {s: [] for s in steps}
    
    for trial in range(n_trials):
        t0 = time.time()
        
        # Fresh ER
        n, edges, adj, degree = make_ER(N, target_k, seed=1000*trial + 42)
        
        prev_attempts = 0
        for step_attempts in steps:
            # Apply incremental triangle-closing
            additional = step_attempts - prev_attempts
            if additional > 0:
                edges, adj, swaps = triangle_closing_rewire(
                    edges, adj, N, additional, seed=5000*trial + step_attempts)
            
            # Measure
            r = full_analysis(N, edges, adj, seed=2000*trial + 7)
            all_results[step_attempts].append(r)
            prev_attempts = step_attempts
        
        print(f"  Trial {trial+1}: {time.time()-t0:.1f}s")
    
    # ── TABLE ──
    print(f"\n{'═'*85}")
    print(f"  {'attempts':<10} {'⟨tri⟩':<8} {'var(tri)':<9} {'⟨|α|⟩':<8} "
          f"{'Δ':<8} {'S':<8} {'π₁':<7} {'μ₁':<7} {'μ₂':<7} {'D':<7} {'⟨k⟩':<6}")
    print(f"  {'─'*83}")
    
    summary = []
    for step in steps:
        trials = all_results[step]
        row = {}
        for key in trials[0]:
            row[key] = np.mean([t[key] for t in trials])
        row['step'] = step
        summary.append(row)
        
        print(f"  {step:<10} {row['tri_mean']:<8.2f} {row['tri_var']:<9.2f} "
              f"{row['alpha_mean']:<8.3f} {row['delta']:<8.3f} "
              f"{row['S']:<8.3f} {row['pi1']:<7.2f} "
              f"{row['mu1']:<7.2f} {row['mu2']:<7.2f} "
              f"{row['D']:<7.2f} {row['k_mean']:<6.1f}")
    
    # ── ANALYSIS ──
    print(f"\n{'═'*85}")
    print("  ANALYSE")
    print(f"{'═'*85}\n")
    
    tris = [s['tri_mean'] for s in summary]
    Ss = [s['S'] for s in summary]
    pi1s = [s['pi1'] for s in summary]
    deltas_arr = [s['delta'] for s in summary]
    alphas = [s['alpha_mean'] for s in summary]
    
    # Does S track the creation of T/S separation?
    # Expected: S starts low (no separation), rises as triangles are added
    
    S_initial = Ss[0]
    S_max = max(Ss)
    S_max_idx = Ss.index(S_max)
    tri_at_max_S = tris[S_max_idx]
    
    print(f"  S initial (ER) = {S_initial:.3f}")
    print(f"  S maximum      = {S_max:.3f} (at ⟨tri⟩ = {tri_at_max_S:.2f})")
    print(f"  S final        = {Ss[-1]:.3f} (at ⟨tri⟩ = {tris[-1]:.2f})")
    print()
    
    # π₁ trajectory
    print(f"  π₁ trajectory (proto-spatial population):")
    for i, s in enumerate(summary):
        bar = '█' * int(s['pi1'] * 40)
        print(f"    ⟨tri⟩={s['tri_mean']:.2f}  π₁={s['pi1']:.2f}  {bar}")
    
    # ── VERDICT ──
    print(f"\n{'═'*85}")
    print("  VERDICT")
    print(f"{'═'*85}\n")
    
    if S_max > S_initial * 2 and S_max > 0.05:
        # Check: was there NO separation at start, and separation at peak?
        if pi1s[0] < 0.08 and max(pi1s) > 0.12:
            print("  ✅ PREUVE MIROIR CONFIRMÉE")
            print()
            print("  ER initial : polarisation uniforme, π₁ ≈ 0, S ≈ 0")
            print(f"  → pas de population proto-spatiale")
            print()
            print(f"  Après triangle-closing : π₁ = {max(pi1s):.2f}, "
                  f"S = {S_max:.3f}")
            print(f"  → une population proto-spatiale ÉMERGE")
            print()
            print("  Les triangles ne corrèlent pas avec la séparation.")
            print("  Ils la CRÉENT.")
        elif S_max > S_initial * 1.5:
            print("  ⚠️ SIGNAL PARTIEL")
            print(f"  S augmente ({S_initial:.3f} → {S_max:.3f})")
            print(f"  mais π₁ reste faible (max={max(pi1s):.2f})")
        else:
            print("  ⚠️ RÉSULTAT AMBIGU")
    else:
        print("  ❌ PAS DE PREUVE MIROIR")
        print(f"  S n'augmente pas significativement")
        print(f"  ({S_initial:.3f} → {S_max:.3f})")
    
    # Comparison: destruction vs creation
    print(f"\n  ── Symétrie destruction/création ──")
    print(f"  Destruction (RGG → rewire) : S pic à ⟨tri⟩ ≈ 2")
    print(f"  Création (ER → tri-close)  : S pic at ⟨tri⟩ ≈ {tri_at_max_S:.1f}")
    
    return summary


if __name__ == "__main__":
    t0 = time.time()
    results = run()
    print(f"\n  Temps total: {time.time()-t0:.0f}s")
