#!/usr/bin/env python3
"""
TCGE — Contrôles causaux pour le rewiring inverse

Contrôle 1 : rewiring aléatoire degree-preserving (pas de critère triangle)
  → Si S reste ≈ 0 : confirme que c'est bien les triangles, pas le rewiring.

Contrôle 2 : triangle-closing MAIS protecteur désactivé (C_protect = 0)
  → Si S reste ≈ 0 : confirme que triangles seuls ne suffisent pas —
    il faut le terme protecteur dans le coût.

Comparaison avec le résultat principal (triangle-closing + protecteur actif).

Author: Zeek (TCGE project) — AI-assisted implementation
Date: 2026-02-20
"""

import numpy as np
from collections import defaultdict
from scipy.special import betaln
import time


# =============================================================================
# BETA MIXTURE (compact)
# =============================================================================

def beta_pdf(x, a, b):
    x = np.clip(x, 1e-10, 1 - 1e-10)
    return np.exp((a-1)*np.log(x) + (b-1)*np.log(1-x) - betaln(a, b))

def fit_beta_mixture(data, max_iter=100, tol=1e-6):
    data = np.clip(data, 1e-6, 1 - 1e-6)
    med = np.median(data)
    mask_lo, mask_hi = data <= med, data > med
    def mom(x):
        m, v = np.mean(x), np.var(x)
        if v >= m*(1-m): v = m*(1-m)*0.9
        c = m*(1-m)/v - 1
        return max(m*c, 0.5), max((1-m)*c, 0.5)
    if sum(mask_lo)>2 and sum(mask_hi)>2:
        a1,b1 = mom(data[mask_lo]); a2,b2 = mom(data[mask_hi])
    else:
        a1,b1,a2,b2 = 2,5,5,2
    pi1 = 0.5
    for _ in range(max_iter):
        p1, p2 = beta_pdf(data,a1,b1), beta_pdf(data,a2,b2)
        d = pi1*p1 + (1-pi1)*p2 + 1e-300
        gamma = pi1*p1/d
        pi1_n = np.clip(np.mean(gamma), 0.01, 0.99)
        def wb(g):
            w = g/(g.sum()+1e-300); m = np.sum(w*data)
            v = np.sum(w*(data-m)**2); m = np.clip(m,0.01,0.99)
            if v >= m*(1-m): v = m*(1-m)*0.9
            if v < 1e-10: v = 1e-4
            c = m*(1-m)/v - 1
            return max(m*c,0.1), max((1-m)*c,0.1)
        a1,b1 = wb(gamma); a2,b2 = wb(1-gamma)
        if abs(pi1_n-pi1) < tol: break
        pi1 = pi1_n
    mu1,mu2 = a1/(a1+b1), a2/(a2+b2)
    v1 = a1*b1/((a1+b1)**2*(a1+b1+1)); v2 = a2*b2/((a2+b2)**2*(a2+b2+1))
    if mu1>mu2: mu1,mu2=mu2,mu1; v1,v2=v2,v1; pi1=1-pi1
    return mu1,mu2,np.sqrt(v1),np.sqrt(v2),pi1

def S_score(data):
    mu1,mu2,s1,s2,pi1 = fit_beta_mixture(data)
    dist = abs(mu2-mu1)
    return dist * 2*min(pi1,1-pi1), mu1, mu2, pi1


# =============================================================================
# GRAPH + REWIRING
# =============================================================================

def make_ER(N, target_k, seed):
    np.random.seed(seed)
    p = target_k/(N-1)
    edges, adj = [], defaultdict(set)
    for i in range(N):
        for j in range(i+1, N):
            if np.random.random() < p:
                edges.append((i,j)); adj[i].add(j); adj[j].add(i)
    return edges, adj

def random_rewire(edges, adj, N, n_swaps, seed):
    """Standard degree-preserving rewiring (no triangle preference)."""
    np.random.seed(seed)
    edges = list(edges); edge_set = set(edges)
    adj = defaultdict(set)
    for i,j in edges: adj[i].add(j); adj[j].add(i)
    n_e = len(edges); done = 0; att = 0
    while done < n_swaps and att < n_swaps*10:
        att += 1
        i1,i2 = np.random.randint(n_e), np.random.randint(n_e)
        if i1==i2: continue
        a,b = edges[i1]; c,d = edges[i2]
        if np.random.random()<0.5:
            n1,n2 = (min(a,d),max(a,d)), (min(c,b),max(c,b))
        else:
            n1,n2 = (min(a,c),max(a,c)), (min(b,d),max(b,d))
        if n1[0]==n1[1] or n2[0]==n2[1]: continue
        if n1 in edge_set or n2 in edge_set or n1==n2: continue
        edge_set.discard(edges[i1]); edge_set.discard(edges[i2])
        adj[a].discard(b); adj[b].discard(a)
        adj[c].discard(d); adj[d].discard(c)
        edges[i1]=n1; edges[i2]=n2
        edge_set.add(n1); edge_set.add(n2)
        adj[n1[0]].add(n1[1]); adj[n1[1]].add(n1[0])
        adj[n2[0]].add(n2[1]); adj[n2[1]].add(n2[0])
        done += 1
    return edges, adj

def triangle_closing_rewire(edges, adj, N, n_attempts, seed):
    """Only accepts swaps that increase triangles."""
    np.random.seed(seed)
    edges = list(edges); edge_set = set(edges)
    adj = defaultdict(set)
    for i,j in edges: adj[i].add(j); adj[j].add(i)
    n_e = len(edges)
    for _ in range(n_attempts):
        i1,i2 = np.random.randint(n_e), np.random.randint(n_e)
        if i1==i2: continue
        a,b = edges[i1]; c,d = edges[i2]
        tri_before = len(adj[a]&adj[b]) + len(adj[c]&adj[d])
        if np.random.random()<0.5:
            n1,n2 = (min(a,d),max(a,d)), (min(c,b),max(c,b))
        else:
            n1,n2 = (min(a,c),max(a,c)), (min(b,d),max(b,d))
        if n1[0]==n1[1] or n2[0]==n2[1]: continue
        if n1 in edge_set or n2 in edge_set or n1==n2: continue
        adj[a].discard(b); adj[b].discard(a)
        adj[c].discard(d); adj[d].discard(c)
        adj[n1[0]].add(n1[1]); adj[n1[1]].add(n1[0])
        adj[n2[0]].add(n2[1]); adj[n2[1]].add(n2[0])
        tri_after = len(adj[n1[0]]&adj[n1[1]]) + len(adj[n2[0]]&adj[n2[1]])
        if tri_after > tri_before:
            edge_set.discard(edges[i1]); edge_set.discard(edges[i2])
            edges[i1]=n1; edges[i2]=n2
            edge_set.add(n1); edge_set.add(n2)
        else:
            adj[n1[0]].discard(n1[1]); adj[n1[1]].discard(n1[0])
            adj[n2[0]].discard(n2[1]); adj[n2[1]].discard(n2[0])
            adj[a].add(b); adj[b].add(a)
            adj[c].add(d); adj[d].add(c)
    return edges, adj


# =============================================================================
# BIPHASAGE (with optional C_protect)
# =============================================================================

def run_biphasage(N, edges, adj, C_protect=0.5, seed=None):
    if seed is not None: np.random.seed(seed)
    n_e = len(edges)
    if n_e < 10:
        return np.array([]), 0, 0
    degree = np.zeros(N)
    for i,j in edges: degree[i]+=1; degree[j]+=1
    tri = np.zeros(n_e, dtype=int)
    for idx,(i,j) in enumerate(edges):
        tri[idx] = len(adj[i] & adj[j])
    
    alpha_e = np.random.randn(n_e) * 0.01
    delta_N = np.array([degree[i]-degree[j] for i,j in edges])
    tri_f = tri.astype(float)
    
    for _ in range(2500):
        grad = (2*alpha_e*(-0.25 + C_protect*tri_f + 0.01) + 0.5*delta_N)
        alpha_e -= 0.005 * grad
        alpha_e = np.clip(alpha_e, -0.999, 0.999)
    
    return np.abs(alpha_e), np.mean(tri), np.mean(degree)


# =============================================================================
# MAIN
# =============================================================================

def run():
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  Contrôles causaux                                        ║")
    print("║  C1: rewire aléatoire (pas de tri)                        ║")
    print("║  C2: tri-closing + protecteur OFF (C_protect=0)           ║")
    print("║  Réf: tri-closing + protecteur ON  (C_protect=0.5)        ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    N = 3000
    target_k = 8
    n_trials = 5
    n_attempts = 500000  # same budget for all conditions
    
    conditions = {
        'ER (baseline)': {
            'rewire': None, 'C_protect': 0.5
        },
        'C1: random rewire': {
            'rewire': 'random', 'C_protect': 0.5
        },
        'C2: tri-close + prot OFF': {
            'rewire': 'triangle', 'C_protect': 0.0
        },
        'REF: tri-close + prot ON': {
            'rewire': 'triangle', 'C_protect': 0.5
        },
    }
    
    results = {}
    
    for name, cfg in conditions.items():
        t0 = time.time()
        all_S, all_mu1, all_mu2, all_pi1 = [], [], [], []
        all_tri, all_alpha = [], []
        
        for trial in range(n_trials):
            edges, adj = make_ER(N, target_k, seed=1000*trial + 42)
            
            if cfg['rewire'] == 'random':
                edges, adj = random_rewire(edges, adj, N, n_attempts, 
                                           seed=5000*trial + 1)
            elif cfg['rewire'] == 'triangle':
                edges, adj = triangle_closing_rewire(edges, adj, N, n_attempts,
                                                      seed=5000*trial + 2)
            
            aa, tri_m, k_m = run_biphasage(N, edges, adj, 
                                            C_protect=cfg['C_protect'],
                                            seed=2000*trial + 7)
            
            if len(aa) > 0:
                S, mu1, mu2, pi1 = S_score(aa)
                all_S.append(S)
                all_mu1.append(mu1)
                all_mu2.append(mu2)
                all_pi1.append(pi1)
                all_alpha.append(np.mean(aa))
                all_tri.append(tri_m)
        
        results[name] = {
            'S': np.mean(all_S), 'S_std': np.std(all_S),
            'mu1': np.mean(all_mu1), 'mu2': np.mean(all_mu2),
            'pi1': np.mean(all_pi1),
            'tri': np.mean(all_tri), 'alpha': np.mean(all_alpha),
        }
        
        r = results[name]
        print(f"\n  {name}")
        print(f"    ⟨tri⟩={r['tri']:.2f}  ⟨|α|⟩={r['alpha']:.3f}  "
              f"S={r['S']:.3f}±{r['S_std']:.3f}  "
              f"μ₁={r['mu1']:.2f}  μ₂={r['mu2']:.2f}  π₁={r['pi1']:.2f}  "
              f"({time.time()-t0:.1f}s)")
    
    # ── TABLE ──
    print(f"\n\n{'═'*75}")
    print(f"  {'Condition':<28} {'⟨tri⟩':<7} {'⟨|α|⟩':<7} {'S':<10} "
          f"{'μ₁':<6} {'μ₂':<6} {'π₁':<6}")
    print(f"  {'─'*73}")
    for name in conditions:
        r = results[name]
        print(f"  {name:<28} {r['tri']:<7.2f} {r['alpha']:<7.3f} "
              f"{r['S']:<10.3f} {r['mu1']:<6.2f} {r['mu2']:<6.2f} "
              f"{r['pi1']:<6.2f}")
    
    # ── VERDICT ──
    print(f"\n\n{'═'*75}")
    print("  VERDICT")
    print(f"{'═'*75}\n")
    
    S_base = results['ER (baseline)']['S']
    S_c1 = results['C1: random rewire']['S']
    S_c2 = results['C2: tri-close + prot OFF']['S']
    S_ref = results['REF: tri-close + prot ON']['S']
    
    c1_pass = S_c1 < S_ref * 0.3
    c2_pass = S_c2 < S_ref * 0.3
    
    print(f"  C1 (random rewire, no tri):     S = {S_c1:.3f}  "
          f"{'✅ S reste bas' if c1_pass else '❌ S augmente'}")
    print(f"  C2 (tri-close, protecteur OFF): S = {S_c2:.3f}  "
          f"{'✅ S reste bas' if c2_pass else '❌ S augmente'}")
    print(f"  REF (tri-close, protecteur ON): S = {S_ref:.3f}")
    
    if c1_pass and c2_pass:
        print(f"\n  ✅ DOUBLE CONTRÔLE PASSÉ")
        print(f"     La séparation T/S nécessite BOTH :")
        print(f"       1. Triangles (cohésion locale)")
        print(f"       2. Protecteur actif dans le coût (C_protect > 0)")
        print(f"     Ni le rewiring seul, ni les triangles seuls ne suffisent.")
    elif c1_pass:
        print(f"\n  ⚠️ C1 passe mais C2 échoue")
        print(f"     Les triangles + coût sans protecteur créent de la séparation.")
    elif c2_pass:
        print(f"\n  ⚠️ C2 passe mais C1 échoue")
        print(f"     Le rewiring aléatoire crée de la séparation inattendue.")
    else:
        print(f"\n  ❌ Les deux contrôles échouent.")


if __name__ == "__main__":
    t0 = time.time()
    run()
    print(f"\n  Temps total: {time.time()-t0:.0f}s")
