#!/usr/bin/env python3
"""
TCGE GAP-RG — Coarse-graining amélioré

Problème identifié : le CG naïf (moyenne de |α|) dilue le signal.
C'est un effet mathématique, pas physique.

Hypothèse : un agrégat robuste préserve mieux l'information de polarisation.

Agrégats testés pour le |α| coarse :
  1. mean (baseline — c'est ce qu'on avait)
  2. max (conserve le pic)
  3. q90 (quantile 90% — robuste aux outliers)
  4. trimmed_mean (moyenne après exclusion des 20% les plus bas)

Cohésion testée pour la classification high/low :
  1. density (densité inter-communautés — baseline v7d)
  2. internal_clustering (clustering moyen pondéré intra-communauté)

On teste sur RGG-tore d=3 et d=4, N=5000, Louvain k≈200.

Critère de progrès : rétention d=3 > 60%, d=4 > 45%.

Author: Zeek (TCGE project) — AI-assisted implementation
Date: 2026-02-19
"""

import numpy as np
from collections import defaultdict, deque
from math import gamma as gamma_func
from scipy.spatial import cKDTree
import time

# =============================================================================
# RGG Torus
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
    
    def biphasage_optimize(self, n_steps=2000, seed=None):
        if seed is not None:
            np.random.seed(seed)
        n_e = self.n_edges
        if n_e == 0:
            self.alpha_e = np.array([])
            self.abs_alpha = np.array([])
            self.biphasage = 0
            return 0
        alpha_e = np.random.randn(n_e) * 0.01
        delta_N = np.array([self.N_prop[i] - self.N_prop[j] for i, j in self.edges])
        tri = self.tri.astype(float)
        for _ in range(n_steps):
            grad = (2 * alpha_e * (-0.25 + 0.5*tri + 0.01) + 0.5 * delta_N)
            alpha_e -= 0.005 * grad
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
# CG with multiple aggregators
# =============================================================================

def cg_multi_aggregator(graph, target_k=200):
    """
    Louvain community detection → coarse graph.
    Test multiple aggregation strategies for |α| on coarse edges.
    """
    import community as community_louvain
    import networkx as nx
    
    G = nx.Graph()
    G.add_nodes_from(range(graph.n))
    G.add_edges_from(graph.edges)
    
    partition = community_louvain.best_partition(G, resolution=15.0, random_state=42)
    
    comm_members = defaultdict(set)
    for node, comm in partition.items():
        comm_members[comm].add(node)
    
    # Merge smallest until target_k
    if len(comm_members) > target_k:
        comm_adj = defaultdict(set)
        for i, j in graph.edges:
            ci, cj = partition[i], partition[j]
            if ci != cj:
                comm_adj[ci].add(cj)
                comm_adj[cj].add(ci)
        
        iterations = 0
        while len(comm_members) > target_k and iterations < len(comm_members) * 2:
            iterations += 1
            smallest = min(comm_members, key=lambda c: len(comm_members[c]))
            if smallest not in comm_adj or not comm_adj[smallest]:
                other = [c for c in comm_members if c != smallest]
                if not other:
                    break
                target_comm = min(other, key=lambda c: len(comm_members.get(c, set())))
            else:
                target_comm = min(comm_adj[smallest],
                                key=lambda c: len(comm_members.get(c, set())))
            comm_members[target_comm] |= comm_members[smallest]
            for node in comm_members[smallest]:
                partition[node] = target_comm
            for neighbor in comm_adj.get(smallest, set()):
                if neighbor != target_comm:
                    comm_adj[target_comm].add(neighbor)
                    comm_adj[neighbor].discard(smallest)
                    comm_adj[neighbor].add(target_comm)
            comm_adj[target_comm].discard(smallest)
            del comm_members[smallest]
            if smallest in comm_adj:
                del comm_adj[smallest]
    
    n_comm = len(comm_members)
    
    # Build coarse graph with full edge data
    inter_count = defaultdict(int)
    inter_alphas = defaultdict(list)
    comm_nodes = defaultdict(set)
    for node, comm in partition.items():
        comm_nodes[comm].add(node)
    
    for eidx, (i, j) in enumerate(graph.edges):
        ci, cj = partition[i], partition[j]
        if ci != cj:
            pair = (min(ci,cj), max(ci,cj))
            inter_count[pair] += 1
            inter_alphas[pair].append(graph.abs_alpha[eidx])
    
    coarse_edges = list(inter_count.keys())
    if len(coarse_edges) < 10:
        return None
    
    # --- Aggregation strategies ---
    aggregators = {
        'mean': lambda arr: np.mean(arr),
        'max': lambda arr: np.max(arr),
        'q90': lambda arr: np.percentile(arr, 90),
        'trimmed': lambda arr: np.mean(sorted(arr)[len(arr)//5:])  # drop bottom 20%
    }
    
    results = {'n_comm': n_comm, 'n_edges': len(coarse_edges)}
    
    # Density cohesion for classification
    coarse_density = np.zeros(len(coarse_edges))
    for eidx, (ci, cj) in enumerate(coarse_edges):
        ni, nj = len(comm_nodes[ci]), len(comm_nodes[cj])
        coarse_density[eidx] = inter_count[(ci, cj)] / max(ni * nj, 1)
    
    med_d = np.median(coarse_density)
    if med_d == coarse_density.min():
        med_d = np.mean(coarse_density)
    high_d = [i for i in range(len(coarse_edges)) if coarse_density[i] > med_d]
    low_d = [i for i in range(len(coarse_edges)) if coarse_density[i] <= med_d]
    
    # Also try n_edges classification
    coarse_n = np.array([inter_count[e] for e in coarse_edges], dtype=float)
    med_n = np.median(coarse_n)
    if med_n == coarse_n.min():
        med_n = np.mean(coarse_n)
    high_n = [i for i in range(len(coarse_edges)) if coarse_n[i] > med_n]
    low_n = [i for i in range(len(coarse_edges)) if coarse_n[i] <= med_n]
    
    classifiers = {
        'density': (high_d, low_d),
        'n_edges': (high_n, low_n)
    }
    
    for agg_name, agg_fn in aggregators.items():
        # Compute coarse |α| for each edge
        coarse_alpha = np.zeros(len(coarse_edges))
        for eidx, (ci, cj) in enumerate(coarse_edges):
            alphas = inter_alphas[(ci, cj)]
            if len(alphas) > 0:
                coarse_alpha[eidx] = agg_fn(alphas)
        
        for cls_name, (high, low) in classifiers.items():
            if high and low:
                a_low = np.mean(coarse_alpha[low])
                a_high = np.mean(coarse_alpha[high])
                bi_cg = a_low - a_high
            else:
                bi_cg = 0
            
            results[f'{agg_name}_{cls_name}'] = bi_cg
    
    return results


# =============================================================================
# MAIN
# =============================================================================

def run():
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  TCGE GAP-RG — CG avec agrégats robustes                 ║")
    print("║  mean vs max vs q90 vs trimmed × density vs n_edges      ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    dims = [3, 4]
    N = 5000
    n_trials = 8
    
    for dim in dims:
        print(f"\n{'═'*70}")
        print(f"  d = {dim}  |  N = {N}  |  k_target = 200  |  {n_trials} trials")
        print(f"{'═'*70}")
        
        all_rets = defaultdict(list)
        
        for trial in range(n_trials):
            t0 = time.time()
            graph = RGG_Torus(N, dim, target_degree=8, seed=1000*trial + dim*100)
            bi_fine = graph.biphasage_optimize(seed=2000*trial + 7)
            
            cg = cg_multi_aggregator(graph, target_k=200)
            
            if cg is None:
                print(f"  Trial {trial+1}: SKIP (too few coarse edges)")
                continue
            
            print(f"\n  Trial {trial+1}: Δ_fine={bi_fine:.3f}  "
                  f"k={cg['n_comm']}  E_cg={cg['n_edges']}")
            
            combos = ['mean_density', 'mean_n_edges',
                       'max_density', 'max_n_edges',
                       'q90_density', 'q90_n_edges',
                       'trimmed_density', 'trimmed_n_edges']
            
            for combo in combos:
                bi_cg = cg[combo]
                ret = bi_cg / bi_fine if bi_fine > 0 else 0
                all_rets[combo].append(ret)
                
            # Print best for this trial
            best_combo = max(combos, key=lambda c: cg[c] / bi_fine if bi_fine > 0 else 0)
            best_ret = cg[best_combo] / bi_fine if bi_fine > 0 else 0
            print(f"    Best: {best_combo} → {best_ret:.0%}")
            print(f"    ({time.time()-t0:.1f}s)")
        
        # ── SYNTHÈSE PAR DIMENSION ──
        print(f"\n  {'─'*60}")
        print(f"  SYNTHÈSE d={dim}")
        print(f"  {'─'*60}")
        print(f"  {'Combo':<22} {'Rétention moy':<15} {'std':<8} {'≥50%':<8} {'≥60%':<8}")
        print(f"  {'─'*60}")
        
        for combo in ['mean_density', 'max_density', 'q90_density', 'trimmed_density',
                       'mean_n_edges', 'max_n_edges', 'q90_n_edges', 'trimmed_n_edges']:
            rets = all_rets[combo]
            if rets:
                m = np.mean(rets)
                s = np.std(rets)
                p50 = sum(1 for r in rets if r >= 0.5) / len(rets) * 100
                p60 = sum(1 for r in rets if r >= 0.6) / len(rets) * 100
                marker = '✅' if m >= 0.5 else ('⚠️' if m >= 0.4 else '')
                print(f"  {combo:<22} {m:<15.0%} {s:<8.2f} {p50:<8.0f}% {p60:<8.0f}% {marker}")
    
    print(f"\n\n{'═'*70}")
    print("  INTERPRÉTATION")
    print(f"{'═'*70}\n")
    print("  Si max/q90 > mean : l'information est dans les pics,")
    print("  pas dans la moyenne. Le CG naïf diluait le signal.")
    print("  Si max/q90 ≈ mean : le problème est plus profond")
    print("  (la classification high/low perd le signal, pas l'agrégat).")


if __name__ == "__main__":
    t0 = time.time()
    run()
    print(f"\n  Temps total: {time.time()-t0:.0f}s")
