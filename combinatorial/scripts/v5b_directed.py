#!/usr/bin/env python3
"""
TCGE — Combinatorial Reformulation v5b: Slice-Based Separation
================================================================
Building on v5's breakthrough (20/20 directional, 6/20 Lorentzian).

Problem identified in v5:
  Spatial neighbors often differ in rank → t(i,j) ≠ 0 → s² < 0
  even for genuinely spatial pairs. The global rank "bleeds" into
  spatial directions.

Fix (conceptually clean):
  1. Use per-edge Δ_v3 = directed[i,j] - directed[j,i] ∈ {-1,0,+1}
  2. Classify each edge as:
     - TEMPORAL: |Δ_v3| = 1 (one direction constrains the other)
     - SPATIAL:  |Δ_v3| = 0 (no directed asymmetry on this edge)
  3. Build spatial graph G_spat from ONLY Δ=0 edges
  4. Build temporal graph G_temp from ONLY |Δ|=1 edges
  5. Define:
     - x(i,j) = distance in G_spat (spatial separation)
     - t(i,j) = distance in G_temp (temporal separation)
  6. s²(i,j) = x² - λt²

This is cleaner than rank-based t: the spatial and temporal
distances are computed on DISJOINT subgraphs, so there's no
contamination by construction.

All quantities remain integers. No real numbers.
"""

import numpy as np
from collections import defaultdict, deque
from itertools import combinations
import json
import sys

sys.setrecursionlimit(50000)


def create_directed_universe(nx=3, ny=3, nz=3, nt=5,
                              temporal_dir_rate=0.7,
                              spatial_dir_rate=0.05,
                              seed=42):
    """Same as v5: lattice + boolean directed constraints."""
    rng = np.random.RandomState(seed)
    
    points = []
    for t in range(nt):
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    points.append((t, x, y, z))
    
    n = len(points)
    coord_to_idx = {pt: i for i, pt in enumerate(points)}
    
    compatible = np.zeros((n, n), dtype=bool)
    edge_direction = {}
    
    directions = [(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)]
    dir_labels = ['t', 'x', 'y', 'z']
    
    for i, pt in enumerate(points):
        for d_idx, delta in enumerate(directions):
            nb = tuple(pt[k] + delta[k] for k in range(4))
            j = coord_to_idx.get(nb)
            if j is not None:
                compatible[i, j] = True
                compatible[j, i] = True
                edge_direction[(i, j)] = dir_labels[d_idx]
                edge_direction[(j, i)] = dir_labels[d_idx]
    
    directed = np.zeros((n, n), dtype=bool)
    
    for i, pt in enumerate(points):
        for d_idx, delta in enumerate(directions):
            nb = tuple(pt[k] + delta[k] for k in range(4))
            j = coord_to_idx.get(nb)
            if j is None:
                continue
            d = dir_labels[d_idx]
            rate = temporal_dir_rate if d == 't' else spatial_dir_rate
            if rng.random() < rate:
                directed[i, j] = True
    
    return points, compatible, directed, edge_direction, n


def compute_per_edge_delta(directed, compatible, n):
    """Δ_v3(i,j) = directed[i,j] - directed[j,i] ∈ {-1, 0, +1}"""
    delta = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if compatible[i, j]:
                delta[i, j] = int(directed[i, j]) - int(directed[j, i])
    return delta


def build_split_graphs(compatible, delta, n):
    """
    Split compatibility graph into spatial and temporal subgraphs
    based on per-edge Δ.
    
    Spatial: edges where |Δ(i,j)| = 0 (no directed asymmetry)
    Temporal: edges where |Δ(i,j)| = 1 (one-way constraint)
    """
    g_spatial = np.zeros((n, n), dtype=bool)
    g_temporal = np.zeros((n, n), dtype=bool)
    
    for i in range(n):
        for j in range(n):
            if compatible[i, j]:
                if delta[i, j] == 0:
                    g_spatial[i, j] = True
                else:
                    g_temporal[i, j] = True
    
    return g_spatial, g_temporal


def bfs_distance(adj, n):
    """Shortest path distance on a boolean adjacency matrix. Integer."""
    dist = np.full((n, n), n + 1, dtype=int)
    np.fill_diagonal(dist, 0)
    for source in range(n):
        visited = {source}
        queue = deque([(source, 0)])
        while queue:
            node, d = queue.popleft()
            for nb in range(n):
                if adj[node, nb] and nb not in visited:
                    visited.add(nb)
                    dist[source, nb] = d + 1
                    queue.append((nb, d + 1))
    return dist


def classify_and_test(points, compatible, edge_direction, 
                       x_dist, t_dist, n, lam):
    """
    Compute s² = x² - λt² and classify per lattice direction.
    """
    INF = n + 1
    by_dir = {'t': [], 'x': [], 'y': [], 'z': []}
    
    for i in range(n):
        for j in range(i + 1, n):
            if not compatible[i, j]:
                continue
            d = edge_direction.get((i, j), '?')
            if d not in by_dir:
                continue
            
            xi = x_dist[i, j]
            ti = t_dist[i, j]
            
            # Skip if disconnected in either subgraph
            if xi >= INF and ti >= INF:
                continue
            
            # If disconnected in spatial graph → purely temporal
            if xi >= INF:
                xi = 0
            # If disconnected in temporal graph → purely spatial
            if ti >= INF:
                ti = 0
            
            s2 = int(xi)**2 - lam * int(ti)**2
            by_dir[d].append(s2)
    
    results = {}
    for d in ['t', 'x', 'y', 'z']:
        vals = np.array(by_dir[d]) if by_dir[d] else np.array([0])
        results[d] = {
            "n": len(vals),
            "mean": round(float(np.mean(vals)), 2),
            "neg_frac": round(float(np.mean(vals < 0)), 3),
            "pos_frac": round(float(np.mean(vals > 0)), 3),
            "zero_frac": round(float(np.mean(vals == 0)), 3),
        }
    
    t_neg = results['t']['neg_frac']
    s_pos = np.mean([results[d]['pos_frac'] for d in ['x', 'y', 'z']])
    results['t_neg'] = round(float(t_neg), 3)
    results['s_pos'] = round(float(s_pos), 3)
    results['lorentzian'] = bool(t_neg > 0.4 and s_pos > 0.4)
    results['score'] = round(float(t_neg + s_pos - 1.0), 3)
    
    return results


def run_v5b():
    print("=" * 72)
    print("  TCGE — v5b: SLICE-BASED SPATIAL/TEMPORAL SEPARATION")
    print("  Per-edge Δ classifies edges → disjoint subgraphs")
    print("  x = distance in G_spatial, t = distance in G_temporal")
    print("=" * 72)
    
    # ==============================================================
    # DETAILED RUN
    # ==============================================================
    points, compatible, directed, edge_dir, n = \
        create_directed_universe(seed=42)
    
    delta = compute_per_edge_delta(directed, compatible, n)
    g_spat, g_temp = build_split_graphs(compatible, delta, n)
    
    n_spat = int(g_spat.sum()) // 2
    n_temp = int(g_temp.sum()) // 2
    n_total = int(compatible.sum()) // 2
    
    print(f"\n  Lattice: {n} atoms, {n_total} compatible pairs")
    print(f"  Edge classification by Δ_v3:")
    print(f"    Spatial  (|Δ|=0): {n_spat} edges ({100*n_spat/n_total:.0f}%)")
    print(f"    Temporal (|Δ|=1): {n_temp} edges ({100*n_temp/n_total:.0f}%)")
    
    # Check: how does classification align with lattice direction?
    class_by_dir = {'t': {'spat': 0, 'temp': 0}, 
                     'x': {'spat': 0, 'temp': 0},
                     'y': {'spat': 0, 'temp': 0}, 
                     'z': {'spat': 0, 'temp': 0}}
    
    for i in range(n):
        for j in range(i+1, n):
            if not compatible[i, j]:
                continue
            d = edge_dir.get((i, j), '?')
            if d not in class_by_dir:
                continue
            if g_spat[i, j]:
                class_by_dir[d]['spat'] += 1
            elif g_temp[i, j]:
                class_by_dir[d]['temp'] += 1
    
    print(f"\n  Edge classification vs lattice direction:")
    print(f"  {'dir':>4} {'spatial':>8} {'temporal':>8} {'%temporal':>10}")
    print(f"  {'─'*4} {'─'*8} {'─'*8} {'─'*10}")
    for d in ['t', 'x', 'y', 'z']:
        s = class_by_dir[d]['spat']
        t = class_by_dir[d]['temp']
        pct = 100 * t / max(1, s + t)
        marker = " ← correctly temporal" if d == 't' and pct > 50 else ""
        print(f"  {d:>4} {s:>8} {t:>8} {pct:>9.0f}%{marker}")
    
    # Compute distances on split graphs
    print(f"\n  Computing spatial distances (G_spatial)...")
    x_dist = bfs_distance(g_spat, n)
    print(f"  Computing temporal distances (G_temporal)...")
    t_dist = bfs_distance(g_temp, n)
    
    # Statistics
    spatial_edges_x = []
    spatial_edges_t = []
    temporal_edges_x = []
    temporal_edges_t = []
    INF = n + 1
    
    for i in range(n):
        for j in range(i+1, n):
            if not compatible[i, j]:
                continue
            d = edge_dir.get((i, j), '?')
            xi = x_dist[i, j] if x_dist[i, j] < INF else -1
            ti = t_dist[i, j] if t_dist[i, j] < INF else -1
            
            if d == 't':
                if xi >= 0: temporal_edges_x.append(xi)
                if ti >= 0: temporal_edges_t.append(ti)
            elif d in ['x', 'y', 'z']:
                if xi >= 0: spatial_edges_x.append(xi)
                if ti >= 0: spatial_edges_t.append(ti)
    
    print(f"\n  Distance statistics:")
    print(f"    Temporal edges: mean x_spat={np.mean(temporal_edges_x):.2f}, "
          f"mean t_temp={np.mean(temporal_edges_t):.2f}")
    print(f"    Spatial edges:  mean x_spat={np.mean(spatial_edges_x):.2f}, "
          f"mean t_temp={np.mean(spatial_edges_t):.2f}")
    
    # The KEY prediction:
    # Temporal edges: high t, low x → s² < 0 (timelike)
    # Spatial edges:  low t, high x → s² > 0 (spacelike)
    # But this only works if the split graphs separate directions correctly.
    
    # Lambda scan
    print(f"\n  λ calibration (s² = x_spat² − λ·t_temp²):")
    print(f"  {'λ':>4} {'t_neg':>8} {'s_pos':>8} {'score':>8}")
    print(f"  {'─'*4} {'─'*8} {'─'*8} {'─'*8}")
    
    best_lam = 1
    best_score = -999
    
    for lam in range(1, 20):
        sig = classify_and_test(points, compatible, edge_dir,
                                x_dist, t_dist, n, lam)
        score = sig['score']
        marker = ""
        if score > best_score:
            best_score = score
            best_lam = lam
            marker = " ←"
        print(f"  {lam:>4} {sig['t_neg']:>7.1%} {sig['s_pos']:>7.1%} "
              f"{score:>+7.3f}{marker}")
    
    # Detailed results at best λ
    sig_best = classify_and_test(points, compatible, edge_dir,
                                  x_dist, t_dist, n, best_lam)
    
    print(f"\n  SIGNATURE at λ={best_lam}:")
    for d in ['t', 'x', 'y', 'z']:
        r = sig_best[d]
        print(f"    {d}: mean s²={r['mean']:>8.1f}, "
              f"neg={r['neg_frac']:.1%}, pos={r['pos_frac']:.1%}")
    
    lor = sig_best['lorentzian']
    print(f"\n    LORENTZIAN: {'✅ YES' if lor else '❌ NO'}")
    if lor:
        print(f"    🎯 Temporal={sig_best['t_neg']:.0%} timelike, "
              f"Spatial={sig_best['s_pos']:.0%} spacelike")
    
    # ==============================================================
    # ROBUSTNESS: 20 seeds
    # ==============================================================
    print(f"\n" + "=" * 72)
    print(f"  ROBUSTNESS: 20 seeds")
    print(f"=" * 72)
    
    lor_count = 0
    dir_count = 0
    fol_count = 0
    scores = []
    
    for seed in range(20):
        pts, comp, dire, edir, nn = create_directed_universe(seed=seed)
        delt = compute_per_edge_delta(dire, comp, nn)
        gs, gt = build_split_graphs(comp, delt, nn)
        
        xd = bfs_distance(gs, nn)
        td = bfs_distance(gt, nn)
        
        # Check directionality
        by_d = {'t': [], 'x': [], 'y': [], 'z': []}
        for i in range(nn):
            for j in range(i+1, nn):
                if comp[i, j]:
                    d = edir.get((i, j), '?')
                    if d in by_d:
                        by_d[d].append(abs(delt[i, j]))
        t_delta = np.mean(by_d['t']) if by_d['t'] else 0
        s_delta = np.mean([np.mean(by_d[d]) for d in ['x','y','z'] if by_d[d]])
        directional = t_delta > s_delta * 1.2
        if directional:
            dir_count += 1
        
        # Find best λ and test signature
        best_l = 1
        best_s = -999
        for lam in range(1, 15):
            sig = classify_and_test(pts, comp, edir, xd, td, nn, lam)
            if sig['score'] > best_s:
                best_s = sig['score']
                best_l = lam
        
        sig = classify_and_test(pts, comp, edir, xd, td, nn, best_l)
        scores.append(best_s)
        
        if sig['lorentzian']:
            lor_count += 1
        
        # Foliation from directed graph rank
        from collections import deque as dq
        adj_d = defaultdict(list)
        in_d = np.zeros(nn, dtype=int)
        for i in range(nn):
            for j in range(nn):
                if dire[i, j]:
                    adj_d[i].append(j)
                    in_d[j] += 1
        rank = np.zeros(nn, dtype=int)
        queue = deque([i for i in range(nn) if in_d[i] == 0])
        while queue:
            node = queue.popleft()
            for nb in adj_d[node]:
                rank[nb] = max(rank[nb], rank[node] + 1)
                in_d[nb] -= 1
                if in_d[nb] == 0:
                    queue.append(nb)
        if len(set(rank)) > 3:
            fol_count += 1
        
        sd = "✅" if directional else "❌"
        sl = "✅" if sig['lorentzian'] else "❌"
        print(f"    seed={seed:>2}: Δ={sd} | λ={best_l:>2} "
              f"t_neg={sig['t_neg']:.0%} s_pos={sig['s_pos']:.0%} "
              f"score={best_s:+.3f} {sl}")
    
    print(f"\n  SUMMARY:")
    print(f"    Δ directional:        {dir_count}/20 ({100*dir_count/20:.0f}%)")
    print(f"    Lorentzian signature: {lor_count}/20 ({100*lor_count/20:.0f}%)")
    print(f"    Foliated structure:   {fol_count}/20 ({100*fol_count/20:.0f}%)")
    print(f"    Mean best score:      {np.mean(scores):+.3f}")
    
    # ==============================================================
    # ISOTROPIC CONTROL
    # ==============================================================
    print(f"\n" + "=" * 72)
    print(f"  CONTROL: ISOTROPIC")
    print(f"=" * 72)
    
    iso_lor = 0
    for seed in range(10):
        pts, comp, dire, edir, nn = create_directed_universe(
            nx=3, ny=3, nz=3, nt=5,
            temporal_dir_rate=0.3, spatial_dir_rate=0.3, seed=seed)
        delt = compute_per_edge_delta(dire, comp, nn)
        gs, gt = build_split_graphs(comp, delt, nn)
        xd = bfs_distance(gs, nn)
        td = bfs_distance(gt, nn)
        
        best_l, best_s = 1, -999
        for lam in range(1, 15):
            sig = classify_and_test(pts, comp, edir, xd, td, nn, lam)
            if sig['score'] > best_s:
                best_s = sig['score']
                best_l = lam
        
        sig = classify_and_test(pts, comp, edir, xd, td, nn, best_l)
        if sig['lorentzian']:
            iso_lor += 1
        print(f"    seed={seed}: score={best_s:+.3f} "
              f"{'✅' if sig['lorentzian'] else '❌'}")
    
    print(f"\n  Isotropic Lorentzian: {iso_lor}/10")
    print(f"  → {'Good: no spurious signature' if iso_lor <= 2 else '⚠️ Some spurious'}")
    
    # ==============================================================
    # VARYING ANISOTROPY
    # ==============================================================
    print(f"\n" + "=" * 72)
    print(f"  PHASE DIAGRAM: varying temporal constraint rate")
    print(f"=" * 72)
    
    print(f"\n  {'t_rate':>7} {'n_temp':>7} {'n_spat':>7} "
          f"{'best_λ':>7} {'t_neg':>7} {'s_pos':>7} {'score':>7} {'Lor':>5}")
    print(f"  {'─'*7} {'─'*7} {'─'*7} {'─'*7} {'─'*7} {'─'*7} {'─'*7} {'─'*5}")
    
    for t_rate in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        pts, comp, dire, edir, nn = create_directed_universe(
            temporal_dir_rate=t_rate, spatial_dir_rate=0.05, seed=42)
        delt = compute_per_edge_delta(dire, comp, nn)
        gs, gt = build_split_graphs(comp, delt, nn)
        
        nt_edges = int(gt.sum()) // 2
        ns_edges = int(gs.sum()) // 2
        
        xd = bfs_distance(gs, nn)
        td = bfs_distance(gt, nn)
        
        best_l, best_s = 1, -999
        for lam in range(1, 20):
            sig = classify_and_test(pts, comp, edir, xd, td, nn, lam)
            if sig['score'] > best_s:
                best_s = sig['score']
                best_l = lam
        
        sig = classify_and_test(pts, comp, edir, xd, td, nn, best_l)
        lor = "✅" if sig['lorentzian'] else "❌"
        
        print(f"  {t_rate:>7.1f} {nt_edges:>7} {ns_edges:>7} "
              f"{best_l:>7} {sig['t_neg']:>6.0%} {sig['s_pos']:>6.0%} "
              f"{best_s:>+6.3f} {lor:>5}")
    
    # ==============================================================
    # CONCLUSION
    # ==============================================================
    print(f"\n" + "=" * 72)
    print(f"  FINAL ASSESSMENT")
    print(f"=" * 72)
    print()
    print(f"  v5b vs v5 improvement:")
    print(f"    v5:  6/20 Lorentzian (30%)")
    print(f"    v5b: {lor_count}/20 Lorentzian ({100*lor_count/20:.0f}%)")
    print()
    
    if lor_count > 6:
        print(f"  ✅ IMPROVEMENT: Slice-based separation works better")
    elif lor_count >= 6:
        print(f"  ≈  COMPARABLE to v5 — separation method alone is not enough")
    else:
        print(f"  ❌ REGRESSION — slice method loses information")
    
    print()
    print(f"  ALL DATA TYPES: boolean (directed[i,j]) and integer (distances)")
    print(f"  ZERO real numbers in the substrate")

    # Save
    output = {
        "robustness": {
            "directional": dir_count,
            "lorentzian": lor_count, 
            "foliation": fol_count,
            "mean_score": round(float(np.mean(scores)), 3)
        }
    }
    with open("tcge_directed_v5b_results.json", "w") as f:
        json.dump(output, f, indent=2)


if __name__ == "__main__":
    run_v5b()
