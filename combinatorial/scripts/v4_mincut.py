#!/usr/bin/env python3
"""
TCGE — Combinatorial Reformulation v4: Min-Cut Weights
=========================================================
Implements Sorkin-response strategy:

  w(i,j) = λ_edge(i,j) = min edge-cut between i and j
         = max number of edge-disjoint paths (Menger's theorem)

This is:
  - Purely combinatorial (integer-valued)
  - Canonical (standard graph theory, not ad hoc)
  - Physically meaningful: "how many independent compatibility 
    paths support coexistence of i and j"

Test plan:
  1. Metric axioms with min-cut distance
  2. Signature from s² = w_sym² - Δw² with combinatorial quantities
  3. Anisotropic hypergraph: does Lorentzian signature emerge
     when H_hard is structurally different along one direction?
  4. Foliated causal structure
  5. Robustness across seeds and lattice sizes

Key hypothesis to test:
  If H_hard is ANISOTROPIC (more constraints along one direction),
  the min-cut itself becomes direction-dependent, creating asymmetry
  WITHOUT any ad hoc multiplier β.
"""

import numpy as np
from collections import defaultdict, deque
from itertools import combinations
import json
import sys
import time

sys.setrecursionlimit(50000)


# ===================================================================
# GRAPH CONSTRUCTION
# ===================================================================

def create_anisotropic_hypergraph(nx=4, ny=4, nz=4, nt=4,
                                   spatial_constraint_rate=0.05,
                                   temporal_constraint_rate=0.30,
                                   seed=42):
    """
    Create a 4D lattice hypergraph where constraints are ANISOTROPIC:
    temporal transitions have more hard constraints than spatial ones.
    
    This is the physical content: the universe has a direction along
    which coexistence is harder. That direction becomes time.
    
    All data: booleans and integers. Zero real numbers.
    """
    rng = np.random.RandomState(seed)
    
    points = []
    for t in range(nt):
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    points.append((t, x, y, z))
    
    n = len(points)
    coord_to_idx = {pt: i for i, pt in enumerate(points)}
    
    # Adjacency: nearest neighbors only
    adj = np.zeros((n, n), dtype=bool)
    edge_direction = {}  # (i,j) -> 't','x','y','z'
    
    directions = [(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)]
    dir_labels = ['t', 'x', 'y', 'z']
    
    for i, pt in enumerate(points):
        for d_idx, delta in enumerate(directions):
            neighbor = tuple(pt[k] + delta[k] for k in range(4))
            j = coord_to_idx.get(neighbor)
            if j is not None:
                adj[i, j] = True
                adj[j, i] = True
                edge_direction[(i, j)] = dir_labels[d_idx]
                edge_direction[(j, i)] = dir_labels[d_idx]
    
    # Now remove edges (hard constraints) anisotropically
    compatible = adj.copy()
    removed_edges = {'t': 0, 'x': 0, 'y': 0, 'z': 0}
    total_edges = {'t': 0, 'x': 0, 'y': 0, 'z': 0}
    
    for i in range(n):
        for j in range(i+1, n):
            if not adj[i, j]:
                continue
            d = edge_direction.get((i, j), '?')
            total_edges[d] = total_edges.get(d, 0) + 1
            
            if d == 't':
                rate = temporal_constraint_rate
            else:
                rate = spatial_constraint_rate
            
            if rng.random() < rate:
                compatible[i, j] = False
                compatible[j, i] = False
                removed_edges[d] += 1
    
    return points, compatible, adj, edge_direction, n, removed_edges, total_edges


def create_isotropic_hypergraph(nx=4, ny=4, nz=4, nt=4,
                                 constraint_rate=0.15, seed=42):
    """Control: isotropic constraints (same rate in all directions)."""
    return create_anisotropic_hypergraph(
        nx, ny, nz, nt,
        spatial_constraint_rate=constraint_rate,
        temporal_constraint_rate=constraint_rate,
        seed=seed
    )


# ===================================================================
# MIN-CUT COMPUTATION (Menger / max-flow for unit capacity)
# ===================================================================

def min_edge_cut(compatible, n, source, target):
    """
    Compute min edge-cut between source and target.
    = max number of edge-disjoint paths (Menger's theorem).
    Uses max-flow with unit capacities (BFS augmenting paths).
    
    Returns an INTEGER.
    """
    if source == target:
        return 0
    if not compatible[source, target] and not _bfs_connected(compatible, n, source, target):
        return 0
    
    # Build residual graph (unit capacities)
    residual = compatible.astype(int).copy()
    
    max_flow = 0
    while True:
        # BFS to find augmenting path
        path = _bfs_path_residual(residual, n, source, target)
        if path is None:
            break
        
        # Augment along path (capacity 1)
        for k in range(len(path) - 1):
            u, v = path[k], path[k+1]
            residual[u, v] -= 1
            residual[v, u] += 1
        
        max_flow += 1
    
    return max_flow


def _bfs_path_residual(residual, n, source, target):
    """BFS on residual graph."""
    prev = {source: None}
    queue = deque([source])
    while queue:
        node = queue.popleft()
        if node == target:
            path = []
            while node is not None:
                path.append(node)
                node = prev[node]
            return path[::-1]
        for nb in range(n):
            if residual[node, nb] > 0 and nb not in prev:
                prev[nb] = node
                queue.append(nb)
    return None


def _bfs_connected(compatible, n, source, target):
    """Check if source and target are connected."""
    visited = {source}
    queue = deque([source])
    while queue:
        node = queue.popleft()
        if node == target:
            return True
        for nb in range(n):
            if compatible[node, nb] and nb not in visited:
                visited.add(nb)
                queue.append(nb)
    return False


def compute_all_neighbor_mincuts(points, compatible, adj, edge_direction, n):
    """
    Compute min-cut for all adjacent pairs.
    Returns w_sym (symmetric) and categorized by direction.
    """
    mincuts_by_dir = {'t': [], 'x': [], 'y': [], 'z': []}
    mincut_matrix = np.zeros((n, n), dtype=int)
    
    pairs_computed = 0
    for i in range(n):
        for j in range(i+1, n):
            if not adj[i, j]:  # only nearest neighbors
                continue
            
            mc = min_edge_cut(compatible, n, i, j)
            mincut_matrix[i, j] = mc
            mincut_matrix[j, i] = mc
            
            d = edge_direction.get((i, j), '?')
            if d in mincuts_by_dir:
                mincuts_by_dir[d].append(mc)
            
            pairs_computed += 1
    
    return mincut_matrix, mincuts_by_dir, pairs_computed


def compute_degree(compatible, n):
    """N(i) = degree in compatibility graph. Integer."""
    return compatible.sum(axis=1).astype(int)


def compute_graph_distance(compatible, n):
    """Shortest path distance. Integer."""
    dist = np.full((n, n), n + 1, dtype=int)
    np.fill_diagonal(dist, 0)
    
    for source in range(n):
        visited = {source}
        queue = deque([(source, 0)])
        while queue:
            node, d = queue.popleft()
            for nb in range(n):
                if compatible[node, nb] and nb not in visited:
                    visited.add(nb)
                    dist[source, nb] = d + 1
                    queue.append((nb, d + 1))
    
    return dist


# ===================================================================
# SIGNATURE AND METRIC TESTS
# ===================================================================

def test_mincut_signature(points, adj, edge_direction, mincut_matrix, N, n):
    """
    Test: does min-cut asymmetry create Lorentzian-like separation?
    
    For each adjacent pair (i,j):
      w_sym = mincut(i,j)                    (integer, symmetric)
      Δw = N(i) - N(j)                       (integer, asymmetric)
      s² = w_sym² - Δw²                      (integer)
    
    Also test the PURE min-cut version:
      If min-cut is direction-dependent (lower along temporal),
      this itself creates the signature without needing Δw.
    """
    s2_by_dir = {'t': [], 'x': [], 'y': [], 'z': []}
    
    for i in range(n):
        for j in range(i+1, n):
            if not adj[i, j]:
                continue
            
            d = edge_direction.get((i, j), '?')
            if d not in s2_by_dir:
                continue
            
            w_sym = mincut_matrix[i, j]
            delta_w = N[i] - N[j]
            s_squared = w_sym**2 - delta_w**2
            
            s2_by_dir[d].append(int(s_squared))
    
    results = {}
    for d in ['t', 'x', 'y', 'z']:
        vals = np.array(s2_by_dir[d]) if s2_by_dir[d] else np.array([0])
        results[d] = {
            "n": len(vals),
            "mean_s2": round(float(np.mean(vals)), 2),
            "negative_frac": round(float(np.mean(vals < 0)), 3),
            "positive_frac": round(float(np.mean(vals > 0)), 3),
            "zero_frac": round(float(np.mean(vals == 0)), 3),
        }
    
    # Lorentzian check: temporal should have lower/more negative s² than spatial
    t_mean = results['t']['mean_s2']
    s_mean = np.mean([results[d]['mean_s2'] for d in ['x', 'y', 'z']])
    
    results['temporal_lower'] = bool(t_mean < s_mean)
    results['lorentzian_strict'] = bool(
        results['t']['negative_frac'] > 0.3 and 
        np.mean([results[d]['positive_frac'] for d in ['x', 'y', 'z']]) > 0.3
    )
    
    return results


def test_mincut_as_metric(mincut_matrix, adj, n):
    """
    Test: does 1/mincut work as a distance? Or mincut_max - mincut?
    Note: min-cut itself is a CONNECTIVITY measure (higher = closer).
    So "distance" = some decreasing function of min-cut.
    
    We test: d(i,j) = K - mincut(i,j) where K = max(mincut) + 1
    This should satisfy metric axioms on adjacent pairs.
    """
    K = mincut_matrix.max() + 1
    
    # Only test adjacent pairs
    pairs = [(i, j) for i in range(n) for j in range(i+1, n) if adj[i, j]]
    
    positivity = True
    symmetry = True
    
    for i, j in pairs:
        d_ij = K - mincut_matrix[i, j]
        d_ji = K - mincut_matrix[j, i]
        if d_ij < 0:
            positivity = False
        if d_ij != d_ji:
            symmetry = False
    
    return {"positivity": positivity, "symmetry": symmetry, "K": int(K)}


def test_foliation(compatible, N, n):
    """Test foliated causal structure from N(i)."""
    classes = defaultdict(list)
    for i in range(n):
        classes[int(N[i])].append(i)
    
    # Build directed graph
    adj_dir = defaultdict(list)
    ties = 0
    total = 0
    
    for i in range(n):
        for j in range(i+1, n):
            if compatible[i, j]:
                total += 1
                if N[i] > N[j]:
                    adj_dir[i].append(j)
                elif N[j] > N[i]:
                    adj_dir[j].append(i)
                else:
                    ties += 1
    
    # Cycle detection
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {i: WHITE for i in range(n)}
    cycles = 0
    
    def dfs(u):
        nonlocal cycles
        color[u] = GRAY
        for v in adj_dir[u]:
            if color[v] == GRAY:
                cycles += 1
            elif color[v] == WHITE:
                dfs(v)
        color[u] = BLACK
    
    for u in range(n):
        if color[u] == WHITE:
            dfs(u)
    
    return {
        "unique_N": len(classes),
        "ties": ties,
        "tie_rate": round(ties / max(1, total), 3),
        "cycles": cycles,
        "is_dag": cycles == 0
    }


# ===================================================================
# DISCRETE INTERVAL WITH CALIBRATION
# ===================================================================

def test_discrete_interval(points, adj, edge_direction, N, dist, n):
    """
    Test s²(i,j) = d(i,j)² - λ·t(i,j)²
    where:
      d(i,j) = graph distance (integer)
      t(i,j) = |N(i) - N(j)| (integer)
      λ ∈ {1, 2, 3, ...} (integer calibration)
    
    Find λ that best separates temporal from spatial pairs.
    """
    temporal_pairs = []
    spatial_pairs = []
    
    for i in range(n):
        for j in range(i+1, n):
            if not adj[i, j]:
                continue
            
            d_ij = dist[i, j] if dist[i, j] < n + 1 else 0
            t_ij = abs(int(N[i]) - int(N[j]))
            direction = edge_direction.get((i, j), '?')
            
            if direction == 't':
                temporal_pairs.append((d_ij, t_ij))
            elif direction in ['x', 'y', 'z']:
                spatial_pairs.append((d_ij, t_ij))
    
    results = {}
    best_lambda = 1
    best_separation = -999
    
    for lam in range(1, 10):
        temp_s2 = [d**2 - lam * t**2 for d, t in temporal_pairs]
        spat_s2 = [d**2 - lam * t**2 for d, t in spatial_pairs]
        
        if not temp_s2 or not spat_s2:
            continue
        
        temp_neg = np.mean(np.array(temp_s2) < 0)
        spat_pos = np.mean(np.array(spat_s2) > 0)
        separation = temp_neg + spat_pos - 1  # 1.0 = perfect
        
        results[lam] = {
            "temp_neg_frac": round(float(temp_neg), 3),
            "spat_pos_frac": round(float(spat_pos), 3),
            "separation": round(float(separation), 3)
        }
        
        if separation > best_separation:
            best_separation = separation
            best_lambda = lam
    
    return results, best_lambda, best_separation


# ===================================================================
# MAIN
# ===================================================================

def run_mincut_tests():
    print("=" * 72)
    print("  TCGE — MIN-CUT COMBINATORIAL WEIGHTS")
    print("  w(i,j) = edge connectivity = max edge-disjoint paths")
    print("  All quantities: integers. Zero real numbers in substrate.")
    print("=" * 72)
    
    # ---------------------------------------------------------------
    # ANISOTROPIC MODEL
    # ---------------------------------------------------------------
    print("\n" + "=" * 72)
    print("  PART 1: ANISOTROPIC HYPERGRAPH")
    print("  (temporal constraints denser than spatial)")
    print("=" * 72)
    
    t0 = time.time()
    points, compatible, adj, edge_dir, n, removed, total = \
        create_anisotropic_hypergraph(nx=3, ny=3, nz=3, nt=4,
                                       spatial_constraint_rate=0.05,
                                       temporal_constraint_rate=0.40,
                                       seed=42)
    
    print(f"\n  Lattice: 4×3×3×3 = {n} atoms")
    print(f"  Edge removal rates:")
    for d in ['t', 'x', 'y', 'z']:
        rem = removed.get(d, 0)
        tot = total.get(d, 0)
        pct = 100 * rem / tot if tot > 0 else 0
        print(f"    {d}-direction: {rem}/{tot} removed ({pct:.0f}%)")
    
    N = compute_degree(compatible, n)
    dist = compute_graph_distance(compatible, n)
    print(f"  N(i) range: [{N.min()}, {N.max()}]")
    
    # Min-cut computation
    print(f"\n  Computing min-cuts for all adjacent pairs...")
    mincut_mat, mincuts_by_dir, n_pairs = \
        compute_all_neighbor_mincuts(points, compatible, adj, edge_dir, n)
    elapsed = time.time() - t0
    print(f"  Computed {n_pairs} pairs in {elapsed:.1f}s")
    
    print(f"\n  MIN-CUT BY DIRECTION (the key test):")
    print(f"  " + "-" * 55)
    for d in ['t', 'x', 'y', 'z']:
        vals = mincuts_by_dir.get(d, [])
        if vals:
            arr = np.array(vals)
            print(f"    {d}: mean={np.mean(arr):.2f}, median={np.median(arr):.0f}, "
                  f"range=[{arr.min()}, {arr.max()}], n={len(arr)}")
        else:
            print(f"    {d}: no pairs")
    
    t_mean_mc = np.mean(mincuts_by_dir['t']) if mincuts_by_dir['t'] else 0
    s_means_mc = [np.mean(mincuts_by_dir[d]) for d in ['x','y','z'] if mincuts_by_dir[d]]
    s_mean_mc = np.mean(s_means_mc) if s_means_mc else 0
    
    print(f"\n    Temporal mean min-cut: {t_mean_mc:.2f}")
    print(f"    Spatial mean min-cut:  {s_mean_mc:.2f}")
    mc_aniso = t_mean_mc < s_mean_mc
    print(f"    Temporal < Spatial:    {'✅ YES — anisotropy detected!' if mc_aniso else '❌ NO'}")
    
    if mc_aniso:
        print(f"    → Min-cut is NATURALLY LOWER along temporal direction")
        print(f"    → No ad hoc β needed — asymmetry emerges from H_hard structure")
    
    # Signature test
    print(f"\n  SIGNED INTERVAL s² = mincut² − ΔN²:")
    print(f"  " + "-" * 55)
    sig = test_mincut_signature(points, adj, edge_dir, mincut_mat, N, n)
    
    for d in ['t', 'x', 'y', 'z']:
        r = sig[d]
        print(f"    {d}: mean s²={r['mean_s2']:>7.1f}, "
              f"neg={r['negative_frac']:.1%}, pos={r['positive_frac']:.1%}, "
              f"zero={r['zero_frac']:.1%}")
    
    print(f"\n    Temporal systematically lower: {'✅' if sig['temporal_lower'] else '❌'}")
    print(f"    Strict Lorentzian (neg/pos):   {'✅' if sig['lorentzian_strict'] else '❌'}")
    
    # Discrete interval test
    print(f"\n  DISCRETE INTERVAL s² = d² − λ·t² (integer λ calibration):")
    print(f"  " + "-" * 55)
    intervals, best_lam, best_sep = test_discrete_interval(
        points, adj, edge_dir, N, dist, n)
    
    print(f"    {'λ':>4} {'temp s²<0':>10} {'spat s²>0':>10} {'separation':>12}")
    print(f"    {'─'*4} {'─'*10} {'─'*10} {'─'*12}")
    for lam in sorted(intervals.keys()):
        r = intervals[lam]
        marker = " ← best" if lam == best_lam else ""
        print(f"    {lam:>4} {r['temp_neg_frac']:>9.1%} {r['spat_pos_frac']:>9.1%} "
              f"{r['separation']:>+11.3f}{marker}")
    
    print(f"\n    Best λ = {best_lam} (separation = {best_sep:+.3f})")
    
    # Foliation
    print(f"\n  CAUSAL FOLIATION:")
    print(f"  " + "-" * 55)
    fol = test_foliation(compatible, N, n)
    print(f"    N(i) classes: {fol['unique_N']} temporal slices")
    print(f"    Tie rate: {fol['tie_rate']:.1%}")
    print(f"    Cycles: {fol['cycles']}")
    print(f"    DAG: {'✅' if fol['is_dag'] else '❌'}")
    
    # ---------------------------------------------------------------
    # ISOTROPIC CONTROL
    # ---------------------------------------------------------------
    print(f"\n\n" + "=" * 72)
    print(f"  PART 2: ISOTROPIC CONTROL")
    print(f"  (same constraint rate in all directions)")
    print(f"=" * 72)
    
    pts_iso, comp_iso, adj_iso, edir_iso, n_iso, rem_iso, tot_iso = \
        create_isotropic_hypergraph(nx=3, ny=3, nz=3, nt=4,
                                     constraint_rate=0.15, seed=42)
    
    N_iso = compute_degree(comp_iso, n_iso)
    
    print(f"\n  Computing min-cuts (isotropic)...")
    mc_iso, mcdir_iso, _ = compute_all_neighbor_mincuts(
        pts_iso, comp_iso, adj_iso, edir_iso, n_iso)
    
    print(f"  MIN-CUT BY DIRECTION (isotropic):")
    for d in ['t', 'x', 'y', 'z']:
        vals = mcdir_iso.get(d, [])
        if vals:
            print(f"    {d}: mean={np.mean(vals):.2f}")
    
    t_iso = np.mean(mcdir_iso['t']) if mcdir_iso['t'] else 0
    s_iso = np.mean([np.mean(mcdir_iso[d]) for d in ['x','y','z'] if mcdir_iso[d]])
    print(f"\n    Temporal: {t_iso:.2f}, Spatial: {s_iso:.2f}")
    print(f"    Anisotropy: {'❌ NONE (as expected)' if abs(t_iso - s_iso) < 0.5 else '⚠️ unexpected'}")
    
    sig_iso = test_mincut_signature(pts_iso, adj_iso, edir_iso, mc_iso, N_iso, n_iso)
    print(f"\n  Signature (isotropic):")
    for d in ['t', 'x', 'y', 'z']:
        r = sig_iso[d]
        print(f"    {d}: mean s²={r['mean_s2']:>7.1f}")
    print(f"    Temporal lower: {'YES' if sig_iso['temporal_lower'] else 'NO'}")
    print(f"    → {'No Lorentzian distinction — confirms anisotropy is needed' if not sig_iso['temporal_lower'] else 'Surprising distinction even isotropically'}")
    
    # ---------------------------------------------------------------
    # ROBUSTNESS: 10 seeds (anisotropic)
    # ---------------------------------------------------------------
    print(f"\n\n" + "=" * 72)
    print(f"  PART 3: ROBUSTNESS (10 seeds, anisotropic)")
    print(f"=" * 72)
    
    mc_aniso_count = 0
    sig_lower_count = 0
    foliation_count = 0
    
    for seed in range(10):
        pts, comp, adj_t, edir, nn, _, _ = create_anisotropic_hypergraph(
            nx=3, ny=3, nz=3, nt=4,
            spatial_constraint_rate=0.05,
            temporal_constraint_rate=0.40,
            seed=seed)
        N_t = compute_degree(comp, nn)
        mc_t, mcdir_t, _ = compute_all_neighbor_mincuts(pts, comp, adj_t, edir, nn)
        
        # Min-cut anisotropy
        t_mc = np.mean(mcdir_t['t']) if mcdir_t['t'] else 0
        s_mc = np.mean([np.mean(mcdir_t[d]) for d in ['x','y','z'] if mcdir_t[d]])
        if t_mc < s_mc:
            mc_aniso_count += 1
        
        # Signature
        sig_t = test_mincut_signature(pts, adj_t, edir, mc_t, N_t, nn)
        if sig_t['temporal_lower']:
            sig_lower_count += 1
        
        # Foliation
        fol_t = test_foliation(comp, N_t, nn)
        if fol_t['unique_N'] > 3:
            foliation_count += 1
        
        print(f"    seed={seed}: mc_t={t_mc:.1f} mc_s={s_mc:.1f} "
              f"{'✅' if t_mc < s_mc else '❌'} | "
              f"sig={'✅' if sig_t['temporal_lower'] else '❌'} | "
              f"fol={fol_t['unique_N']}cls")
    
    print(f"\n  Min-cut anisotropy:    {mc_aniso_count}/10 ({100*mc_aniso_count/10:.0f}%)")
    print(f"  Temporal sig. lower:   {sig_lower_count}/10 ({100*sig_lower_count/10:.0f}%)")
    print(f"  Foliated structure:    {foliation_count}/10 ({100*foliation_count/10:.0f}%)")
    
    # ---------------------------------------------------------------
    # VARYING ANISOTROPY
    # ---------------------------------------------------------------
    print(f"\n\n" + "=" * 72)
    print(f"  PART 4: ANISOTROPY PHASE DIAGRAM")
    print(f"  (varying temporal constraint rate, spatial fixed at 5%)")
    print(f"=" * 72)
    
    print(f"\n  {'temp_rate':>10} {'mc_temp':>8} {'mc_spat':>8} {'ratio':>8} {'sig_lower':>10}")
    print(f"  {'─'*10} {'─'*8} {'─'*8} {'─'*8} {'─'*10}")
    
    for t_rate in [0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60]:
        pts, comp, adj_t, edir, nn, _, _ = create_anisotropic_hypergraph(
            nx=3, ny=3, nz=3, nt=4,
            spatial_constraint_rate=0.05,
            temporal_constraint_rate=t_rate,
            seed=42)
        N_t = compute_degree(comp, nn)
        mc_t, mcdir_t, _ = compute_all_neighbor_mincuts(pts, comp, adj_t, edir, nn)
        
        t_mc = np.mean(mcdir_t['t']) if mcdir_t['t'] else 0
        s_mcs = [np.mean(mcdir_t[d]) for d in ['x','y','z'] if mcdir_t[d]]
        s_mc = np.mean(s_mcs) if s_mcs else 0
        ratio = t_mc / s_mc if s_mc > 0 else 0
        
        sig_t = test_mincut_signature(pts, adj_t, edir, mc_t, N_t, nn)
        sl = "✅" if sig_t['temporal_lower'] else "❌"
        
        print(f"  {t_rate:>10.2f} {t_mc:>8.2f} {s_mc:>8.2f} {ratio:>8.2f} {sl:>10}")
    
    # ---------------------------------------------------------------
    # CONCLUSION
    # ---------------------------------------------------------------
    print(f"\n\n" + "=" * 72)
    print(f"  CONCLUSION")
    print(f"=" * 72)
    print()
    print(f"  WHAT WE TESTED:")
    print(f"    w(i,j) = min edge-cut = max edge-disjoint paths")
    print(f"    Purely integer, canonical (Menger's theorem)")
    print(f"    No real numbers, no ad hoc parameters")
    print()
    print(f"  KEY FINDING:")
    if mc_aniso_count >= 7:
        print(f"    ✅ Min-cut is NATURALLY ANISOTROPIC when H_hard is anisotropic")
        print(f"       Temporal min-cut < Spatial min-cut in {mc_aniso_count}/10 seeds")
        print(f"       The 'weights' need not be postulated — they emerge from")
        print(f"       the combinatorial structure of the constraint hypergraph.")
    else:
        print(f"    ⚠️  Min-cut anisotropy detected in {mc_aniso_count}/10 seeds")
        print(f"       Effect is present but may need refinement")
    print()
    print(f"  ON SORKIN'S OBJECTION:")
    print(f"    The 'apparatus of real numbers' is not primitive in TCGE.")
    print(f"    Integer-valued min-cuts on the compatibility graph provide")
    print(f"    all the quantitative information needed. Real numbers")
    print(f"    emerge only in the continuum limit (N → ∞), as ratios")
    print(f"    of large integers — exactly as in statistical mechanics.")
    print()
    print(f"  WHAT REMAINS OPEN:")
    print(f"    - Clean Lorentzian sign flip (g_tt < 0) requires the")
    print(f"      anisotropy to be strong enough — this is physics, not math")
    print(f"    - Continuum limit convergence (Gromov-Hausdorff)")
    print(f"    - Full Einstein equations from combinatorial min-cut")
    
    # Save
    output = {
        "anisotropic": {
            "mincut_temporal": round(float(t_mean_mc), 3),
            "mincut_spatial": round(float(s_mean_mc), 3),
            "anisotropy_detected": bool(mc_aniso),
            "signature": sig,
            "foliation": fol,
            "best_lambda": int(best_lam),
        },
        "robustness": {
            "mc_aniso_rate": mc_aniso_count / 10,
            "sig_lower_rate": sig_lower_count / 10,
            "foliation_rate": foliation_count / 10,
        }
    }
    
    with open("tcge_mincut_results.json", "w") as f:
        json.dump(output, f, indent=2, default=lambda x: 
                  int(x) if isinstance(x, (np.integer,)) else
                  float(x) if isinstance(x, (np.floating,)) else
                  bool(x) if isinstance(x, (np.bool_,)) else x)
    
    print(f"\n  Results saved to tcge_mincut_results.json")


if __name__ == "__main__":
    run_mincut_tests()
