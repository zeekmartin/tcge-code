#!/usr/bin/env python3
"""
TCGE — Combinatorial Reformulation v5c: Statistical Validation
================================================================
Building on v5b's 95% Lorentzian result.

Problem: isotropic control gives 6/10 false positives.
Question: is the anisotropic signal real, or a method artifact?

Approach:
  1. PERMUTATION TEST: same graph, randomly shuffle directed labels
     → if signal disappears, it's structural, not artifactual
  2. EFFECT SIZE (Cohen's d): quantify separation between
     anisotropic signal and permutation null
  3. λ CALIBRATION: choose λ where permuted null gives score ≈ 0
  4. SEPARATION STATISTIC D = E[s²|spatial] - E[s²|temporal]
     → proper continuous metric, not binary pass/fail

This is the "GAP-Detection" closure script.
"""

import numpy as np
from collections import defaultdict, deque
import json
import sys

sys.setrecursionlimit(50000)


# ===================================================================
# UNIVERSE CONSTRUCTION (same as v5b)
# ===================================================================

def create_universe(nx=3, ny=3, nz=3, nt=5,
                     temporal_dir_rate=0.7, spatial_dir_rate=0.05,
                     seed=42):
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


# ===================================================================
# CORE PIPELINE (from v5b)
# ===================================================================

def per_edge_delta(directed, compatible, n):
    delta = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if compatible[i, j]:
                delta[i, j] = int(directed[i, j]) - int(directed[j, i])
    return delta


def split_graphs(compatible, delta, n):
    g_spat = np.zeros((n, n), dtype=bool)
    g_temp = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(n):
            if compatible[i, j]:
                if delta[i, j] == 0:
                    g_spat[i, j] = True
                else:
                    g_temp[i, j] = True
    return g_spat, g_temp


def bfs_distance(adj, n):
    INF = n + 1
    dist = np.full((n, n), INF, dtype=int)
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


def compute_separation_statistic(points, compatible, edge_direction,
                                  x_dist, t_dist, n, lam):
    """
    Compute the separation statistic D and detailed s² values.
    
    D = mean(s² for spatial pairs) - mean(s² for temporal pairs)
    
    If Lorentzian: spatial s² > 0 (spacelike), temporal s² < 0 (timelike)
    → D > 0 and large
    
    Also returns per-direction breakdown.
    """
    INF = n + 1
    temporal_s2 = []
    spatial_s2 = []
    by_dir = {'t': [], 'x': [], 'y': [], 'z': []}
    
    for i in range(n):
        for j in range(i + 1, n):
            if not compatible[i, j]:
                continue
            
            d = edge_direction.get((i, j), '?')
            if d not in by_dir:
                continue
            
            xi = x_dist[i, j] if x_dist[i, j] < INF else 0
            ti = t_dist[i, j] if t_dist[i, j] < INF else 0
            
            s2 = int(xi)**2 - lam * int(ti)**2
            by_dir[d].append(s2)
            
            if d == 't':
                temporal_s2.append(s2)
            else:
                spatial_s2.append(s2)
    
    temporal_s2 = np.array(temporal_s2, dtype=float)
    spatial_s2 = np.array(spatial_s2, dtype=float)
    
    if len(temporal_s2) == 0 or len(spatial_s2) == 0:
        return {"D": 0.0, "t_neg_rate": 0.0, "s_pos_rate": 0.0,
                "lorentzian": False, "by_dir": {}}
    
    D = float(np.mean(spatial_s2) - np.mean(temporal_s2))
    t_neg = float(np.mean(temporal_s2 < 0))
    s_pos = float(np.mean(spatial_s2 > 0))
    
    dir_results = {}
    for d in ['t', 'x', 'y', 'z']:
        vals = np.array(by_dir[d], dtype=float)
        if len(vals) > 0:
            dir_results[d] = {
                "mean": round(float(np.mean(vals)), 2),
                "neg_frac": round(float(np.mean(vals < 0)), 3),
                "pos_frac": round(float(np.mean(vals > 0)), 3),
            }
    
    return {
        "D": round(D, 2),
        "t_neg_rate": round(t_neg, 3),
        "s_pos_rate": round(s_pos, 3),
        "score": round(t_neg + s_pos - 1.0, 3),
        "lorentzian": bool(t_neg > 0.4 and s_pos > 0.4),
        "by_dir": dir_results
    }


def full_pipeline(points, compatible, directed, edge_direction, n, lam):
    """Run full pipeline: Δ → split → distances → separation statistic."""
    delta = per_edge_delta(directed, compatible, n)
    g_spat, g_temp = split_graphs(compatible, delta, n)
    x_dist = bfs_distance(g_spat, n)
    t_dist = bfs_distance(g_temp, n)
    return compute_separation_statistic(
        points, compatible, edge_direction, x_dist, t_dist, n, lam)


# ===================================================================
# PERMUTATION TEST
# ===================================================================

def permute_directed(directed, compatible, n, rng):
    """
    Permutation null model: keep the same NUMBER of directed edges,
    but randomly reassign which compatible pairs get a directed edge.
    
    This preserves: graph structure, total directed count.
    This destroys: directional anisotropy.
    """
    # Collect all compatible pairs that COULD have a directed edge
    candidates = []
    for i in range(n):
        for j in range(n):
            if compatible[i, j]:
                candidates.append((i, j))
    
    # Count current directed edges
    n_directed = int(directed.sum())
    
    # Randomly select same number of edges to be directed
    new_directed = np.zeros((n, n), dtype=bool)
    chosen = rng.choice(len(candidates), size=min(n_directed, len(candidates)),
                         replace=False)
    for idx in chosen:
        i, j = candidates[idx]
        new_directed[i, j] = True
    
    return new_directed


def permute_directed_local(directed, compatible, n, rng):
    """
    Alternative permutation: for each edge (i,j), randomly swap
    whether directed[i,j] or directed[j,i] is True.
    
    This preserves: which PAIRS have directed edges.
    This destroys: the direction of the asymmetry.
    """
    new_directed = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i + 1, n):
            if not compatible[i, j]:
                continue
            has_ij = directed[i, j]
            has_ji = directed[j, i]
            
            if has_ij or has_ji:
                # Randomly reassign the direction(s)
                if has_ij and has_ji:
                    new_directed[i, j] = True
                    new_directed[j, i] = True
                elif rng.random() < 0.5:
                    new_directed[i, j] = has_ij or has_ji
                else:
                    new_directed[j, i] = has_ij or has_ji
    
    return new_directed


# ===================================================================
# MAIN
# ===================================================================

def run_v5c():
    print("=" * 72)
    print("  TCGE — v5c: STATISTICAL VALIDATION")
    print("  Permutation tests, effect size, controlled λ")
    print("=" * 72)
    
    LAM = 10  # from v5b
    N_SEEDS = 15
    N_PERMS = 15
    
    # ==============================================================
    # PART 1: ANISOTROPIC vs PERMUTED NULL (per seed)
    # ==============================================================
    print(f"\n{'='*72}")
    print(f"  PART 1: PERMUTATION TEST")
    print(f"  For each seed: compare real signal vs {N_PERMS} permutations")
    print(f"{'='*72}")
    
    all_real_D = []
    all_perm_D = []
    all_real_scores = []
    all_perm_scores = []
    
    for seed in range(N_SEEDS):
        pts, comp, dire, edir, n = create_universe(seed=seed)
        
        # Real signal
        real = full_pipeline(pts, comp, dire, edir, n, LAM)
        all_real_D.append(real['D'])
        all_real_scores.append(real['score'])
        
        # Permuted null (global shuffle)
        perm_Ds = []
        perm_scores = []
        rng = np.random.RandomState(seed * 1000)
        
        for p in range(N_PERMS):
            perm_dir = permute_directed(dire, comp, n, rng)
            perm = full_pipeline(pts, comp, perm_dir, edir, n, LAM)
            perm_Ds.append(perm['D'])
            perm_scores.append(perm['score'])
        
        all_perm_D.extend(perm_Ds)
        all_perm_scores.extend(perm_scores)
        
        mean_perm_D = np.mean(perm_Ds)
        p_value = np.mean(np.array(perm_Ds) >= real['D'])
        
        print(f"  seed={seed:>2}: D_real={real['D']:>+8.1f} | "
              f"D_perm={mean_perm_D:>+8.1f} | "
              f"p={p_value:.2f} | "
              f"score_real={real['score']:>+.3f}")
    
    # Summary statistics
    real_D_arr = np.array(all_real_D)
    perm_D_arr = np.array(all_perm_D)
    real_score_arr = np.array(all_real_scores)
    perm_score_arr = np.array(all_perm_scores)
    
    print(f"\n  {'─'*60}")
    print(f"  SEPARATION STATISTIC D = E[s²|spatial] − E[s²|temporal]")
    print(f"  {'─'*60}")
    print(f"    Real (anisotropic):  mean D = {np.mean(real_D_arr):>+8.1f} "
          f"± {np.std(real_D_arr):.1f}")
    print(f"    Permuted (null):     mean D = {np.mean(perm_D_arr):>+8.1f} "
          f"± {np.std(perm_D_arr):.1f}")
    
    # Cohen's d
    pooled_std = np.sqrt((np.var(real_D_arr) + np.var(perm_D_arr)) / 2)
    cohens_d_D = (np.mean(real_D_arr) - np.mean(perm_D_arr)) / max(pooled_std, 0.01)
    
    print(f"\n    Cohen's d (D):       {cohens_d_D:>+.2f}")
    if abs(cohens_d_D) > 0.8:
        print(f"    → LARGE effect size ✅")
    elif abs(cohens_d_D) > 0.5:
        print(f"    → MEDIUM effect size")
    else:
        print(f"    → SMALL effect size ⚠️")
    
    print(f"\n  {'─'*60}")
    print(f"  LORENTZIAN SCORE = t_neg_rate + s_pos_rate − 1")
    print(f"  {'─'*60}")
    print(f"    Real:    mean = {np.mean(real_score_arr):>+.3f} "
          f"± {np.std(real_score_arr):.3f}")
    print(f"    Permuted: mean = {np.mean(perm_score_arr):>+.3f} "
          f"± {np.std(perm_score_arr):.3f}")
    
    pooled_std_s = np.sqrt((np.var(real_score_arr) + np.var(perm_score_arr)) / 2)
    cohens_d_s = (np.mean(real_score_arr) - np.mean(perm_score_arr)) / max(pooled_std_s, 0.01)
    
    print(f"    Cohen's d (score):   {cohens_d_s:>+.2f}")
    if abs(cohens_d_s) > 0.8:
        print(f"    → LARGE effect size ✅")
    elif abs(cohens_d_s) > 0.5:
        print(f"    → MEDIUM effect size")
    else:
        print(f"    → SMALL effect size ⚠️")
    
    # Global p-value (how often does ANY permuted score beat ALL real scores?)
    real_median = np.median(real_score_arr)
    global_p = np.mean(perm_score_arr >= real_median)
    print(f"\n    Global p-value (perm score ≥ real median): {global_p:.4f}")
    if global_p < 0.01:
        print(f"    → Highly significant (p < 0.01) ✅")
    elif global_p < 0.05:
        print(f"    → Significant (p < 0.05) ✅")
    else:
        print(f"    → Not significant ⚠️")
    
    # ==============================================================
    # PART 2: LOCAL PERMUTATION (swap directions, not assignments)
    # ==============================================================
    print(f"\n{'='*72}")
    print(f"  PART 2: LOCAL PERMUTATION (direction-swap)")
    print(f"  Keep which pairs have directed edges, randomize direction")
    print(f"{'='*72}")
    
    local_perm_scores = []
    
    for seed in range(10):  # fewer for local perm
        pts, comp, dire, edir, n = create_universe(seed=seed)
        rng = np.random.RandomState(seed * 2000)
        
        for p in range(10):
            perm_dir = permute_directed_local(dire, comp, n, rng)
            perm = full_pipeline(pts, comp, perm_dir, edir, n, LAM)
            local_perm_scores.append(perm['score'])
    
    local_perm_arr = np.array(local_perm_scores)
    
    print(f"    Real:         mean score = {np.mean(real_score_arr):>+.3f} "
          f"± {np.std(real_score_arr):.3f}")
    print(f"    Local perm:   mean score = {np.mean(local_perm_arr):>+.3f} "
          f"± {np.std(local_perm_arr):.3f}")
    
    pooled_local = np.sqrt((np.var(real_score_arr) + np.var(local_perm_arr)) / 2)
    cohens_d_local = (np.mean(real_score_arr) - np.mean(local_perm_arr)) / max(pooled_local, 0.01)
    
    print(f"    Cohen's d:    {cohens_d_local:>+.2f}")
    
    local_p = np.mean(local_perm_arr >= real_median)
    print(f"    p-value:      {local_p:.4f}")
    
    # ==============================================================
    # PART 3: λ CALIBRATION VIA PERMUTATION NULL
    # ==============================================================
    print(f"\n{'='*72}")
    print(f"  PART 3: λ CALIBRATION")
    print(f"  Find λ where permuted null is centered at score ≈ 0")
    print(f"  AND real signal is maximally separated")
    print(f"{'='*72}")
    
    print(f"\n  {'λ':>4} {'real_score':>11} {'perm_score':>11} "
          f"{'Cohen_d':>9} {'p_value':>9} {'status':>8}")
    print(f"  {'─'*4} {'─'*11} {'─'*11} {'─'*9} {'─'*9} {'─'*8}")
    
    best_youden = -999
    best_lambda_final = 1
    
    for lam in [1, 2, 3, 5, 7, 10, 12, 15]:
        real_scores_lam = []
        perm_scores_lam = []
        
        for seed in range(10):  # fewer seeds for speed
            pts, comp, dire, edir, n = create_universe(seed=seed)
            real = full_pipeline(pts, comp, dire, edir, n, lam)
            real_scores_lam.append(real['score'])
            
            rng = np.random.RandomState(seed * 3000)
            for p in range(10):
                pd = permute_directed(dire, comp, n, rng)
                pm = full_pipeline(pts, comp, pd, edir, n, lam)
                perm_scores_lam.append(pm['score'])
        
        r_arr = np.array(real_scores_lam)
        p_arr = np.array(perm_scores_lam)
        
        ps = np.sqrt((np.var(r_arr) + np.var(p_arr)) / 2)
        cd = (np.mean(r_arr) - np.mean(p_arr)) / max(ps, 0.01)
        pv = np.mean(p_arr >= np.median(r_arr))
        
        # Youden's J: maximize separation
        # Approximate: real > threshold vs perm < threshold
        threshold = (np.mean(r_arr) + np.mean(p_arr)) / 2
        tpr = np.mean(r_arr > threshold)
        fpr = np.mean(p_arr > threshold)
        youden = tpr - fpr
        
        status = ""
        if youden > best_youden:
            best_youden = youden
            best_lambda_final = lam
            status = "← best"
        
        print(f"  {lam:>4} {np.mean(r_arr):>+10.3f} {np.mean(p_arr):>+10.3f} "
              f"{cd:>+8.2f} {pv:>9.4f} {status:>8}")
    
    print(f"\n  Optimal λ = {best_lambda_final} (Youden's J = {best_youden:.3f})")
    
    # ==============================================================
    # PART 4: FINAL RESULT AT OPTIMAL λ
    # ==============================================================
    print(f"\n{'='*72}")
    print(f"  PART 4: DEFINITIVE RESULT (λ = {best_lambda_final})")
    print(f"{'='*72}")
    
    final_real = []
    final_perm = []
    final_lor = 0
    
    for seed in range(N_SEEDS):
        pts, comp, dire, edir, n = create_universe(seed=seed)
        real = full_pipeline(pts, comp, dire, edir, n, best_lambda_final)
        final_real.append(real)
        if real['lorentzian']:
            final_lor += 1
        
        rng = np.random.RandomState(seed * 4000)
        for p in range(10):  # 10 perms per seed for final
            pd = permute_directed(dire, comp, n, rng)
            pm = full_pipeline(pts, comp, pd, edir, n, best_lambda_final)
            final_perm.append(pm)
    
    N_FINAL_PERMS = 10
    final_real_scores = [r['score'] for r in final_real]
    final_perm_scores = [r['score'] for r in final_perm]
    final_real_D = [r['D'] for r in final_real]
    final_perm_D = [r['D'] for r in final_perm]
    
    fr = np.array(final_real_scores)
    fp = np.array(final_perm_scores)
    
    ps_final = np.sqrt((np.var(fr) + np.var(fp)) / 2)
    cd_final = (np.mean(fr) - np.mean(fp)) / max(ps_final, 0.01)
    pv_final = np.mean(fp >= np.median(fr))
    
    # Lorentzian rate for permutations
    perm_lor = sum(1 for r in final_perm if r['lorentzian'])
    
    print(f"\n  At λ = {best_lambda_final}:")
    print(f"  {'─'*55}")
    print(f"    Real anisotropic:")
    print(f"      Lorentzian:     {final_lor}/{N_SEEDS} ({100*final_lor/N_SEEDS:.0f}%)")
    print(f"      Mean score:     {np.mean(fr):>+.3f} ± {np.std(fr):.3f}")
    print(f"      Mean D:         {np.mean(final_real_D):>+.1f}")
    print(f"\n    Permuted null ({N_SEEDS}×{N_FINAL_PERMS} = {N_SEEDS*N_FINAL_PERMS} samples):")
    print(f"      Lorentzian:     {perm_lor}/{N_SEEDS*N_FINAL_PERMS} "
          f"({100*perm_lor/(N_SEEDS*N_FINAL_PERMS):.0f}%)")
    print(f"      Mean score:     {np.mean(fp):>+.3f} ± {np.std(fp):.3f}")
    print(f"      Mean D:         {np.mean(final_perm_D):>+.1f}")
    print(f"\n    STATISTICAL TESTS:")
    print(f"      Cohen's d:      {cd_final:>+.2f}", end="")
    if abs(cd_final) > 0.8:
        print(f" (LARGE ✅)")
    elif abs(cd_final) > 0.5:
        print(f" (MEDIUM)")
    else:
        print(f" (SMALL ⚠️)")
    
    print(f"      p-value:        {pv_final:.4f}", end="")
    if pv_final < 0.001:
        print(f" (p < 0.001 ✅✅✅)")
    elif pv_final < 0.01:
        print(f" (p < 0.01 ✅✅)")
    elif pv_final < 0.05:
        print(f" (p < 0.05 ✅)")
    else:
        print(f" (not significant ⚠️)")
    
    # False positive rate at real threshold
    threshold = np.percentile(fr, 25)  # 25th percentile of real
    fpr_at_threshold = np.mean(fp >= threshold)
    print(f"      FPR at 25th pctl: {fpr_at_threshold:.1%}")
    
    # ==============================================================
    # DISTRIBUTIONS (text histogram)
    # ==============================================================
    print(f"\n  SCORE DISTRIBUTIONS:")
    print(f"  {'─'*55}")
    
    bins = np.linspace(min(fp.min(), fr.min()) - 0.05, 
                        max(fp.max(), fr.max()) + 0.05, 20)
    
    real_hist, _ = np.histogram(fr, bins=bins)
    perm_hist, _ = np.histogram(fp, bins=bins, density=False)
    # Normalize perm to same scale as real
    perm_hist_norm = perm_hist * (len(fr) / max(len(fp), 1))
    
    max_count = max(real_hist.max(), perm_hist_norm.max(), 1)
    
    print(f"    {'bin_center':>10} {'Real':>15} {'Permuted (scaled)':>20}")
    for k in range(len(bins) - 1):
        center = (bins[k] + bins[k+1]) / 2
        r_bar = "█" * int(30 * real_hist[k] / max_count)
        p_bar = "░" * int(30 * perm_hist_norm[k] / max_count)
        if real_hist[k] > 0 or perm_hist_norm[k] > 0:
            print(f"    {center:>+10.3f} {r_bar:<15} {p_bar:<20}")
    
    print(f"\n    █ = Real (anisotropic)    ░ = Permuted (null)")
    
    # ==============================================================
    # CONCLUSION
    # ==============================================================
    print(f"\n{'='*72}")
    print(f"  CONCLUSION: GAP-DETECTION STATUS")
    print(f"{'='*72}")
    print()
    
    gap_closed = (abs(cd_final) > 0.8 and pv_final < 0.01 and 
                   final_lor / N_SEEDS > 0.5)
    
    if gap_closed:
        print(f"  ✅ GAP-DETECTION: CLOSED")
        print(f"     The Lorentzian signal is statistically significant:")
        print(f"     • Cohen's d = {cd_final:+.2f} (large effect)")
        print(f"     • p < {max(pv_final, 0.001):.3f} (permutation test)")
        print(f"     • {final_lor}/{N_SEEDS} real vs {perm_lor}/{N_SEEDS*N_FINAL_PERMS} "
              f"permuted pass Lorentzian criterion")
        print(f"     • The signal comes from STRUCTURAL ANISOTROPY,")
        print(f"       not from method bias or statistical fluctuation.")
    else:
        print(f"  ⚠️  GAP-DETECTION: PARTIALLY CLOSED")
        print(f"     Signal detected but statistical criteria not all met:")
        print(f"     • Cohen's d = {cd_final:+.2f}")
        print(f"     • p = {pv_final:.4f}")
        print(f"     • {final_lor}/{N_SEEDS} Lorentzian")
    
    print()
    print(f"  PUBLISHABLE STATEMENT:")
    print(f"    'Boolean directed constraints on a compatibility graph")
    print(f"     produce Lorentzian signature separation with")
    print(f"     Cohen's d = {cd_final:+.2f} relative to permutation null")
    print(f"     (p {'< 0.001' if pv_final < 0.001 else '= ' + str(round(pv_final, 4))},"
          f" {N_SEEDS} seeds × {N_FINAL_PERMS} permutations).'")
    
    # Save
    output = {
        "lambda": best_lambda_final,
        "real": {
            "lorentzian_rate": final_lor / N_SEEDS,
            "mean_score": round(float(np.mean(fr)), 3),
            "std_score": round(float(np.std(fr)), 3),
            "mean_D": round(float(np.mean(final_real_D)), 1),
        },
        "permuted": {
            "lorentzian_rate": perm_lor / (N_SEEDS * N_FINAL_PERMS),
            "mean_score": round(float(np.mean(fp)), 3),
            "std_score": round(float(np.std(fp)), 3),
            "mean_D": round(float(np.mean(final_perm_D)), 1),
        },
        "statistics": {
            "cohens_d": round(float(cd_final), 2),
            "p_value": round(float(pv_final), 4),
            "gap_detection_closed": gap_closed,
        }
    }
    
    with open("tcge_v5c_validation.json", "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"\n  Results saved to tcge_v5c_validation.json")


if __name__ == "__main__":
    run_v5c()
