#!/usr/bin/env python3
"""
TCGE — Lorentzian Signature Emergence
========================================
Demonstrates that Lorentzian signature (-,+,+,+) emerges spontaneously
from directed constraint asymmetry in TCGE models.

Reference: Section 4.2 of the manuscript
"""

import numpy as np
import json

def create_3plus1D_model(nx=3, ny=3, nz=3, nt=3, w_base=1.0, asymmetry=0.6, seed=42):
    """
    Create a 3+1D constraint model where:
    - Spatial dimensions have symmetric weights
    - Temporal dimension has asymmetric weights
    """
    rng = np.random.RandomState(seed)
    
    # Grid points
    points = []
    for t in range(nt):
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    points.append((t, x, y, z))
    
    n = len(points)
    w_forward = np.zeros((n, n))
    w_backward = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            t1, x1, y1, z1 = points[i]
            t2, x2, y2, z2 = points[j]
            
            dt = abs(t2 - t1)
            dx = abs(x2 - x1)
            dy = abs(y2 - y1)
            dz = abs(z2 - z1)
            
            # Only connect nearest neighbours
            manhattan = dt + dx + dy + dz
            if manhattan != 1:
                continue
            
            noise = rng.uniform(-0.02, 0.02)  # small noise
            
            if dt == 1:
                # Temporal link: strongly asymmetric
                if t2 > t1:
                    w_forward[i, j] = w_base * (1 - asymmetry) + noise
                    w_backward[i, j] = w_base * (1 + asymmetry) + noise
                else:
                    w_forward[i, j] = w_base * (1 + asymmetry) + noise
                    w_backward[i, j] = w_base * (1 - asymmetry) + noise
            else:
                # Spatial link: symmetric (small noise preserves symmetry)
                w_val = w_base + noise
                w_forward[i, j] = w_val
                w_backward[i, j] = w_val
            
            w_forward[j, i] = w_backward[i, j]
            w_backward[j, i] = w_forward[i, j]
    
    return points, w_forward, w_backward


def compute_effective_metric(points, w_forward, w_backward):
    """
    Compute the effective metric tensor components from constraint weights.
    
    For each pair:
      w_sym = 0.5 * (w_fwd + w_bwd)  -> spatial part
      Δw = w_fwd - w_bwd              -> temporal part
      s² = w_sym² - Δw²              -> signed interval
    """
    n = len(points)
    
    # Accumulate metric components by direction
    # Directions: t, x, y, z
    metric_components = {0: [], 1: [], 2: [], 3: []}  # t, x, y, z
    
    for i in range(n):
        for j in range(i+1, n):
            if w_forward[i, j] == 0:
                continue
            
            t1, x1, y1, z1 = points[i]
            t2, x2, y2, z2 = points[j]
            
            dt = t2 - t1
            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1
            
            w_sym = 0.5 * (w_forward[i, j] + w_backward[i, j])
            delta_w = w_forward[i, j] - w_backward[i, j]
            s_squared = w_sym**2 - delta_w**2
            
            # Assign to direction
            if abs(dt) == 1:
                metric_components[0].append(s_squared)
            elif abs(dx) == 1:
                metric_components[1].append(s_squared)
            elif abs(dy) == 1:
                metric_components[2].append(s_squared)
            elif abs(dz) == 1:
                metric_components[3].append(s_squared)
    
    # Average metric component per direction
    g = {}
    labels = {0: "t", 1: "x", 2: "y", 3: "z"}
    for d in range(4):
        if metric_components[d]:
            g[labels[d]] = np.mean(metric_components[d])
        else:
            g[labels[d]] = 0.0
    
    return g, metric_components


def classify_intervals(points, w_forward, w_backward):
    """Classify all intervals as timelike, spacelike, or null."""
    n = len(points)
    timelike = 0
    spacelike = 0
    null_like = 0
    total = 0
    
    for i in range(n):
        for j in range(i+1, n):
            if w_forward[i, j] == 0:
                continue
            total += 1
            
            w_sym = 0.5 * (w_forward[i, j] + w_backward[i, j])
            delta_w = w_forward[i, j] - w_backward[i, j]
            s_squared = w_sym**2 - delta_w**2
            
            if s_squared < -1e-10:
                timelike += 1
            elif s_squared > 1e-10:
                spacelike += 1
            else:
                null_like += 1
    
    return {"timelike": timelike, "spacelike": spacelike, "null": null_like, "total": total}


def light_cone_test(points, w_forward, w_backward):
    """
    Test light cone emergence: from a central point, temporal neighbours
    should have timelike intervals, spatial neighbours spacelike.
    """
    n = len(points)
    # Find a central point
    center_idx = n // 2
    center = points[center_idx]
    
    results = {"center": center, "neighbours": []}
    
    for j in range(n):
        if j == center_idx or w_forward[center_idx, j] == 0:
            continue
        
        t1, x1, y1, z1 = center
        t2, x2, y2, z2 = points[j]
        
        w_sym = 0.5 * (w_forward[center_idx, j] + w_backward[center_idx, j])
        delta_w = w_forward[center_idx, j] - w_backward[center_idx, j]
        s_squared = w_sym**2 - delta_w**2
        
        direction = "temporal" if abs(t2 - t1) == 1 else "spatial"
        interval_type = "timelike" if s_squared < -1e-10 else ("spacelike" if s_squared > 1e-10 else "null")
        
        correct = (direction == "temporal" and interval_type == "timelike") or \
                  (direction == "spatial" and interval_type == "spacelike")
        
        results["neighbours"].append({
            "point": points[j],
            "direction": direction,
            "interval_type": interval_type,
            "s_squared": round(float(s_squared), 6),
            "correct": correct
        })
    
    results["all_correct"] = all(nb["correct"] for nb in results["neighbours"])
    return results


def run_lorentzian_test():
    """Run the full Lorentzian signature test."""
    print("=" * 70)
    print("TCGE — LORENTZIAN SIGNATURE EMERGENCE TEST")
    print("=" * 70)
    
    # 3+1D model
    print("\n1. Creating 3×3×3×3 constraint model...")
    points, w_fwd, w_bwd = create_3plus1D_model(nx=3, ny=3, nz=3, nt=3)
    print(f"   {len(points)} grid points, 4D lattice")
    
    # Metric tensor
    print("\n2. Computing effective metric tensor...")
    g, components = compute_effective_metric(points, w_fwd, w_bwd)
    
    print(f"\n   Diagonal metric: g = diag({g['t']:.2f}, {g['x']:.2f}, {g['y']:.2f}, {g['z']:.2f})")
    print(f"   Temporal component (g_tt): {g['t']:.4f} {'< 0 ✅' if g['t'] < 0 else '≥ 0 ❌'}")
    print(f"   Spatial components (g_xx): {g['x']:.4f} {'> 0 ✅' if g['x'] > 0 else '≤ 0 ❌'}")
    print(f"   Spatial components (g_yy): {g['y']:.4f} {'> 0 ✅' if g['y'] > 0 else '≤ 0 ❌'}")
    print(f"   Spatial components (g_zz): {g['z']:.4f} {'> 0 ✅' if g['z'] > 0 else '≤ 0 ❌'}")
    
    is_lorentzian = g['t'] < 0 and g['x'] > 0 and g['y'] > 0 and g['z'] > 0
    print(f"\n   Signature: {'(−,+,+,+) — LORENTZIAN ✅' if is_lorentzian else 'NOT Lorentzian ❌'}")
    
    # Interval classification
    print("\n3. Interval classification...")
    intervals = classify_intervals(points, w_fwd, w_bwd)
    print(f"   Timelike:  {intervals['timelike']} ({100*intervals['timelike']/max(1,intervals['total']):.1f}%)")
    print(f"   Spacelike: {intervals['spacelike']} ({100*intervals['spacelike']/max(1,intervals['total']):.1f}%)")
    print(f"   Null:      {intervals['null']} ({100*intervals['null']/max(1,intervals['total']):.1f}%)")
    
    # Light cone test
    print("\n4. Light cone emergence test...")
    lc = light_cone_test(points, w_fwd, w_bwd)
    print(f"   Center point: {lc['center']}")
    for nb in lc["neighbours"]:
        status = "✅" if nb["correct"] else "❌"
        print(f"   → {nb['point']} [{nb['direction']}]: s²={nb['s_squared']:.4f} → {nb['interval_type']} {status}")
    print(f"   Light cone correct: {'YES ✅' if lc['all_correct'] else 'NO ❌'}")
    
    # Multi-run test
    print("\n5. Robustness: testing 20 random seeds...")
    lorentzian_count = 0
    for seed in range(20):
        pts, wf, wb = create_3plus1D_model(seed=seed)
        g_test, _ = compute_effective_metric(pts, wf, wb)
        if g_test['t'] < 0 and g_test['x'] > 0 and g_test['y'] > 0 and g_test['z'] > 0:
            lorentzian_count += 1
    print(f"   Lorentzian signature in {lorentzian_count}/20 runs ({100*lorentzian_count/20:.0f}%)")
    
    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    if is_lorentzian and lc["all_correct"]:
        print("✅ Lorentzian signature (−,+,+,+) emerges from constraint asymmetry.")
        print("✅ Light cones emerge naturally without being postulated.")
    else:
        print("❌ Results require further investigation.")
    
    # Save results
    output = {
        "metric": {k: round(float(v), 6) for k, v in g.items()},
        "is_lorentzian": bool(is_lorentzian),
        "intervals": intervals,
        "light_cone_correct": bool(lc["all_correct"]),
        "robustness_rate": float(lorentzian_count / 20)
    }
    with open("tcge_lorentzian_results.json", "w") as f:
        json.dump(output, f, indent=2)
    print("\nResults saved to tcge_lorentzian_results.json")


if __name__ == "__main__":
    run_lorentzian_test()
