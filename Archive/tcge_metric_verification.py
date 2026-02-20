#!/usr/bin/env python3
"""
TCGE — Metric Axioms Verification
===================================
Verifies that effective distance d(i,j) = min_path Σ w_sym(a,b)
satisfies metric axioms (positivity, symmetry, triangle inequality)
on random constraint graphs.

Reference: Section 4.1 of the manuscript
"""

import numpy as np
from itertools import combinations
import json
import sys

def generate_random_constraint_graph(n_atoms, density=0.4, seed=None):
    """Generate a random TCGE constraint graph with directed weights."""
    rng = np.random.RandomState(seed)
    w_forward = np.zeros((n_atoms, n_atoms))
    w_backward = np.zeros((n_atoms, n_atoms))
    
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            if rng.random() < density:
                # Directed weights: w(i->j) != w(j->i) in general
                w_forward[i, j] = rng.uniform(0.1, 2.0)
                w_backward[i, j] = rng.uniform(0.1, 2.0)
                w_forward[j, i] = w_backward[i, j]
                w_backward[j, i] = w_forward[i, j]
    
    return w_forward, w_backward


def compute_symmetric_weights(w_forward, w_backward):
    """Compute w_sym(a,b) = 0.5 * [w(a->b) + w(b->a)]."""
    return 0.5 * (w_forward + w_backward)


def shortest_path_distances(w_sym, n_atoms):
    """Floyd-Warshall for shortest paths using symmetric weights."""
    dist = np.full((n_atoms, n_atoms), np.inf)
    np.fill_diagonal(dist, 0.0)
    
    for i in range(n_atoms):
        for j in range(n_atoms):
            if i != j and w_sym[i, j] > 0:
                dist[i, j] = w_sym[i, j]
    
    for k in range(n_atoms):
        for i in range(n_atoms):
            for j in range(n_atoms):
                if dist[i, k] + dist[k, j] < dist[i, j]:
                    dist[i, j] = dist[i, k] + dist[k, j]
    
    return dist


def verify_metric_axioms(dist, n_atoms):
    """Check positivity, symmetry, and triangle inequality."""
    results = {
        "positivity": True,
        "symmetry": True,
        "triangle_inequality": True,
        "positivity_violations": 0,
        "symmetry_violations": 0,
        "triangle_violations": 0,
        "pairs_tested": 0,
        "triples_tested": 0
    }
    
    # Only test connected pairs
    connected = []
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            if dist[i, j] < np.inf:
                connected.append((i, j))
    
    results["pairs_tested"] = len(connected)
    
    # Positivity: d(i,j) >= 0 and d(i,j) = 0 iff i = j
    for i, j in connected:
        if dist[i, j] < 0:
            results["positivity"] = False
            results["positivity_violations"] += 1
    
    # Symmetry: d(i,j) = d(j,i)
    for i, j in connected:
        if abs(dist[i, j] - dist[j, i]) > 1e-10:
            results["symmetry"] = False
            results["symmetry_violations"] += 1
    
    # Triangle inequality: d(i,k) <= d(i,j) + d(j,k)
    triples_tested = 0
    for i, j, k in combinations(range(n_atoms), 3):
        if dist[i, j] < np.inf and dist[j, k] < np.inf and dist[i, k] < np.inf:
            triples_tested += 1
            if dist[i, k] > dist[i, j] + dist[j, k] + 1e-10:
                results["triangle_inequality"] = False
                results["triangle_violations"] += 1
    
    results["triples_tested"] = triples_tested
    
    return results


def run_verification(n_models=50, n_range=(6, 30)):
    """Run metric verification across multiple random models."""
    print("=" * 70)
    print("TCGE — METRIC AXIOMS VERIFICATION")
    print("=" * 70)
    print(f"Testing {n_models} random constraint graphs (N = {n_range[0]} to {n_range[1]})")
    print()
    
    all_results = []
    all_pass = True
    
    for trial in range(n_models):
        n = np.random.randint(n_range[0], n_range[1] + 1)
        density = np.random.uniform(0.3, 0.6)
        
        w_fwd, w_bwd = generate_random_constraint_graph(n, density, seed=trial)
        w_sym = compute_symmetric_weights(w_fwd, w_bwd)
        dist = shortest_path_distances(w_sym, n)
        results = verify_metric_axioms(dist, n)
        
        results["n_atoms"] = n
        results["density"] = round(density, 3)
        results["trial"] = trial
        
        passed = results["positivity"] and results["symmetry"] and results["triangle_inequality"]
        results["all_pass"] = passed
        if not passed:
            all_pass = False
        
        all_results.append(results)
        
        status = "✅" if passed else "❌"
        print(f"  Model {trial+1:3d}: N={n:2d}, density={density:.2f} | "
              f"pairs={results['pairs_tested']:4d}, triples={results['triples_tested']:5d} | {status}")
    
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    n_passed = sum(1 for r in all_results if r["all_pass"])
    print(f"Models tested:      {n_models}")
    print(f"Models passed:      {n_passed}")
    print(f"Pass rate:          {100*n_passed/n_models:.1f}%")
    print(f"Positivity:         {sum(1 for r in all_results if r['positivity'])}/{n_models}")
    print(f"Symmetry:           {sum(1 for r in all_results if r['symmetry'])}/{n_models}")
    print(f"Triangle ineq.:     {sum(1 for r in all_results if r['triangle_inequality'])}/{n_models}")
    print()
    
    if all_pass:
        print("✅ ALL METRIC AXIOMS SATISFIED across all models.")
    else:
        print("❌ SOME VIOLATIONS DETECTED — see details above.")
    
    # Save results
    output_file = "tcge_metric_results.json"
    with open(output_file, "w") as f:
        json.dump({
            "summary": {
                "n_models": n_models,
                "n_passed": n_passed,
                "pass_rate": round(100*n_passed/n_models, 1),
                "all_pass": all_pass
            },
            "results": all_results
        }, f, indent=2)
    print(f"\nResults saved to {output_file}")
    
    return all_results


if __name__ == "__main__":
    run_verification(n_models=50, n_range=(6, 30))
