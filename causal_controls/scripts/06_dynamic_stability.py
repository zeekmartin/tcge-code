#!/usr/bin/env python3
"""
TCGE GAP 3 — Dynamic Stability
================================

Test 1: Perturbation resilience
  Start from equilibrium. Add noise ε to all weights. 
  Re-optimize. Measure how far the new equilibrium is from the old.
  
Test 2: Relaxation dynamics
  Track f_S as a function of optimization step (not just final state).
  Show convergence to phase separation.

Test 3: Basin structure
  Start from different initial conditions (random, all-polarized, all-zero).
  Show they converge to the same equilibrium.
"""

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time


def edge_triangles(G):
    tri = {}
    for u, v in G.edges():
        common = len(set(G.neighbors(u)) & set(G.neighbors(v)))
        tri[(u, v)] = common
        tri[(v, u)] = common
    return tri


def optimize_tracked(G, C_protect=0.5, A=1.0, beta_reg=0.01, n_iter=3000,
                     seed=42, w_init=None, track_every=50):
    """Optimizer that returns trajectory of f_S and ⟨|α|⟩."""
    rng = np.random.RandomState(seed)
    edges = list(G.edges())
    M = len(edges)
    W = 1.0
    
    tri_dict = edge_triangles(G)
    tri_arr = np.array([tri_dict.get((u, v), 0) for u, v in edges], dtype=float)
    
    if w_init is not None:
        w = w_init.copy()
    else:
        w = rng.uniform(-0.15, 0.15, M)
    
    lr = 0.02
    trajectory = []
    
    for it in range(n_iter):
        eff = -A + beta_reg + C_protect * tri_arr / (W**2)
        grad = 2 * eff * w
        w -= lr * grad
        
        if it < n_iter * 0.7:
            w += rng.normal(0, 0.02 * (1 - it/n_iter), M)
        
        w = np.clip(w, -W, W)
        
        if it == n_iter // 2:
            lr *= 0.3
        
        if it % track_every == 0:
            alpha = np.abs(w) / W
            trajectory.append({
                'step': it,
                'fS': float(np.mean(alpha < 0.5)),
                'mean_alpha': float(np.mean(alpha)),
                'std_alpha': float(np.std(alpha))
            })
    
    alpha_final = np.abs(w) / W
    return w, alpha_final, tri_arr, trajectory


# ═══════════════════════════════════════════════════════════
# TEST 1: PERTURBATION RESILIENCE
# ═══════════════════════════════════════════════════════════

def test_perturbation():
    print("=" * 60)
    print("  TEST 1: Perturbation resilience")
    print("=" * 60)
    
    N = 80
    G = nx.watts_strogatz_graph(N, 8, 0.1, seed=1000)
    
    # Reach equilibrium
    w_eq, alpha_eq, tri_arr, _ = optimize_tracked(G, C_protect=0.5, seed=42, n_iter=3000)
    fS_eq = np.mean(alpha_eq < 0.5)
    
    # Perturb at different noise levels
    epsilons = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]
    results = []
    
    for eps in epsilons:
        recoveries = []
        for trial in range(10):
            rng = np.random.RandomState(5000 + trial)
            w_perturbed = w_eq + rng.normal(0, eps, len(w_eq))
            w_perturbed = np.clip(w_perturbed, -1, 1)
            
            # Re-optimize from perturbed state (shorter, no noise)
            w_new, alpha_new, _, _ = optimize_tracked(
                G, C_protect=0.5, seed=6000+trial, n_iter=1500,
                w_init=w_perturbed)
            
            fS_new = np.mean(alpha_new < 0.5)
            # Correlation between old and new alpha
            corr = np.corrcoef(alpha_eq, alpha_new)[0, 1]
            recoveries.append({'fS': fS_new, 'corr': corr})
        
        mean_fS = np.mean([r['fS'] for r in recoveries])
        mean_corr = np.mean([r['corr'] for r in recoveries])
        print(f"  ε={eps:.2f}  →  f_S={mean_fS:.3f} (eq: {fS_eq:.3f})  "
              f"corr(α_old, α_new)={mean_corr:.3f}")
        results.append({'eps': eps, 'fS': mean_fS, 'fS_eq': fS_eq, 'corr': mean_corr,
                       'fS_std': np.std([r['fS'] for r in recoveries])})
    
    return results


# ═══════════════════════════════════════════════════════════
# TEST 2: RELAXATION DYNAMICS
# ═══════════════════════════════════════════════════════════

def test_relaxation():
    print("\n" + "=" * 60)
    print("  TEST 2: Relaxation dynamics (trajectory)")
    print("=" * 60)
    
    N = 80
    G = nx.watts_strogatz_graph(N, 8, 0.1, seed=1000)
    
    trajectories = {}
    
    # From random initial (standard)
    _, _, _, traj = optimize_tracked(G, C_protect=0.5, seed=42, n_iter=4000, track_every=20)
    trajectories['random_init'] = traj
    
    # From all-polarized (w = ±0.95)
    rng = np.random.RandomState(99)
    w_pol = rng.choice([-0.95, 0.95], size=G.number_of_edges())
    _, _, _, traj = optimize_tracked(G, C_protect=0.5, seed=42, n_iter=4000,
                                     w_init=w_pol, track_every=20)
    trajectories['polarized_init'] = traj
    
    # From all-zero
    w_zero = np.zeros(G.number_of_edges())
    _, _, _, traj = optimize_tracked(G, C_protect=0.5, seed=42, n_iter=4000,
                                     w_init=w_zero, track_every=20)
    trajectories['zero_init'] = traj
    
    # Without protector (control)
    _, _, _, traj = optimize_tracked(G, C_protect=0.0, seed=42, n_iter=4000, track_every=20)
    trajectories['no_protect'] = traj
    
    for name, traj in trajectories.items():
        final = traj[-1]
        print(f"  {name:<20}  final f_S={final['fS']:.3f}  ⟨|α|⟩={final['mean_alpha']:.3f}")
    
    return trajectories


# ═══════════════════════════════════════════════════════════
# FIGURE
# ═══════════════════════════════════════════════════════════

def make_figure(perturb, trajectories):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Dynamic Stability: Perturbation Resilience and Relaxation',
                 fontsize=13, fontweight='bold', y=1.02)
    
    # ── A: Perturbation resilience ──
    ax = axes[0]
    eps_vals = [r['eps'] for r in perturb]
    corrs = [r['corr'] for r in perturb]
    fS_vals = [r['fS'] for r in perturb]
    fS_eq = perturb[0]['fS_eq']
    
    ax.plot(eps_vals, corrs, 'o-', color='#2563eb', linewidth=2.5, markersize=8,
            label='corr(α_old, α_new)')
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.3)
    ax.set_xlabel('Perturbation amplitude ε', fontsize=11)
    ax.set_ylabel('Correlation with equilibrium', fontsize=11, color='#2563eb')
    ax.tick_params(axis='y', labelcolor='#2563eb')
    ax.set_ylim(0.5, 1.05)
    
    ax2 = ax.twinx()
    ax2.errorbar(eps_vals, fS_vals, yerr=[r['fS_std'] for r in perturb],
                 fmt='s--', color='#16a34a', linewidth=1.5, markersize=6, capsize=3)
    ax2.axhline(fS_eq, color='#16a34a', linestyle=':', alpha=0.5)
    ax2.set_ylabel('f_S after recovery', fontsize=10, color='#16a34a')
    ax2.tick_params(axis='y', labelcolor='#16a34a')
    ax2.set_ylim(0.5, 1.05)
    
    ax.set_title('A. Perturbation resilience', fontweight='bold', fontsize=12)
    ax.legend(fontsize=9, loc='lower left')
    
    # ── B: Relaxation — f_S trajectories ──
    ax = axes[1]
    colors_traj = {'random_init': '#2563eb', 'polarized_init': '#ef4444',
                   'zero_init': '#16a34a', 'no_protect': '#6b7280'}
    labels_traj = {'random_init': 'Random init', 'polarized_init': 'All-polarized init',
                   'zero_init': 'All-zero init', 'no_protect': 'No protector (control)'}
    
    for name, traj in trajectories.items():
        steps = [t['step'] for t in traj]
        fS = [t['fS'] for t in traj]
        ax.plot(steps, fS, linewidth=2.5 if name != 'no_protect' else 1.5,
                color=colors_traj[name], label=labels_traj[name],
                linestyle='--' if name == 'no_protect' else '-')
    
    ax.set_xlabel('Optimization step', fontsize=11)
    ax.set_ylabel('f_S  (fraction proto-spatial)', fontsize=11)
    ax.set_title('B. Relaxation dynamics', fontweight='bold', fontsize=12)
    ax.legend(fontsize=8.5, loc='center right')
    ax.set_ylim(-0.05, 1.05)
    
    # ── C: Relaxation — ⟨|α|⟩ trajectories ──
    ax = axes[2]
    for name, traj in trajectories.items():
        steps = [t['step'] for t in traj]
        alpha_m = [t['mean_alpha'] for t in traj]
        ax.plot(steps, alpha_m, linewidth=2.5 if name != 'no_protect' else 1.5,
                color=colors_traj[name], label=labels_traj[name],
                linestyle='--' if name == 'no_protect' else '-')
    
    ax.set_xlabel('Optimization step', fontsize=11)
    ax.set_ylabel('⟨|α|⟩  (mean polarization)', fontsize=11)
    ax.set_title('C. Convergence of polarization', fontweight='bold', fontsize=12)
    ax.legend(fontsize=8.5, loc='center right')
    ax.set_ylim(-0.05, 1.15)
    
    fig.text(0.5, -0.02,
             'WS graph (N=80, k=8, p=0.1). Perturbation: 10 trials/ε. '
             'Relaxation: single graph, 4000 steps.',
             ha='center', fontsize=8.5, color='gray', style='italic')
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    plt.savefig('/home/claude/fig6_dynamics.png', dpi=200, bbox_inches='tight', facecolor='white')
    print('\nFigure saved → fig6_dynamics.png')


if __name__ == '__main__':
    t0 = time.time()
    perturb = test_perturbation()
    trajectories = test_relaxation()
    make_figure(perturb, trajectories)
    print(f'\nTotal: {time.time()-t0:.0f}s')
