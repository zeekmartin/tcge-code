#!/usr/bin/env python3
"""
TCGE Reviewer Checklist — Dose-response + Triangle sweep
=========================================================

Test 1: Dose-response on C_protect (0 → 1) at fixed high-tri (WS)
  → f_S should rise monotonically with C_protect

Test 2: Triangle level sweep (vary p_ws from 0 to 0.9)  
  → At protect ON: f_S drops as clustering drops
  → At protect OFF: f_S stays at 0 regardless

Test 3: Seed stability (already done: 25 seeds in double dissociation)

N=80, 15 seeds per point. Should run in <90s total.
"""

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import json


def edge_triangles(G):
    tri = {}
    for u, v in G.edges():
        common = len(set(G.neighbors(u)) & set(G.neighbors(v)))
        tri[(u, v)] = common
        tri[(v, u)] = common
    return tri


def optimize(G, C_protect=0.5, A=1.0, beta_reg=0.01, n_iter=3000, seed=42):
    rng = np.random.RandomState(seed)
    edges = list(G.edges())
    M = len(edges)
    W = 1.0
    tri_dict = edge_triangles(G)
    tri_arr = np.array([tri_dict.get((u, v), 0) for u, v in edges], dtype=float)
    w = rng.uniform(-0.15, 0.15, M)
    lr = 0.02
    noise_scale = 0.02
    for it in range(n_iter):
        eff = -A + beta_reg + C_protect * tri_arr / (W**2)
        grad = 2 * eff * w
        w -= lr * grad
        if it < n_iter * 0.7:
            noise = rng.normal(0, noise_scale * (1 - it / n_iter), M)
            w += noise
        w = np.clip(w, -W, W)
        if it == n_iter // 2:
            lr *= 0.3
    alpha = np.abs(w) / W
    return alpha, tri_arr


def compute_metrics(alpha, tri_arr):
    frac_spatial = float(np.mean(alpha < 0.5))
    mu1 = float(np.std(alpha))
    mean_alpha = float(np.mean(alpha))
    mean_tri = float(np.mean(tri_arr))
    
    median_tri = np.median(tri_arr)
    if median_tri == 0:
        low_mask = tri_arr == 0
        high_mask = tri_arr > 0
    else:
        low_mask = tri_arr <= median_tri
        high_mask = tri_arr > median_tri
    
    S = float(np.mean(alpha[low_mask]) - np.mean(alpha[high_mask])) if np.any(high_mask) and np.any(low_mask) else 0.0
    
    return {'S': S, 'mu1': mu1, 'mean_alpha': mean_alpha, 
            'mean_tri': mean_tri, 'frac_spatial': frac_spatial}


# ═══════════════════════════════════════════════════════════════
# TEST 1: DOSE-RESPONSE on C_protect
# ═══════════════════════════════════════════════════════════════

def test_dose_response(N=80, n_seeds=15):
    """Sweep C_protect from 0 to 1.5 on WS graph."""
    print("═" * 60)
    print("  TEST 1: Dose-response (C_protect sweep)")
    print("═" * 60)
    
    C_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5]
    results = []
    
    for C in C_values:
        fs_list, alpha_list, mu1_list = [], [], []
        for seed in range(n_seeds):
            G = nx.watts_strogatz_graph(N, 8, 0.1, seed=1000 + seed)
            alpha, tri_arr = optimize(G, C_protect=C, seed=3000 + seed)
            m = compute_metrics(alpha, tri_arr)
            fs_list.append(m['frac_spatial'])
            alpha_list.append(m['mean_alpha'])
            mu1_list.append(m['mu1'])
        
        r = {'C': C, 
             'fs_mean': np.mean(fs_list), 'fs_std': np.std(fs_list),
             'alpha_mean': np.mean(alpha_list), 'alpha_std': np.std(alpha_list),
             'mu1_mean': np.mean(mu1_list), 'mu1_std': np.std(mu1_list)}
        results.append(r)
        print(f"  C={C:.2f}  →  f_S={r['fs_mean']:.3f}±{r['fs_std']:.3f}  "
              f"⟨|α|⟩={r['alpha_mean']:.3f}  μ₁={r['mu1_mean']:.3f}")
    
    return results


# ═══════════════════════════════════════════════════════════════
# TEST 2: TRIANGLE LEVEL SWEEP
# ═══════════════════════════════════════════════════════════════

def test_triangle_sweep(N=80, n_seeds=15):
    """Sweep p_ws from 0.0 to 0.9 (controls clustering level).
    p_ws=0 → max clustering, p_ws=1 → random-like.
    Compare protect ON vs OFF.
    """
    print("\n" + "═" * 60)
    print("  TEST 2: Triangle level sweep (p_ws)")
    print("═" * 60)
    
    p_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9]
    results_on = []
    results_off = []
    
    for p_ws in p_values:
        # Protect ON
        fs_on, tri_on = [], []
        for seed in range(n_seeds):
            G = nx.watts_strogatz_graph(N, 8, p_ws, seed=1000 + seed)
            alpha, tri_arr = optimize(G, C_protect=0.5, seed=3000 + seed)
            m = compute_metrics(alpha, tri_arr)
            fs_on.append(m['frac_spatial'])
            tri_on.append(m['mean_tri'])
        
        # Protect OFF
        fs_off = []
        for seed in range(n_seeds):
            G = nx.watts_strogatz_graph(N, 8, p_ws, seed=1000 + seed)
            alpha, tri_arr = optimize(G, C_protect=0.0, seed=3000 + seed)
            m = compute_metrics(alpha, tri_arr)
            fs_off.append(m['frac_spatial'])
        
        r_on = {'p_ws': p_ws, 'fs_mean': np.mean(fs_on), 'fs_std': np.std(fs_on),
                'tri_mean': np.mean(tri_on)}
        r_off = {'p_ws': p_ws, 'fs_mean': np.mean(fs_off), 'fs_std': np.std(fs_off),
                 'tri_mean': np.mean(tri_on)}
        results_on.append(r_on)
        results_off.append(r_off)
        
        print(f"  p_ws={p_ws:.2f}  ⟨tri⟩={r_on['tri_mean']:.2f}  "
              f"ON: f_S={r_on['fs_mean']:.3f}±{r_on['fs_std']:.3f}  "
              f"OFF: f_S={r_off['fs_mean']:.3f}")
    
    return results_on, results_off


# ═══════════════════════════════════════════════════════════════
# FIGURE
# ═══════════════════════════════════════════════════════════════

def make_figure(dose_results, tri_on, tri_off):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Reviewer Checklist: Dose-Response and Triangle Sweep',
                 fontsize=13, fontweight='bold', y=1.02)
    
    # ── Panel A: Dose-response ──
    ax = axes[0]
    Cs = [r['C'] for r in dose_results]
    fs = [r['fs_mean'] for r in dose_results]
    fs_err = [r['fs_std'] for r in dose_results]
    alphas = [r['alpha_mean'] for r in dose_results]
    
    line1, = ax.plot(Cs, fs, 'o-', color='#2563eb', linewidth=2.5, markersize=7,
                     label='f_S (proto-spatial)')
    ax.fill_between(Cs, [f - e for f, e in zip(fs, fs_err)],
                    [f + e for f, e in zip(fs, fs_err)], alpha=0.2, color='#2563eb')
    ax.set_xlabel('C_protect', fontsize=12)
    ax.set_ylabel('f_S  (fraction proto-spatial)', fontsize=11, color='#2563eb')
    ax.tick_params(axis='y', labelcolor='#2563eb')
    ax.set_ylim(-0.05, 1.05)
    
    ax2 = ax.twinx()
    line2, = ax2.plot(Cs, alphas, 's--', color='#dc2626', linewidth=2, markersize=6,
                      label='⟨|α|⟩')
    ax2.set_ylabel('⟨|α|⟩  (mean polarization)', fontsize=11, color='#dc2626')
    ax2.tick_params(axis='y', labelcolor='#dc2626')
    ax2.set_ylim(-0.05, 1.15)
    
    ax.legend([line1, line2], ['f_S (proto-spatial)', '⟨|α|⟩ (polarization)'],
              loc='center right', fontsize=9)
    ax.set_title('A. Dose-response: C_protect', fontweight='bold', fontsize=12)
    
    # Annotate threshold region
    ax.axvspan(0.1, 0.2, alpha=0.08, color='green')
    ax.text(0.15, 0.95, 'onset', ha='center', fontsize=8, color='green', 
            transform=ax.get_xaxis_transform())
    
    # ── Panel B: Triangle sweep ──
    ax = axes[1]
    
    # Use ⟨tri⟩ as x-axis instead of p_ws (more intuitive)
    tri_x = [r['tri_mean'] for r in tri_on]
    fs_on_y = [r['fs_mean'] for r in tri_on]
    fs_on_err = [r['fs_std'] for r in tri_on]
    fs_off_y = [r['fs_mean'] for r in tri_off]
    
    ax.errorbar(tri_x, fs_on_y, yerr=fs_on_err, fmt='o-', color='#2563eb',
                linewidth=2.5, markersize=7, capsize=4, label='Protect ON')
    ax.plot(tri_x, fs_off_y, 's--', color='#ef4444', linewidth=2, markersize=6,
            label='Protect OFF')
    
    ax.set_xlabel('⟨tri⟩  (mean triangles per edge)', fontsize=12)
    ax.set_ylabel('f_S  (fraction proto-spatial)', fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=10)
    ax.set_title('B. Triangle level × Protector', fontweight='bold', fontsize=12)
    
    # Annotate
    ax.annotate('Triangles carry the signal,\nprotector translates it',
                xy=(2.5, 0.7), fontsize=9, style='italic', color='#1e40af',
                ha='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#dbeafe', edgecolor='#93c5fd'))
    
    plt.tight_layout()
    plt.savefig('/home/claude/fig_reviewer_checklist.png', dpi=180, bbox_inches='tight',
                facecolor='white')
    print('\nFigure saved → fig_reviewer_checklist.png')


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    t0 = time.time()
    
    dose = test_dose_response()
    tri_on, tri_off = test_triangle_sweep()
    
    make_figure(dose, tri_on, tri_off)
    
    # Save all data
    with open('/home/claude/reviewer_checklist.json', 'w') as f:
        json.dump({'dose_response': dose, 'tri_sweep_on': tri_on, 'tri_sweep_off': tri_off}, 
                  f, default=float)
    
    print(f'\nTotal time: {time.time() - t0:.0f}s')
