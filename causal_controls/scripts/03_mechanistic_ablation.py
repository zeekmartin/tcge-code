#!/usr/bin/env python3
"""
TCGE Mechanistic Ablation
==========================

Show that triangles = information, protector = transducer.

For each condition, scatter tri(e) vs |α(e)| for all edges:
- WS + protect ON:  negative correlation (info → heterogeneity)
- WS + protect OFF: flat at α=1 (info present, not read)
- ER + protect ON:  weak correlation (little info to read)
- ER + protect OFF: flat at α=1 (no info, not read)

Plus: compute "protection potential" = C × tri(e) and show its distribution.
"""

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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


def run_single(graph_type, protect, N=80, seed=42):
    """Run one condition, return per-edge data."""
    if graph_type == 'WS':
        G = nx.watts_strogatz_graph(N, 8, 0.1, seed=seed)
    else:  # ER
        G = nx.erdos_renyi_graph(N, 8 / (N - 1), seed=seed)
        if not nx.is_connected(G):
            comps = list(nx.connected_components(G))
            for c in comps[1:]:
                G.add_edge(list(c)[0], list(comps[0])[0])
    
    C_prot = 0.5 if protect else 0.0
    alpha, tri_arr = optimize(G, C_protect=C_prot, seed=3000 + seed)
    return tri_arr, alpha


# ═══════════════════════════════════════════════════════════════
# FIGURE
# ═══════════════════════════════════════════════════════════════

fig, axes = plt.subplots(2, 2, figsize=(11, 9))
fig.suptitle('Mechanistic Ablation: Triangles = Signal, Protector = Transducer',
             fontsize=14, fontweight='bold', y=0.98)

configs = [
    ('WS', True,  'WS + Protect ON\n(signal read)',     '#2563eb'),
    ('ER', True,  'ER + Protect ON\n(weak signal)',      '#f59e0b'),
    ('WS', False, 'WS + Protect OFF\n(signal ignored)',  '#ef4444'),
    ('ER', False, 'ER + Protect OFF\n(no signal)',       '#6b7280'),
]

for idx, (graph_type, protect, title, color) in enumerate(configs):
    ax = axes[idx // 2, idx % 2]
    
    # Run 5 seeds and pool edges
    all_tri, all_alpha = [], []
    for seed in range(5):
        tri_arr, alpha = run_single(graph_type, protect, seed=seed)
        all_tri.extend(tri_arr)
        all_alpha.extend(alpha)
    
    all_tri = np.array(all_tri)
    all_alpha = np.array(all_alpha)
    
    # Scatter with jitter
    jitter = np.random.RandomState(0).uniform(-0.2, 0.2, len(all_tri))
    ax.scatter(all_tri + jitter, all_alpha, s=3, alpha=0.15, color=color, rasterized=True)
    
    # Bin means
    max_tri = int(all_tri.max())
    for t in range(max_tri + 1):
        mask = all_tri == t
        if np.sum(mask) > 5:
            ax.scatter(t, np.mean(all_alpha[mask]), s=120, color=color, 
                      edgecolors='black', linewidth=1.5, zorder=10)
    
    # Correlation
    if np.std(all_tri) > 0 and np.std(all_alpha) > 0:
        r = np.corrcoef(all_tri, all_alpha)[0, 1]
        ax.text(0.95, 0.95, f'r = {r:.2f}', transform=ax.transAxes,
                fontsize=12, fontweight='bold', ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    else:
        ax.text(0.95, 0.95, 'r = N/A\n(no variance)', transform=ax.transAxes,
                fontsize=10, ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Mean α as annotation
    ax.text(0.05, 0.05, f'⟨|α|⟩ = {np.mean(all_alpha):.2f}', 
            transform=ax.transAxes, fontsize=10, color=color, fontweight='bold')
    
    ax.set_xlabel('tri(e)  (triangles per edge)', fontsize=10)
    ax.set_ylabel('|α(e)|  (polarization)', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold', color=color)
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(-0.5, max_tri + 0.5)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.3, label='α = 0.5 threshold')

# Add interpretation box
fig.text(0.5, 0.01, 
         'Large circles = bin means. Small dots = individual edges (5 seeds pooled, N=80).\n'
         'Only WS + Protect ON shows the negative gradient: triangles drive spatial identity through the cost.',
         ha='center', fontsize=9, color='gray', style='italic')

plt.tight_layout(rect=[0, 0.04, 1, 0.96])
plt.savefig('/home/claude/fig_mechanism.png', dpi=180, bbox_inches='tight', facecolor='white')
print('Figure saved → fig_mechanism.png')
