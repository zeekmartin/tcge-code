#!/usr/bin/env python3
"""
Figure 3 — Mechanistic Ablation (polished)
With P(e) = C_protect * tri(e) as intermediate variable.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx


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
    for it in range(n_iter):
        eff = -A + beta_reg + C_protect * tri_arr / (W**2)
        grad = 2 * eff * w
        w -= lr * grad
        if it < n_iter * 0.7:
            noise = rng.normal(0, 0.02 * (1 - it / n_iter), M)
            w += noise
        w = np.clip(w, -W, W)
        if it == n_iter // 2:
            lr *= 0.3
    return np.abs(w) / W, tri_arr


def run_pool(graph_type, protect, N=80, n_seeds=5):
    all_tri, all_alpha = [], []
    for seed in range(n_seeds):
        if graph_type == 'WS':
            G = nx.watts_strogatz_graph(N, 8, 0.1, seed=1000+seed)
        else:
            G = nx.erdos_renyi_graph(N, 8/(N-1), seed=1000+seed)
            if not nx.is_connected(G):
                comps = list(nx.connected_components(G))
                for c in comps[1:]:
                    G.add_edge(list(c)[0], list(comps[0])[0])
        C = 0.5 if protect else 0.0
        alpha, tri_arr = optimize(G, C_protect=C, seed=3000+seed)
        all_tri.extend(tri_arr)
        all_alpha.extend(alpha)
    return np.array(all_tri), np.array(all_alpha)


# ═══════════════════════════════════════════
# FIGURE
# ═══════════════════════════════════════════

fig, axes = plt.subplots(2, 2, figsize=(11, 9))
fig.suptitle('Mechanistic Ablation: P(e) = C·tri(e) as Protection Potential',
             fontsize=13, fontweight='bold', y=0.99)

configs = [
    ('WS', True,  'A. Tri↑ + Prot ON: signal transduced',  '#2563eb'),
    ('ER', True,  'B. Tri↓ + Prot ON: weak signal',         '#f59e0b'),
    ('WS', False, 'C. Tri↑ + Prot OFF: signal ignored',     '#ef4444'),
    ('ER', False, 'D. Tri↓ + Prot OFF: baseline',           '#6b7280'),
]

C_val = 0.5
A_val = 1.0
beta_val = 0.01
threshold = A_val - beta_val  # P(e) > threshold → edge protected

for idx, (gtype, protect, title, color) in enumerate(configs):
    ax = axes[idx // 2, idx % 2]
    tri_arr, alpha = run_pool(gtype, protect)
    
    # X-axis: P(e) = C * tri(e) (or 0 * tri(e) = 0 when protect OFF)
    C_used = C_val if protect else 0.0
    P_e = C_used * tri_arr
    
    # Jitter for visibility
    jitter = np.random.RandomState(42).uniform(-0.03, 0.03, len(P_e))
    ax.scatter(P_e + jitter, alpha, s=3, alpha=0.12, color=color, rasterized=True)
    
    # Bin means by integer tri
    max_tri = int(tri_arr.max())
    for t in range(max_tri + 1):
        mask = tri_arr == t
        if np.sum(mask) > 5:
            p_val = C_used * t
            ax.scatter(p_val, np.mean(alpha[mask]), s=130, color=color,
                      edgecolors='black', linewidth=1.5, zorder=10)
    
    # Threshold line
    if protect:
        ax.axvline(threshold, color='#16a34a', linestyle=':', linewidth=1.5, alpha=0.6)
        ax.text(threshold + 0.05, 0.55, 'P(e) = A−β\n(protection\nthreshold)',
                fontsize=7.5, color='#16a34a', va='center')
    
    # Correlation
    if np.std(P_e) > 0 and np.std(alpha) > 0:
        r = np.corrcoef(P_e, alpha)[0, 1]
        ax.text(0.95, 0.95, f'r = {r:.2f}', transform=ax.transAxes,
                fontsize=12, fontweight='bold', ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='#d1d5db'))
    else:
        ax.text(0.95, 0.95, 'P(e) = 0\n∀ edges', transform=ax.transAxes,
                fontsize=10, ha='right', va='top', color='#6b7280',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
    
    ax.text(0.05, 0.05, f'⟨|α|⟩ = {np.mean(alpha):.2f}', transform=ax.transAxes,
            fontsize=10, color=color, fontweight='bold')
    
    ax.set_xlabel('P(e) = C_protect · tri(e)', fontsize=10)
    ax.set_ylabel('|α(e)|  (polarization)', fontsize=10)
    ax.set_title(title, fontsize=10.5, fontweight='bold', color=color)
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(-0.15, max(C_val * max_tri, 0.5) + 0.3)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.25)

fig.text(0.5, 0.005,
         'P(e) = C·tri(e): local protection potential. Edge is protected (spatial) when P(e) > A−β ≈ 0.99.\n'
         'Large circles = bin means over 5 seeds (N=80). Lorentzian-like separation = negative gradient P(e) → |α(e)|.',
         ha='center', fontsize=8.5, color='gray', style='italic')

plt.tight_layout(rect=[0, 0.04, 1, 0.97])
plt.savefig('/home/claude/fig3_mechanism_polished.png', dpi=200, bbox_inches='tight', facecolor='white')
print('saved fig3')
