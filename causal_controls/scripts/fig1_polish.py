#!/usr/bin/env python3
"""
Figure 1 — Double Dissociation (polished, reviewer-proof)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json

with open('/home/claude/dd_results.json') as f:
    data = json.load(f)

conditions = ['ER', 'C1', 'C2', 'REF']

# Compute means and SEMs
n = len(data['ER']['frac_spatial'])
stats = {}
for c in conditions:
    fs = np.array(data[c]['frac_spatial'])
    al = np.array(data[c]['mean_alpha'])
    mu = np.array(data[c]['mu1'])
    stats[c] = {
        'fs_m': np.mean(fs), 'fs_sem': np.std(fs)/np.sqrt(n),
        'al_m': np.mean(al), 'al_sem': np.std(al)/np.sqrt(n),
        'mu_m': np.mean(mu), 'mu_sem': np.std(mu)/np.sqrt(n),
        'tri_m': np.mean(data[c]['mean_tri']),
    }

# ═══════════════════════════════════════════════
# FIGURE: 1×3 layout (compact, slide-friendly)
# ═══════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(14, 5))
fig.suptitle('Double Dissociation: Triangles × Protector → Lorentzian-like Phase Separation',
             fontsize=13, fontweight='bold', y=1.01)

x = np.arange(4)
bar_w = 0.55
nice = ['Tri↑\nProt ON', 'Tri↓\nProt ON', 'Tri↑\nProt OFF', 'Tri↓\nProt OFF']
colors = ['#2563eb', '#f59e0b', '#ef4444', '#6b7280']

# ── A: f_S ──
ax = axes[0]
fs_vals = [stats[c]['fs_m'] for c in conditions]
fs_errs = [stats[c]['fs_sem'] for c in conditions]
bars = ax.bar(x, fs_vals, bar_w, yerr=fs_errs, capsize=4,
              color=colors, edgecolor='black', linewidth=0.7, alpha=0.85)
for i, (v, e) in enumerate(zip(fs_vals, fs_errs)):
    ax.text(i, v + e + 0.03, f'{v:.2f}', ha='center', va='bottom', 
            fontsize=11, fontweight='bold')
ax.set_ylabel('f_S  (fraction |α| < 0.5)', fontsize=11)
ax.set_title('A. Proto-spatial fraction', fontweight='bold', fontsize=12)
ax.set_xticks(x); ax.set_xticklabels(nice, fontsize=9)
ax.set_ylim(0, 1.1)

# ── B: ⟨|α|⟩ ──
ax = axes[1]
al_vals = [stats[c]['al_m'] for c in conditions]
al_errs = [stats[c]['al_sem'] for c in conditions]
bars = ax.bar(x, al_vals, bar_w, yerr=al_errs, capsize=4,
              color=colors, edgecolor='black', linewidth=0.7, alpha=0.85)
for i, (v, e) in enumerate(zip(al_vals, al_errs)):
    ax.text(i, v + e + 0.03, f'{v:.2f}', ha='center', va='bottom',
            fontsize=11, fontweight='bold')
ax.set_ylabel('⟨|α|⟩  (mean polarization)', fontsize=11)
ax.set_title('B. Polarization level', fontweight='bold', fontsize=12)
ax.set_xticks(x); ax.set_xticklabels(nice, fontsize=9)
ax.set_ylim(0, 1.2)

# ── C: 2×2 Heatmap ──
ax = axes[2]
mat = np.array([[stats['ER']['fs_m'], stats['C2']['fs_m']],
                [stats['C1']['fs_m'], stats['REF']['fs_m']]])
im = ax.imshow(mat, cmap='RdYlBu', vmin=0, vmax=1, aspect='auto')
ax.set_xticks([0, 1]); ax.set_xticklabels(['Protector ON', 'Protector OFF'], fontsize=10)
ax.set_yticks([0, 1]); ax.set_yticklabels(['Triangles ↑\n(WS, ⟨tri⟩≈3.3)', 
                                            'Triangles ↓\n(ER, ⟨tri⟩≈0.8)'], fontsize=9)
ax.set_title('C. Interaction (f_S)', fontweight='bold', fontsize=12)
for i in range(2):
    for j in range(2):
        color = 'white' if mat[i, j] > 0.5 else 'black'
        ax.text(j, i, f'{mat[i, j]:.2f}', ha='center', va='center',
                fontsize=18, fontweight='bold', color=color)
plt.colorbar(im, ax=ax, label='f_S', shrink=0.8, pad=0.08)

# ANOVA annotation
ax.text(0.5, -0.22, 'ANOVA interaction: F(1,96) = 4497, p < 10⁻⁶⁰\n'
        'η²(interaction) = 0.23, η²(residual) = 0.005',
        transform=ax.transAxes, ha='center', fontsize=8.5, color='#374151',
        style='italic',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f3f4f6', edgecolor='#d1d5db'))

fig.text(0.5, -0.02, 
         'WS: Watts-Strogatz (N=80, k=8, p=0.1) | ER: Erdős-Rényi (same ⟨k⟩) | '
         'n=25 seeds/condition, bars=±SEM',
         ha='center', fontsize=8.5, color='gray', style='italic')

plt.tight_layout(rect=[0, 0.02, 1, 0.97])
plt.savefig('/home/claude/fig1_dd_polished.png', dpi=200, bbox_inches='tight', facecolor='white')
print('saved fig1')
