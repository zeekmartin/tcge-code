#!/usr/bin/env python3
"""
Figure 2 — Dose-response & Triangle sweep (polished)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json

with open('/home/claude/reviewer_checklist.json') as f:
    rc = json.load(f)

dose = rc['dose_response']
tri_on = rc['tri_sweep_on']
tri_off = rc['tri_sweep_off']

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle('Causal Controls: Dose-Response and Triangle Sweep',
             fontsize=13, fontweight='bold', y=1.01)

# ── A: Dose-response ──
ax = axes[0]
Cs = [r['C'] for r in dose]
fs = [r['fs_mean'] for r in dose]
fs_err = [r['fs_std']/np.sqrt(15) for r in dose]  # SEM
alphas = [r['alpha_mean'] for r in dose]

line1 = ax.errorbar(Cs, fs, yerr=fs_err, fmt='o-', color='#2563eb', linewidth=2.5,
                    markersize=7, capsize=3, label='f_S (proto-spatial)')
ax.fill_between(Cs, [f-e for f,e in zip(fs,fs_err)], [f+e for f,e in zip(fs,fs_err)],
                alpha=0.15, color='#2563eb')
ax.set_xlabel('C_protect  (protector strength)', fontsize=11)
ax.set_ylabel('f_S  (fraction proto-spatial)', fontsize=11, color='#2563eb')
ax.tick_params(axis='y', labelcolor='#2563eb')
ax.set_ylim(-0.05, 1.08)

ax2 = ax.twinx()
line2, = ax2.plot(Cs, alphas, 's--', color='#dc2626', linewidth=2, markersize=6,
                  label='⟨|α|⟩')
ax2.set_ylabel('⟨|α|⟩  (mean polarization)', fontsize=11, color='#dc2626')
ax2.tick_params(axis='y', labelcolor='#dc2626')
ax2.set_ylim(-0.05, 1.15)

# Threshold annotation
ax.axvline(0.175, color='#16a34a', linestyle=':', linewidth=1.5, alpha=0.7)
ax.annotate('onset\nC_crit ≈ (A−β)/⟨tri⟩', xy=(0.175, 0.10), xytext=(0.45, 0.15),
            fontsize=8.5, color='#16a34a', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#16a34a', lw=1.2))

ax.legend([line1, line2], ['f_S (proto-spatial)', '⟨|α|⟩ (polarization)'],
          loc='center right', fontsize=9)
ax.set_title('A. Dose-response: protector strength', fontweight='bold', fontsize=12)

# ── B: Triangle sweep ──
ax = axes[1]
tri_x = [r['tri_mean'] for r in tri_on]
fs_on = [r['fs_mean'] for r in tri_on]
fs_on_sem = [r['fs_std']/np.sqrt(15) for r in tri_on]
fs_off = [r['fs_mean'] for r in tri_off]

ax.errorbar(tri_x, fs_on, yerr=fs_on_sem, fmt='o-', color='#2563eb',
            linewidth=2.5, markersize=7, capsize=4, label='Protector ON (C=0.5)')
ax.plot(tri_x, fs_off, 's--', color='#ef4444', linewidth=2, markersize=6,
        label='Protector OFF (C=0)', zorder=5)

ax.set_xlabel('⟨tri⟩  (mean triangles per edge)', fontsize=11)
ax.set_ylabel('f_S  (fraction proto-spatial)', fontsize=11)
ax.set_ylim(-0.05, 1.08)
ax.legend(fontsize=10, loc='center left')
ax.set_title('B. Triangle level × protector', fontweight='bold', fontsize=12)

# Annotation
ax.annotate('Triangles carry the signal;\nprotector transduces it\ninto polarization heterogeneity',
            xy=(2.2, 0.6), fontsize=8.5, style='italic', color='#1e40af', ha='center',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#eff6ff', edgecolor='#93c5fd'))

# Flat-line emphasis
ax.annotate('OFF = flat at 0\nregardless of ⟨tri⟩',
            xy=(2.5, 0.0), xytext=(3.5, 0.15), fontsize=8, color='#ef4444',
            arrowprops=dict(arrowstyle='->', color='#ef4444', lw=1.2))

fig.text(0.5, -0.02,
         'Watts-Strogatz graphs (N=80, k=8), p_ws swept 0→0.9 | n=15 seeds/point, bars=±SEM',
         ha='center', fontsize=8.5, color='gray', style='italic')

plt.tight_layout(rect=[0, 0.02, 1, 0.97])
plt.savefig('/home/claude/fig2_checklist_polished.png', dpi=200, bbox_inches='tight', facecolor='white')
print('saved fig2')
