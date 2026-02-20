#!/usr/bin/env python3
"""
Phase B — Rewiring control figure (recreated from session data).
Shows the dual role of the triangle protector.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Data from dim_causal.py session
rewire = [0, 5, 10, 15, 20, 30, 50, 70, 90]
delta =  [0.369, 0.452, 0.502, 0.486, 0.473, 0.299, 0.211, 0.157, 0.139]
alpha_mean = [0.433, 0.549, 0.648, 0.729, 0.790, 0.876, 0.954, 0.978, 0.982]
tri_mean = [3.75, 2.78, 2.06, 1.52, 1.14, 0.64, 0.20, 0.08, 0.04]

fig, ax1 = plt.subplots(figsize=(9, 5.5))

ax1.axvspan(-1, 22, alpha=0.06, color='#22c55e', zorder=0)
ax1.axvspan(22, 55, alpha=0.06, color='#f59e0b', zorder=0)
ax1.axvspan(55, 95, alpha=0.06, color='#ef4444', zorder=0)

ax1.text(10, 0.03, 'protection\nactive', ha='center', fontsize=8, color='#166534', style='italic')
ax1.text(38, 0.03, 'Δ collapses', ha='center', fontsize=8, color='#92400e', style='italic')
ax1.text(75, 0.03, 'uniform\npolarization', ha='center', fontsize=8, color='#991b1b', style='italic')

line1, = ax1.plot(rewire, delta, 'o-', color='#2563eb', linewidth=2.5, markersize=7,
                  label='Δ (phase separation)', zorder=5)
ax1.set_xlabel('Rewiring fraction (%)', fontsize=12)
ax1.set_ylabel('Δ  (tri-based separation)', fontsize=11, color='#2563eb')
ax1.tick_params(axis='y', labelcolor='#2563eb')
ax1.set_ylim(0, 0.6)
ax1.set_xlim(-2, 93)

ax2 = ax1.twinx()
line2, = ax2.plot(rewire, alpha_mean, 's--', color='#dc2626', linewidth=2, markersize=6,
                  label='⟨|α|⟩ (polarization)', zorder=4)
ax2.set_ylabel('⟨|α|⟩  (mean polarization)', fontsize=11, color='#dc2626')
ax2.tick_params(axis='y', labelcolor='#dc2626')
ax2.set_ylim(0.3, 1.05)

ax3 = ax1.twinx()
ax3.spines['right'].set_position(('outward', 0))
ax3.set_yticks([])
bars = ax3.bar(rewire, tri_mean, width=3.5, alpha=0.15, color='#9ca3af', zorder=1, label='⟨tri⟩')
ax3.set_ylim(0, 8)

ax1.annotate('Δ peak\n(⟨tri⟩ ≈ 2)', xy=(10, 0.502), xytext=(25, 0.55), fontsize=9,
             arrowprops=dict(arrowstyle='->', color='#1e40af', lw=1.5),
             color='#1e40af', fontweight='bold')

lines = [line1, line2, bars]
labels = ['Δ (structured separation)', '⟨|α|⟩ (total polarization)', '⟨tri⟩ (triangles)']
ax1.legend(lines, labels, loc='center right', fontsize=9, framealpha=0.9, edgecolor='#d1d5db')
ax1.set_title('Triangle protector: structures separation, not polarization',
              fontsize=13, fontweight='bold', pad=12)

plt.tight_layout()
plt.savefig('/home/claude/tcge_package/phaseB_dimensionality/figures/fig_rewiring_control.png',
            dpi=180, bbox_inches='tight')
print('saved')
