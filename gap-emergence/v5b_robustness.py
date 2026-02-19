#!/usr/bin/env python3
"""
TCGE — GAP-Emergence v5b : Tests de consolidation

Test A : Le biphasage persiste-t-il si on remplace tri(e) par d'autres 
         mesures de cohésion locale ? (Jaccard, edge clustering, k-truss, 4-cycles)
         → Si oui : "any local cohesion measure yields the same phase separation"

Test B : Reformulation TCGE-native : R(e) = nombre de compatibilités détruites
         si on polarise e. Le protecteur devient un coût de rupture de cohérence.
         → Transforme "triangles" en "compatibility destruction"

Test C : Scaling N = 30..400. Le biphasage converge-t-il vers un plateau ?
         → Si Δ se stabilise : effet structurel. Si Δ→0 : effet de taille finie.

Author: Zeek (TCGE project) — AI-assisted implementation
Date: 2026-02-19
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import time

# =============================================================================
# GRAPHE + MESURES DE COHÉSION
# =============================================================================

class ConstraintGraph:
    def __init__(self, n_nodes, edge_prob=0.3, seed=None):
        if seed is not None:
            np.random.seed(seed)
        self.n = n_nodes
        self.edges = []
        self.adj = defaultdict(set)
        for i in range(n_nodes):
            for j in range(i + 1, n_nodes):
                if np.random.random() < edge_prob:
                    self.edges.append((i, j))
                    self.adj[i].add(j)
                    self.adj[j].add(i)
        self.n_edges = len(self.edges)
        self.degree = np.array([len(self.adj[i]) for i in range(n_nodes)])
        self.N = self.degree.astype(float)
        self._compute_all_cohesion_measures()

    def _compute_all_cohesion_measures(self):
        """Compute multiple cohesion measures for each edge."""
        n = self.n_edges
        self.measures = {}
        
        # 1. Triangles
        tri = np.zeros(n, dtype=int)
        for idx, (i, j) in enumerate(self.edges):
            tri[idx] = len(self.adj[i] & self.adj[j])
        self.measures['triangles'] = tri.astype(float)
        
        # 2. Jaccard similarity
        jaccard = np.zeros(n)
        for idx, (i, j) in enumerate(self.edges):
            union = len(self.adj[i] | self.adj[j])
            inter = len(self.adj[i] & self.adj[j])
            jaccard[idx] = inter / max(union, 1)
        self.measures['jaccard'] = jaccard
        
        # 3. Edge clustering coefficient
        # ecc(e) = tri(e) / (min(deg(i), deg(j)) - 1)
        ecc = np.zeros(n)
        for idx, (i, j) in enumerate(self.edges):
            min_deg = min(len(self.adj[i]), len(self.adj[j]))
            if min_deg > 1:
                ecc[idx] = tri[idx] / (min_deg - 1)
        self.measures['edge_clustering'] = ecc
        
        # 4. Quadrangles (4-cycles) through edge
        quad = np.zeros(n, dtype=int)
        for idx, (i, j) in enumerate(self.edges):
            # Count paths i-k-l-j where k≠j, l≠i, k-l edge exists
            for k in self.adj[i]:
                if k == j:
                    continue
                for l in self.adj[j]:
                    if l == i or l == k:
                        continue
                    if l in self.adj[k]:
                        quad[idx] += 1
            quad[idx] //= 2  # each quad counted twice
        self.measures['quadrangles'] = quad.astype(float)
        
        # 5. k-truss level (simplified: max k such that edge is in k-truss)
        # k-truss: maximal subgraph where every edge is in ≥(k-2) triangles
        # Simplified: truss(e) = tri(e) + 2
        self.measures['truss'] = tri.astype(float) + 2
        
        # 6. TCGE-native: R(e) = compatibility destruction cost
        # R(e) = number of "coherent triplets" broken if e is polarized
        # A coherent triplet (i,j,k) means i-j, j-k, i-k all exist
        # Polarizing i-j affects all triplets containing i-j
        # R(e) = tri(e) (each triangle = one coherent triplet broken)
        # But we can extend: also count paths of length 2 through e
        # that would lose coherence
        compat_destroy = np.zeros(n)
        for idx, (i, j) in enumerate(self.edges):
            # Direct: triangles destroyed
            n_tri = tri[idx]
            # Indirect: 2-hop coherence lost
            # Paths i-j-k and k-i-j that rely on i-j being symmetric
            n_2hop = len(self.adj[i] - {j}) + len(self.adj[j] - {i})
            # Weighted: triangles count more (tighter coupling)
            compat_destroy[idx] = 2.0 * n_tri + 0.5 * n_2hop
        self.measures['compat_destroy'] = compat_destroy


def classify_edges(measure_values):
    """Split edges into high/low by median."""
    median = np.median(measure_values)
    # Handle case where median = min (many ties)
    if median == np.min(measure_values):
        # Use mean instead
        median = np.mean(measure_values)
    high = [i for i, v in enumerate(measure_values) if v > median]
    low = [i for i, v in enumerate(measure_values) if v <= median]
    return high, low, median


# =============================================================================
# OPTIMIZATION (reused from v5)
# =============================================================================

def optimize(graph, measure_name, A=1.0, C_protect=0.5, gamma=0.5,
             beta_reg=0.01, W=1.0, lr=0.005, n_steps=3000, seed=None):
    """
    Minimize E = Σ_e [A·(W/2)²·(1-α²) + C·cohesion(e)·α² + γ·W·α·ΔN] + β·Σα²
    """
    if seed is not None:
        np.random.seed(seed)
    
    n_e = graph.n_edges
    alpha_e = np.random.randn(n_e) * 0.01
    
    delta_N = np.array([graph.N[i] - graph.N[j] for i, j in graph.edges])
    cohesion = graph.measures[measure_name].copy()
    
    # Normalize cohesion to [0, max_tri_equivalent] for fair comparison
    if cohesion.max() > 0:
        cohesion = cohesion / cohesion.max() * graph.measures['triangles'].max()
    
    high_idx, low_idx, threshold = classify_edges(graph.measures[measure_name])
    
    for step in range(n_steps):
        grad = (2 * alpha_e * (-A * (W/2)**2 + C_protect * cohesion + beta_reg)
                + gamma * W * delta_N)
        alpha_e = alpha_e - lr * grad
        alpha_e = np.clip(alpha_e, -0.999, 0.999)
    
    # Results
    a_high = np.mean(np.abs(alpha_e[high_idx])) if high_idx else 0
    a_low = np.mean(np.abs(alpha_e[low_idx])) if low_idx else 0
    biphasage = a_low - a_high
    
    # Coherence
    coherence = 0
    n_counted = 0
    for idx in range(n_e):
        i, j = graph.edges[idx]
        dN = graph.N[i] - graph.N[j]
        if dN != 0:
            if np.sign(alpha_e[idx]) == -np.sign(dN):
                coherence += 1
            n_counted += 1
    coherence = coherence / max(n_counted, 1)
    
    return {
        'alpha_low': a_low, 'alpha_high': a_high,
        'biphasage': biphasage, 'coherence': coherence,
        'n_high': len(high_idx), 'n_low': len(low_idx),
        'alpha_e': alpha_e
    }


# =============================================================================
# TEST A : ROBUSTESSE AU CHOIX DE MESURE
# =============================================================================

def test_A_robustness(n_nodes=60, edge_prob=0.25, n_trials=25):
    print(f"\n{'═'*70}")
    print("  TEST A — ROBUSTESSE AU CHOIX DE MESURE DE COHÉSION")
    print(f"  Question: le biphasage dépend-il du choix 'triangles' ?")
    print(f"{'═'*70}\n")
    
    measures = ['triangles', 'jaccard', 'edge_clustering', 
                'quadrangles', 'truss', 'compat_destroy']
    
    all_results = {m: [] for m in measures}
    
    for trial in range(n_trials):
        graph = ConstraintGraph(n_nodes, edge_prob, seed=1000*trial + 42)
        
        for measure in measures:
            result = optimize(graph, measure, A=1.0, C_protect=0.5, 
                            gamma=0.5, seed=2000*trial + 7)
            all_results[measure].append(result)
    
    # Summary
    print(f"  {'Mesure':<20} {'|α|_low':<10} {'|α|_high':<10} "
          f"{'Biphasage':<12} {'Cohérence':<10} {'Status':<6}")
    print(f"  {'-'*70}")
    
    summary = {}
    for measure in measures:
        results = all_results[measure]
        bi = np.mean([r['biphasage'] for r in results])
        al = np.mean([r['alpha_low'] for r in results])
        ah = np.mean([r['alpha_high'] for r in results])
        coh = np.mean([r['coherence'] for r in results])
        std_bi = np.std([r['biphasage'] for r in results])
        status = '✅' if bi >= 0.2 else ('⚠️' if bi >= 0.1 else '❌')
        
        print(f"  {measure:<20} {al:<10.4f} {ah:<10.4f} "
              f"{bi:<12.4f} {coh:<10.3f} {status}")
        
        summary[measure] = {
            'biphasage': bi, 'std': std_bi,
            'alpha_low': al, 'alpha_high': ah, 'coherence': coh
        }
    
    # Verdict
    all_pass = all(summary[m]['biphasage'] >= 0.15 for m in measures)
    most_pass = sum(1 for m in measures if summary[m]['biphasage'] >= 0.15) >= 4
    
    print(f"\n  ┌────────────────────────────────────────────────────────┐")
    if all_pass:
        print(f"  │  ✅ TOUTES les mesures donnent un biphasage ≥ 0.15   │")
        print(f"  │  → 'Any local cohesion measure yields phase separation' │")
    elif most_pass:
        print(f"  │  ⚠️  La majorité des mesures donnent un biphasage     │")
        print(f"  │  → Robuste mais pas universel                          │")
    else:
        print(f"  │  ❌ Le biphasage dépend fortement de la mesure         │")
    print(f"  └────────────────────────────────────────────────────────┘")
    
    return summary, all_results


# =============================================================================
# TEST B : REFORMULATION TCGE-NATIVE
# =============================================================================

def test_B_tcge_native(n_nodes=60, edge_prob=0.25, n_trials=25):
    print(f"\n{'═'*70}")
    print("  TEST B — REFORMULATION TCGE-NATIVE")
    print(f"  Question: R(e) = 'compatibilités détruites' donne-t-il le même biphasage ?")
    print(f"{'═'*70}\n")
    
    # Compare: triangles (ad hoc) vs compat_destroy (TCGE-native)
    results_tri = []
    results_native = []
    
    for trial in range(n_trials):
        graph = ConstraintGraph(n_nodes, edge_prob, seed=1000*trial + 42)
        
        r_tri = optimize(graph, 'triangles', A=1.0, C_protect=0.5,
                        gamma=0.5, seed=2000*trial + 7)
        r_native = optimize(graph, 'compat_destroy', A=1.0, C_protect=0.5,
                           gamma=0.5, seed=2000*trial + 7)
        
        results_tri.append(r_tri)
        results_native.append(r_native)
    
    bi_tri = [r['biphasage'] for r in results_tri]
    bi_nat = [r['biphasage'] for r in results_native]
    
    print(f"  Triangles (ad hoc)        : Δ = {np.mean(bi_tri):.4f} ± {np.std(bi_tri):.4f}")
    print(f"  Compat. destroy (native)  : Δ = {np.mean(bi_nat):.4f} ± {np.std(bi_nat):.4f}")
    print(f"  Corrélation trial-by-trial: {np.corrcoef(bi_tri, bi_nat)[0,1]:.3f}")
    
    # Detailed comparison per trial
    print(f"\n  {'Trial':<8} {'Δ_tri':<10} {'Δ_native':<10} {'Diff':<10}")
    print(f"  {'-'*40}")
    for t in range(min(8, n_trials)):
        print(f"  {t+1:<8} {bi_tri[t]:<10.4f} {bi_nat[t]:<10.4f} "
              f"{bi_nat[t]-bi_tri[t]:+.4f}")
    
    equivalent = abs(np.mean(bi_nat) - np.mean(bi_tri)) < 0.05
    
    print(f"\n  ┌────────────────────────────────────────────────────────┐")
    if equivalent:
        print(f"  │  ✅ Le coût 'compatibility destruction' reproduit     │")
        print(f"  │  le même biphasage que les triangles.                  │")
        print(f"  │  → Le protecteur N'EST PAS ad hoc : c'est le coût     │")
        print(f"  │  de rupture de cohérence locale.                       │")
    else:
        print(f"  │  ⚠️ Différence significative entre les deux mesures    │")
    print(f"  └────────────────────────────────────────────────────────┘")
    
    return bi_tri, bi_nat


# =============================================================================
# TEST C : FINITE-SIZE SCALING
# =============================================================================

def test_C_scaling(C_protect=0.5, n_trials=20):
    print(f"\n{'═'*70}")
    print("  TEST C — FINITE-SIZE SCALING")
    print(f"  Question: Δ(N) converge-t-il vers un plateau ou → 0 ?")
    print(f"{'═'*70}\n")
    
    # Adjusted edge_prob to keep mean degree roughly constant
    # mean_degree ≈ (N-1) * edge_prob ≈ 15
    sizes = [30, 50, 70, 100, 150, 200, 300]
    target_mean_deg = 15
    
    scale_data = []
    
    for N in sizes:
        p = min(target_mean_deg / (N - 1), 0.8)
        
        biphasages = []
        a_lows = []
        a_highs = []
        cohs = []
        tri_vars = []
        
        for trial in range(n_trials):
            graph = ConstraintGraph(N, p, seed=1000*trial + 42)
            result = optimize(graph, 'triangles', A=1.0, C_protect=C_protect,
                            gamma=0.5, n_steps=3000, seed=2000*trial + 7)
            biphasages.append(result['biphasage'])
            a_lows.append(result['alpha_low'])
            a_highs.append(result['alpha_high'])
            cohs.append(result['coherence'])
            tri_vars.append(np.var(graph.measures['triangles']))
        
        mean_bi = np.mean(biphasages)
        std_bi = np.std(biphasages)
        mean_tri_var = np.mean(tri_vars)
        
        # Normalized biphasage: Δ / σ(tri) 
        norm_bi = mean_bi / np.sqrt(mean_tri_var) if mean_tri_var > 0 else 0
        
        scale_data.append({
            'N': N, 'p': p, 
            'biphasage': mean_bi, 'std': std_bi,
            'alpha_low': np.mean(a_lows), 'alpha_high': np.mean(a_highs),
            'coherence': np.mean(cohs),
            'tri_variance': mean_tri_var,
            'norm_biphasage': norm_bi,
            'pct_02': 100 * np.mean([b > 0.2 for b in biphasages]),
            'pct_01': 100 * np.mean([b > 0.1 for b in biphasages])
        })
        
        status = '✅' if mean_bi >= 0.2 else ('⚠️' if mean_bi >= 0.1 else '❌')
        print(f"  N={N:<4} p={p:.3f}: Δ={mean_bi:.4f}±{std_bi:.4f}  "
              f"|α|_L={np.mean(a_lows):.3f} |α|_H={np.mean(a_highs):.3f}  "
              f"coh={np.mean(cohs):.2f}  {status}")
    
    # Fit: does Δ decay as 1/sqrt(N), 1/N, or plateau?
    Ns = np.array([s['N'] for s in scale_data])
    Ds = np.array([s['biphasage'] for s in scale_data])
    
    # Log-log fit
    valid = Ds > 0
    if np.sum(valid) >= 3:
        log_N = np.log(Ns[valid])
        log_D = np.log(Ds[valid])
        slope, intercept = np.polyfit(log_N, log_D, 1)
        print(f"\n  Log-log fit: Δ ∝ N^{slope:.3f}")
        print(f"    slope = -0.5 → 1/√N (taille finie)")
        print(f"    slope = -1.0 → 1/N (disparaît)")
        print(f"    slope ≈ 0   → plateau (structurel)")
    else:
        slope = 0
    
    # Verdict
    print(f"\n  ┌────────────────────────────────────────────────────────┐")
    if slope > -0.3:
        print(f"  │  ✅ Slope = {slope:.3f} ≈ 0 → PLATEAU                │")
        print(f"  │  Le biphasage est un effet structurel, pas de taille. │")
    elif slope > -0.6:
        print(f"  │  ⚠️ Slope = {slope:.3f} → décroissance lente          │")
        print(f"  │  Le biphasage persiste mais s'atténue.                 │")
    else:
        print(f"  │  ❌ Slope = {slope:.3f} → le biphasage disparaît       │")
        print(f"  │  C'est un effet de taille finie.                       │")
    print(f"  └────────────────────────────────────────────────────────┘")
    
    return scale_data, slope


# =============================================================================
# VISUALISATION
# =============================================================================

def plot_all(test_a_summary, test_b_tri, test_b_nat, test_c_data, slope,
             filename="/home/claude/tcge_gap_emergence_v5b.png"):
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    fig.suptitle("TCGE — GAP-Emergence v5b : Tests de Consolidation\n"
                 "A: Robustesse mesure | B: TCGE-native | C: Scaling",
                 fontsize=13, fontweight='bold', y=0.99)
    
    # ── TEST A : Bar chart des biphasages par mesure ──
    ax = axes[0, 0]
    measures = list(test_a_summary.keys())
    bis = [test_a_summary[m]['biphasage'] for m in measures]
    stds = [test_a_summary[m]['std'] for m in measures]
    colors = ['#2ecc71' if b >= 0.2 else '#f39c12' if b >= 0.1 else '#e74c3c' 
              for b in bis]
    bars = ax.bar(range(len(measures)), bis, yerr=stds, capsize=4, 
                  color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax.axhline(0.2, color='green', linestyle='--', alpha=0.7, label='Seuil 0.2')
    ax.axhline(0.1, color='orange', linestyle='--', alpha=0.5, label='Seuil 0.1')
    ax.set_xticks(range(len(measures)))
    ax.set_xticklabels([m.replace('_', '\n') for m in measures], fontsize=7)
    ax.set_ylabel('Biphasage Δ')
    ax.set_title('A. Robustesse au choix de mesure', fontweight='bold')
    ax.legend(fontsize=7)
    
    # ── TEST A : α_low vs α_high par mesure ──
    ax = axes[0, 1]
    for i, m in enumerate(measures):
        al = test_a_summary[m]['alpha_low']
        ah = test_a_summary[m]['alpha_high']
        ax.scatter(ah, al, s=100, zorder=5, edgecolors='black', linewidth=0.5)
        ax.annotate(m[:6], (ah, al), fontsize=6, ha='center', va='bottom')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='Pas de biphasage')
    ax.set_xlabel('|α|_high (proto-S)')
    ax.set_ylabel('|α|_low (proto-T)')
    ax.set_title('A. Séparation par mesure', fontweight='bold')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.legend(fontsize=7)
    
    # ── TEST B : Scatter tri vs native ──
    ax = axes[0, 2]
    ax.scatter(test_b_tri, test_b_nat, s=50, alpha=0.7, c='steelblue',
              edgecolors='black', linewidth=0.3)
    lim = max(max(test_b_tri), max(test_b_nat)) * 1.1
    ax.plot([0, lim], [0, lim], 'k--', alpha=0.3, label='y=x')
    ax.set_xlabel('Δ (triangles)')
    ax.set_ylabel('Δ (compat. destroy)')
    ax.set_title('B. TCGE-native vs triangles', fontweight='bold')
    ax.legend(fontsize=8)
    corr = np.corrcoef(test_b_tri, test_b_nat)[0, 1]
    ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
           fontsize=10, va='top', fontweight='bold')
    
    # ── TEST C : Biphasage vs N ──
    ax = axes[1, 0]
    Ns = [s['N'] for s in test_c_data]
    Ds = [s['biphasage'] for s in test_c_data]
    Ds_std = [s['std'] for s in test_c_data]
    ax.errorbar(Ns, Ds, yerr=Ds_std, fmt='o-', color='steelblue', 
                linewidth=2, capsize=5, markersize=8)
    ax.axhline(0.2, color='green', linestyle='--', alpha=0.5)
    ax.axhline(0.1, color='orange', linestyle='--', alpha=0.5)
    ax.set_xlabel('N (nodes)')
    ax.set_ylabel('Biphasage Δ')
    ax.set_title(f'C. Scaling (slope={slope:.3f})', fontweight='bold')
    
    # Power law fit overlay
    if len(Ns) >= 3:
        N_fit = np.linspace(min(Ns), max(Ns), 100)
        D_fit = np.exp(np.log(Ds[0]) + slope * (np.log(N_fit) - np.log(Ns[0])))
        ax.plot(N_fit, D_fit, 'r--', alpha=0.5, label=f'N^{slope:.2f}')
        ax.legend(fontsize=8)
    
    # ── TEST C : |α|_low et |α|_high vs N ──
    ax = axes[1, 1]
    als = [s['alpha_low'] for s in test_c_data]
    ahs = [s['alpha_high'] for s in test_c_data]
    ax.plot(Ns, als, 'o-', color='tab:red', linewidth=2, markersize=8, label='|α|_low (T)')
    ax.plot(Ns, ahs, 's-', color='tab:blue', linewidth=2, markersize=8, label='|α|_high (S)')
    ax.fill_between(Ns, ahs, als, alpha=0.15, color='green')
    ax.set_xlabel('N (nodes)')
    ax.set_ylabel('|α|')
    ax.set_title('C. Séparation vs taille', fontweight='bold')
    ax.legend(fontsize=8)
    
    # ── TEST C : Log-log ──
    ax = axes[1, 2]
    valid = np.array(Ds) > 0
    Ns_v = np.array(Ns)[valid]
    Ds_v = np.array(Ds)[valid]
    ax.plot(np.log10(Ns_v), np.log10(Ds_v), 'o-', color='steelblue',
           linewidth=2, markersize=8)
    # Reference lines
    N_ref = np.linspace(min(np.log10(Ns_v)), max(np.log10(Ns_v)), 50)
    ax.plot(N_ref, np.log10(Ds_v[0]) - 0.5*(N_ref - np.log10(Ns_v[0])),
           'r--', alpha=0.4, label='∝ 1/√N')
    ax.plot(N_ref, np.log10(Ds_v[0]) + 0*(N_ref - np.log10(Ns_v[0])),
           'g--', alpha=0.4, label='plateau')
    ax.set_xlabel('log₁₀(N)')
    ax.set_ylabel('log₁₀(Δ)')
    ax.set_title('C. Log-log scaling', fontweight='bold')
    ax.legend(fontsize=8)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"\nFigure sauvegardée : {filename}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    t0 = time.time()
    
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  TCGE — GAP-Emergence v5b : Tests de Consolidation        ║")
    print("║  A: Robustesse | B: TCGE-native | C: Scaling              ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    # ── TEST A ──
    test_a_summary, test_a_all = test_A_robustness(
        n_nodes=60, edge_prob=0.25, n_trials=25)
    
    # ── TEST B ──
    test_b_tri, test_b_nat = test_B_tcge_native(
        n_nodes=60, edge_prob=0.25, n_trials=25)
    
    # ── TEST C ──
    test_c_data, slope = test_C_scaling(C_protect=0.5, n_trials=20)
    
    # ── VISUALISATION ──
    plot_all(test_a_summary, test_b_tri, test_b_nat, test_c_data, slope)
    
    # ── CONCLUSION ──
    elapsed = time.time() - t0
    
    print(f"\n\n{'═'*70}")
    print("  CONCLUSION FINALE — GAP-Emergence v5b")
    print(f"{'═'*70}")
    
    # Test A verdict
    n_pass_a = sum(1 for m in test_a_summary 
                   if test_a_summary[m]['biphasage'] >= 0.15)
    n_total_a = len(test_a_summary)
    
    # Test B verdict
    corr_b = np.corrcoef(test_b_tri, test_b_nat)[0, 1]
    equiv_b = abs(np.mean(test_b_nat) - np.mean(test_b_tri)) < 0.05
    
    # Test C verdict
    plateau_c = slope > -0.3
    
    print(f"""
  Temps total : {elapsed:.0f}s

  ┌──────────┬────────────────────────────────┬──────────┐
  │ Test     │ Résultat                       │ Verdict  │
  ├──────────┼────────────────────────────────┼──────────┤
  │ A robust │ {n_pass_a}/{n_total_a} mesures Δ≥0.15            │ {'✅' if n_pass_a >= 4 else '❌'}        │
  │ B native │ corr={corr_b:.3f}, equiv={'oui' if equiv_b else 'non'}            │ {'✅' if equiv_b else '❌'}        │
  │ C scale  │ slope={slope:.3f}                    │ {'✅' if plateau_c else '❌'}        │
  └──────────┴────────────────────────────────┴──────────┘

  INTERPRÉTATION :
""")
    
    if n_pass_a >= 4 and equiv_b and plateau_c:
        print("""  ✅✅✅ LES TROIS TESTS PASSENT.
  
  Le biphasage T/S est :
    • Robuste au choix de mesure (pas spécifique aux triangles)
    • Reformulable en coût TCGE-native (rupture de compatibilité)
    • Stable en taille (pas un artefact de taille finie)
  
  → GAP-Emergence(Signature) peut être considéré comme 
    PARTIELLEMENT FERMÉ avec justification solide.
  
  Formulation défendable :
  "Phase separation between polarized (temporal) and symmetric (spatial)
   edges emerges from competition between symmetry-breaking product cost
   and local coherence protection. The protection term is not ad hoc:
   it measures compatibility destruction cost. The separation is robust
   across cohesion metrics and graph sizes."
""")
    else:
        issues = []
        if n_pass_a < 4: issues.append("robustesse mesure insuffisante")
        if not equiv_b: issues.append("reformulation native diverge")
        if not plateau_c: issues.append("effet de taille finie")
        print(f"  ⚠️ Points à consolider : {', '.join(issues)}")
    
    print(f"{'═'*70}")
