# TCGE — Causal Controls for Lorentzian-like Phase Separation

**Status:** Core result. Reviewer-proof package.  
**Date:** February 2025  

## Summary

Double dissociation establishing that Lorentzian-like phase separation
requires BOTH local cohesive structure (triangles) AND a coherence-protection
term in the cost functional. Neither alone suffices.

## Structure

```
causal_controls/
├── text/
│   ├── results_1page.md           ← Main results section (paper-ready)
│   ├── derivation_protector.md    ← Why tri(e) is necessary (3 axioms)
│   └── gaps_theoretical.md        ← Honest gap analysis (4 gaps)
├── figures/
│   ├── fig1_double_dissociation   ← 2×2 factorial + heatmap
│   ├── fig2_dose_response         ← C_protect sweep + triangle sweep
│   ├── fig3_mechanism             ← Per-edge scatter (P(e) → |α|)
│   ├── fig4_robustness            ← Threshold / form / scale invariance
│   ├── fig5_bifurcation_potentials← V(w) for 6 penalty families
│   ├── fig5b_phase_diagram        ← Universal bifurcation diagram
│   └── fig6_dynamic_stability     ← Perturbation + relaxation + basins
├── scripts/
│   ├── 01_double_dissociation.py  ← 2×2 factorial (WS vs ER × ON vs OFF)
│   ├── 02_dose_response_...py     ← C_protect sweep + triangle sweep
│   ├── 03_mechanistic_ablation.py ← Per-edge P(e) vs |α| scatter
│   ├── 04_robustness_blindage.py  ← τ-sweep, form variation, scaling
│   ├── 05_bifurcation_theorem.py  ← Analytical proof + 6 families
│   ├── 06_dynamic_stability.py    ← Perturbation, relaxation, basins
│   └── fig[1-3]_polish.py        ← Polished figure generation
└── data/
    ├── dd_results.json            ← Raw double dissociation data
    └── reviewer_checklist.json    ← Dose-response + sweep data
```

## Key Results

| Test | Result | Implication |
|------|--------|-------------|
| Double dissociation | f_S: 0.85 / 0.17 / 0.00 / 0.00 | Both factors needed |
| ANOVA interaction | F(1,96)=4497, p<10⁻⁶⁰, η²=0.23 | Massive interaction |
| Dose-response | Threshold at C≈0.2, saturation C≈0.5 | Analytically predicted |
| Triangle sweep | ON: monotone ↑ with ⟨tri⟩. OFF: flat at 0 | Triangles=signal |
| Threshold robustness | F>2800 for all τ∈[0.2, 0.8] | Metric-independent |
| Form invariance | w², w⁴, √tri all significant (F>1000) | Not form-specific |
| Scale stability | f_S≈0.85 for N=50→180, F increases | Effect strengthens |
| Bifurcation theorem | Pitchfork for any convex g with g''(0)>0 | Structurally inevitable |
| Dynamic stability | corr>0.94 at ε=0.5, unique basin | Robust equilibrium |

## Gap Status

| Gap | Status | Key evidence |
|-----|--------|-------------|
| Bifurcation universality | ✅ Closed | Theorem + 6 families |
| Dynamic stability | ✅ Closed | 3 initial conditions converge |
| Lorentzian nature | 🟡 Partial | Symmetry correct, cones open |
| Continuum bridge | 🟡 Partial | Landau analogy, CG 40-60% |

## Reproduction

```bash
pip install numpy networkx matplotlib scipy
python scripts/01_double_dissociation.py   # ~5s
python scripts/02_dose_response_...py      # ~18s
python scripts/03_mechanistic_ablation.py  # ~5s
python scripts/04_robustness_blindage.py   # ~23s
python scripts/05_bifurcation_theorem.py   # ~2s
python scripts/06_dynamic_stability.py     # ~2s
```

All scripts are self-contained. Total runtime: ~1 minute.
