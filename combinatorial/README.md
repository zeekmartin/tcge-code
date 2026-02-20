# Combinatorial Reformulation — Sorkin Response

Response to the objection that weighted DAGs are less simple than posets.

## Question

Can TCGE's key results be derived from a purely combinatorial substrate — no real numbers in the foundations?

## Answer

Yes. Lorentzian-like separation (GAP-Signature) emerges from **anisotropic constraint density** alone — a binary edge marking (constrained/unconstrained) with anisotropic distribution. No real numbers required.

However, the temporal arrow (GAP-Arrow) requires additional structure. Orientation is irrelevant for the separation in current models.

## Progression

| Version | Primitive | Lorentzian | Status |
|---------|-----------|:---:|--------|
| v1 (degree) | integer | 0% | ❌ |
| v2 (common neighbors) | integer | 0% | ❌ |
| v3 (gradient) | integer | 0% | ❌ |
| v4 (min-cut) | integer | 0% | ❌ |
| v5b (directed boolean) | boolean | 95% | ✅ |
| v5d (knife-edge) | binary marking | ✅ | ✅ |

## Key Insight

The minimal primitive for signature-like separation is not a directed graph, not real-valued weights, not even a partial order — it is a **binary edge marking with anisotropic distribution**. This is arguably simpler than a causal set's partial order.

Full details in `RESEARCH_LOG.md`.
