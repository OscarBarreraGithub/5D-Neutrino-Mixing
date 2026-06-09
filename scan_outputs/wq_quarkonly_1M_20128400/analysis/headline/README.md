# Headline figures — WQ quark-only 1M r×M_KK scan (run 20128400)

Five large, single-message, publication-quality figures distilled from the
1,000,000-row scan cache (`../../wq_quarkonly_cache.parquet`). Quark-only,
non-custodial RS. All figures read-only from the cache; regenerate with
`_make_headline_figs.py` (matplotlib Agg).

Definitions used (verified against the cache):
- **rigorous strict survival** = all 11 rigorous (`tag=rigorous`) HARD
  constraints pass = `survives_strict` = `n_excl_rigorous == 0`.
- **veto fraction** of a constraint at a given M_KK = fraction of non-skipped
  points where that constraint fails (`*__pass == False`).
- **reach** of a constraint = the highest M_KK grid point (TeV) at which it
  still vetoes >50% of points.
- `partial`-tag constraints (always-fail placeholders, e.g. B034/EW003) are
  **excluded** from survival and the reach plot — they are not part of the gate.
- 121,293 rows (12.1%) are `skipped` (fit failures) and are dropped throughout;
  878,707 non-skipped rows are used. M_KK axes are floored at 4 TeV. Curves are
  never silently pooled across r (Fig 1 shows one line per r); the reach and
  bulk-mass figures aggregate over r/M_KK by design and say so.

| File | One-sentence message |
|------|----------------------|
| `fig1_survival_vs_MKK.png` | Rigorous strict survival is ~0 below ~20 TeV and ~1 by 30 TeV — a sharp, nearly draw-independent floor at 25–30 TeV (only the smallest fit draw r=0.05 already survives at 20 TeV). |
| `fig2_constraint_reach_ranking.png` | **Money plot:** Z→bb (rigorous, reach 20 TeV, true >50% crossing 20–30 TeV) sets the floor and out-reaches every other constraint by ~4×; the next contenders are the proxy collider/EWPT bounds (CR001, EW001) at ~5 TeV, while εK and the rest sit at 1–3 TeV. |
| `fig3_Zbb_veto_vs_MKK.png` | Z→bb vetoes ~100% of points up to 20 TeV and drops to 0 by 30 TeV, dwarfing the KK-gluon collider bound (dead by 7 TeV) and εK (flat at 0 over the whole ≥4 TeV range). |
| `fig4_bulk_mass_localization.png` | The RS geometric flavor hierarchy: gen-1/2 quarks are UV-localized (c>1/2), 3rd-generation c_Q (~0.39) and c_u (top, ~0.33) cross below c=1/2 into IR localization. |
| `fig5_yukawa_hierarchy.png` | Fitted up- and down-type Yukawa singular values rise monotonically SV1<SV2<SV3, spanning the quark-mass hierarchy (log10 scale). |

## Key numbers

- **Z→bb (T010) reach: 20 TeV at grid resolution; true >50% veto crossing
  20–30 TeV.** Veto fraction: 100% (≤15 TeV) → 89.2% (20 TeV) → 0% (30 TeV).
  This is THE binding rigorous constraint: T010 pass-fraction (0.210) equals the
  overall rigorous strict survival fraction exactly.
- **Per-constraint reach ranking** (TeV, highest M_KK with >50% veto;
  R=rigorous, P=proxy):
  - T010 Z→bb (R): **20** (crossing 20–30)
  - EW001 S,T,U oblique (P): 5
  - CR001 KK gluon → tt̄ (P): 5
  - CR012 KK V⁽¹⁾→WW/WZ/ZZ (P): 3
  - CR013 KK graviton →γγ (P): 3
  - T011 Z→bb asymmetry (R): 2
  - B012 B⁰→K*⁰γ radiative (P): 2
  - K001 εK (R), B003 Δm_s (R), B004 φ_s (R), and CR002/3/4/8/10 + B011/B013
    (P): 1
- εK (K001) and the collider CR* / EWPT bounds are **subdominant** — exactly the
  prompt's physics: the floor is the non-custodial RS Zb_L problem.
- Bulk-mass medians (c = M_5/k): c_Q gen1/2/3 ≈ 0.61/0.56/0.39;
  c_u ≈ 0.66/0.52/0.33; c_d ≈ 0.65/0.60/0.55.

Subtitle on every figure: *"quark-only, non-custodial RS | Z→bb floor ~25–30
TeV is the non-custodial RS Zb_L bound"*.
