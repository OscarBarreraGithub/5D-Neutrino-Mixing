# PLOT SPECS — implementation-ready reproduction of each reference figure

Compiled 2026-06-22. For the next agent: each spec says exactly what to draw,
what per-draw quantity it needs, the experimental bound line, and the **shift vs
the paper** we expect. Read `references/REFERENCES.md` for the table view and
`refs.bib` for citations. PDFs in `references/papers/`.

## Scan data access
- Per-draw JSONL: `scan_outputs/wq_quarkonly_20260622T090807/tile-*.jsonl`
  (quark-sector; B/C/K/T/EW constraints). Full-catalog (incl. L001 mu->e gamma):
  `scan_outputs/fix100k_minimal_20260622T080053/tile-*.jsonl[.tmp.*]`.
- Compact matrix: `<dir>/constraint_matrix.parquet`.
- Each row: `params.M_KK` [GeV]; `constraints[ID] = {passes, ratio, severity, tag, proxy_flags, needs_human_physics}`.
  Convention: **ratio = NP / experimental-bound**, so **ratio < 1 = passes**.
- Survival flags: `survives_all_HARD_strict`, `survives_all_HARD_inclusive`.

Common axis helpers: M_KK in GeV -> divide by 1000 for TeV. First gauge KK mass
M_g(1) = x1 * Lambda_IR with x1 ~ 2.4487 (params.xi_KK stores the per-draw x1).

---

## SPEC 1 — epsilon_K (K001)  [PRIMARY TEMPLATE: Bauer et al. 0912.1625 Fig. 4]
- **Reference:** Bauer-Casagrande-Haisch-Neubert 2010, arXiv:0912.1625, Fig. 4 (p.51).
  (Also CFW 0804.1954 Fig. 1 as the anarchy-scatter style; Blanke 0809.1073 Fig. 9 for the tuning view.)
- **Axes:** x = M_KK [TeV], linear, ~1-10 (extend to 30). y = |eps_K| [dimensionless], log, span 1e-5..1e0.
- **Draw:** scatter of all draws, point = (M_KK/1000, |eps_K^total| for that draw).
  Colour: gray = fails Zbb (T010), blue = passes Zbb, orange = passes Zbb AND eps_K.
  Overlay 5/50/95% quantile curves of |eps_K| vs M_KK. Horizontal band = experimental
  |eps_K| = 2.228e-3 +- (allowed NP window).
- **Per-draw quantity needed:** |eps_K^total| (SM+NP). If only `ratio_K001` is stored,
  plot ratio_K001 vs M_KK instead (ratio = |eps_K^NP| / allowed-room). **Flag:** absolute
  |eps_K| not stored -> either reconstruct from ratio*room+SM, or recompute in the K001 adapter.
- **Bound line:** ratio_K001 = 1 (horizontal) on the ratio version; or the
  |eps_K| = 2.228e-3 band on the absolute version.
- **Experimental input — paper vs ours:** paper |eps_K| = (2.229+-0.010)e-3; ours = 2.228e-3
  (PDG). **Essentially identical** -> our scatter/quantile curves should overlay theirs almost
  exactly. Expected M_KK floor ~ 8-10 TeV (median consistent). This is the key methodology-
  validation plot: it should reproduce 0912.1625 Fig. 4 to good accuracy.

## SPEC 2 — D0-D0bar mixing (C001/C002)  [Gedalia et al. 0906.1879 Fig. 1]
- **Reference:** arXiv:0906.1879, Fig. 1 (p.7) — the canonical allowed-region plot.
- **Axes:** x = sin(2 sigma_D) [-1,1] (the NP CP phase combination); y = |x12^NP / x| [0,1].
- **Draw:** grey allowed region from the 1-sigma D-mixing constraints (the funnel: wide at
  sin2sigma=0, pinched to ~0.2 at +-1). Overlay our scan draws as points (one per draw) using
  their (sin 2 sigma_D, |x12^NP/x|); points inside the funnel pass C001/C002.
- **Per-draw quantity needed:** **x12^NP/x** (magnitude of NP contribution to D mixing,
  normalized) AND **sin(2 sigma_D)** (NP phase). **Flag — NOT directly stored.** C001 stores
  `ratio_C001` (QCD-evolved CP-odd M12^NP fraction). C002 stores `proxy_flags`:
  `q_over_p_deviation_proxy`, `phi_d_sine_deviation_proxy`, `central_observed_deviation_proxy` —
  these are 1-D proxies, not the full 2-D (magnitude, phase). **Must recompute** the complex
  M12^NP per draw from the C001/C002 adapter to populate the 2-D plane.
- **Bound:** the funnel boundary (x12^NP <~ 0.012 at the relevant phase).
- **Shift vs paper:** paper x_D=(1.00+-0.25)e-2, y_D=(0.77+-0.18)e-2. **Current (CKM2025):
  x_D=(0.405+-0.043)%, y_D=(0.636+-0.024)%** -> much smaller, better-measured -> the allowed
  funnel is **substantially tighter/narrower** than the 2009 figure. Big expected shift.

## SPEC 3 — Delta m_d (B001)  [no literature figure; build our own scatter]
- **Reference for inputs:** Blanke 0809.1073 Table 3 (DM_d=0.507(5) ps^-1); no Dm_d-vs-M_KK
  exclusion figure exists in the literature.
- **Axes:** x = M_KK [TeV]; y = ratio_B001 (= |M12^NP|_d / allowed-room) [log], or |Dm_d^NP|/Dm_d.
- **Draw:** scatter (M_KK, ratio_B001); horizontal line ratio = 1 = exclusion.
- **Per-draw quantity:** `ratio_B001` (rigorous, stored). Sufficient.
- **Bound line:** ratio_B001 = 1.
- **Shift vs paper:** our DM_d = 0.5069 ps^-1 ~ identical to Blanke. Anarchic floor only
  ~few TeV (B_d much weaker than eps_K). Expect most low-M_KK draws still pass B001.

## SPEC 4 — Delta m_s (B003)  [no literature figure; build our own scatter]
- **Reference for inputs:** Blanke 0809.1073 Table 3 (DM_s=17.77(12) ps^-1, phi_s).
- **Axes:** x = M_KK [TeV]; y = ratio_B003 [log].
- **Draw:** scatter (M_KK, ratio_B003); line ratio = 1.
- **Per-draw quantity:** `ratio_B003` (rigorous, stored). Sufficient. Optionally also plot
  the phi_s phase shift if needed (check B-system proxy fields).
- **Shift vs paper:** our DM_s = 17.766 ps^-1 ~ identical. Floor ~2-3 TeV. Note many draws
  sit right at ratio~1 (the example row had ratio_B003=1.007) -> B_s is a near-saturated
  constraint; the shift vs paper is small.

## SPEC 5 — mu -> e gamma (L001)  [Agashe-Blechman-Petriello hep-ph/0606021 Fig. 8R/9]
- **Reference:** hep-ph/0606021, Fig. 8 (right, p.22) for BR-vs-loc, Fig. 9 (p.23) for
  M_KK reach. Model formula from Perez-Randall 0805.4652 Eq. 30.
- **Axes (our version):** x = M_KK [TeV], log or linear 1-30; y = BR(mu->e gamma) [log,
  1e-15..1e-9].
- **Draw:** scatter (M_KK, BR) over the full-catalog minimal scan (L001 lives there).
  Two horizontal bound lines: **current MEG II 1.5e-13** (scan default) AND **paper-era
  MEGA 1.2e-11** (what 0606021 / Perez-Randall used). Optionally an analytic curve
  BR = 4e-8 * |spurion|^2 * (3 TeV/M_KK)^4 for a reference spurion.
- **Per-draw quantity:** **BR(mu->e gamma)** = `ratio_L001 * 1.5e-13` (since ratio = BR/limit
  and the live limit is 1.5e-13). ratio_L001 is rigorous (LMFV NDA carrier). Stored.
- **Bound lines:** 1.5e-13 (MEG II) and 1.2e-11 (paper era).
- **Shift vs paper:** the limit dropped ~80x (1.2e-11 -> 1.5e-13) -> the allowed region
  shifts up in M_KK by ~80^(1/4) ~ 3x. With C=0.02 (paper) the floor was M_KK > ~21 TeV;
  with the live MEG II C=1.94e-3 it is even higher. Show both lines to make the tightening
  visible.

## SPEC 6 — Z -> b bbar (T010 R_b, T011 A_b)  [Casagrande et al. 0807.4937 Fig. 8]
- **Reference:** arXiv:0807.4937, Fig. 8 (p.50/51) — the (g_L^b, g_R^b) plane.
- **Axes:** x = g_L^b [~ -0.44 .. -0.40]; y = g_R^b [~ 0.05 .. 0.12].
- **Draw:** Z-pole 68/95/99% CL ellipses (from R_b, A_b, A_FB^b fit); RS scatter of draws
  as points at (g_L^b_SM + delta g_L^b, g_R^b_SM + delta g_R^b). SM ref point
  (g_L^b=-0.42114, g_R^b=0.077420). Mark a couple of fixed-M_KK predictions.
- **Per-draw quantity needed:** **delta g_L^b and delta g_R^b** per draw.
  **Flag — NOT stored as numerics.** Only `ratio_T010` (R_b) and `ratio_T011` (A_b) are
  stored. To draw the plane, **recompute delta g_L^b, delta g_R^b** from the Zbb adapter
  (gauge-KK piece + fermion-admixture m_b^2/M_KK^2 piece). The review notes the admixture
  dominates the gauge piece ~23x.
- **Simpler fallback:** scatter (M_KK, ratio_T010) and (M_KK, ratio_T011) with ratio=1 lines.
- **Experimental input — paper vs ours:** paper used **R_b^0=0.21629+-0.00066, A_b=0.923+-0.020,
  A_FB^{0,b}=0.0992+-0.0016** — these are the SAME LEP/SLC values we use. **Therefore the
  ellipses and the RS stripe should look ~IDENTICAL to 0807.4937 Fig. 8.** This is the user's
  explicit sanity-check plot: a near-perfect overlay confirms our Zbb implementation matches
  the canonical one. Any visible shift would flag a discrepancy.

## SPEC 7 — S, T oblique (EW001)  [Casagrande et al. 0807.4937 Fig. 4]
- **Reference:** arXiv:0807.4937, Fig. 4 (p.38/39) — the S-T plane with RS wedge.
  Cross-check: Carena et al. hep-ph/0701055 Fig. 1 (delta g_bL/g_bL vs Delta T).
- **Axes:** x = S [-0.4 .. 0.6]; y = T [-0.4 .. 0.6] (U=0).
- **Draw:** experimental 68/95/99% ellipses; SM stripe; RS prediction. For OUR data, draw the
  blue RS curve/wedge as (Delta S, Delta T) parameterized by M_KK from 1 to 10 TeV, for both
  minimal and custodial models (two curves). Mark M_KK points along it.
- **Per-draw quantity needed:** **actual (S, T) per draw or per M_KK.**
  **Flag — NOT stored as numerics; recompute** from the EW001 proxy in the `needs_human_physics`
  text:
    Delta S = c_S * v^2 / M_KK^2,  with c_S = 30, U = 0.
    Delta T (minimal_rs) = x1^2 * pi * L / (2 c_W^2) * v^2 / M_KK^2.
    Delta T (custodial_rs_plr) = -x1^2 * pi / (4 c_W^2 L) * v^2 / M_KK^2.
  with x1 = first gauge-KK root (params.xi_KK ~ 2.4487), M_KK = x1 * Lambda_IR,
  v=246.21965 (NB: EW001 uses 246 convention), s_W^2=0.23122 (c_W^2=1-s_W^2), L~35.
  Only `ratio_EW001` is stored numerically; (S,T) must be reconstructed.
- **Experimental ellipse — paper vs ours:** paper used **S=0.07+-0.10, T=0.16+-0.10, rho=0.85**.
  **Our fit (PDG 2025): S=0.026+-0.075, T=0.047+-0.066, rho=0.90.** Our ellipse is **tighter and
  more centered on the SM** -> the RS wedge is **excluded at higher M_KK** than in 0807.4937.
  Minimal RS (large +Delta T) is pushed across the band -> floor ~8-10 TeV; custodial RS
  (small/negative Delta T) survives down to ~3 TeV. The shift (tighter ellipse) is the headline.

## SPEC 8 — Collider direct M_KK reach  [CMS-B2G-25-009 arXiv:2603.23454 Fig. 16L]
- **Reference (PRIMARY):** CMS arXiv:2603.23454 (CMS-B2G-25-009), Fig. 16 (left panel) —
  RS KK-gluon -> tt, sigma*BR vs mass, excluded 0.5-5.5 TeV.
  Theory cross-section curve: Lillie-Randall-Wang hep-ph/0701166 Fig. 2.
- **Axes:** x = M(KK gluon) [TeV, ~0.5-7]; y = sigma*BR(g_KK->tt) [pb, log].
- **Draw:** the CMS observed/expected 95% CL upper-limit curve + the RS theory sigma*BR(M_KK)
  curve (from LRW, scaled to our model couplings: light-fermion 0.2 g_s, t_R 4 g_s, Gamma/M~0.17).
  Their intersection = excluded mass range. Overlay a histogram/rug of our scan's
  `params.M_KK` distribution to show what fraction of draws is collider-excluded.
- **Per-draw quantity:** `params.M_KK` [GeV] directly. **This is a MASS CUT, not a ratio** —
  draws with M_KK < 5500 GeV (CMS) are collider-excluded in the tt channel. The CRxxx
  constraints in the catalog already encode the experimental mass edges (CR001 = tt KK-gluon).
- **Cross-checks:** ATLAS arXiv:2512.17856 Fig. 10(c) (g_KK width 30% < 4.1 TeV) and
  Fig. 10(b) (bulk graviton < 1.3 TeV); ATLAS arXiv:2102.13405 Fig. 5(b) (RS1 graviton
  diphoton k/MPl=0.1 < 4.5 TeV — different coupling convention).
- **Shift vs old reference:** LRW 2007 was a 14-TeV *projection* (~5 TeV discoverable);
  the CMS Run-2 *measurement* now *excludes* up to 5.5 TeV with real data. Show the projection
  curve vs the realized exclusion.

---

## Summary of quantities to add to the scan output for full fidelity
If the next agent wants the rich literature figures (not just ratio-vs-M_KK), the
scan must emit these per-draw numerics (currently only in text/proxy form or absent):
1. **delta_g_L^b, delta_g_R^b** (for the Zbb (g_L^b, g_R^b) plane, SPEC 6).
2. **S, T** values (for the oblique ellipse, SPEC 7) — recomputable from M_KK + model coeffs.
3. **x12^NP/x and sin(2 sigma_D)** for D mixing (SPEC 2) — complex M12^NP magnitude + phase.
4. **|eps_K^total|** absolute (SPEC 1) — currently only ratio is guaranteed stored.
Everything else (ratio_<ID>, M_KK, BR(mu->e gamma)=ratio*limit) is directly available.
