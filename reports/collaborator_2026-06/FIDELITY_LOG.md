# Figure-by-figure fidelity log — reproduction vs published figures

Pass date: 2026-06-23. Author: figure-fidelity audit (FIXED code).
Scope: the eight reproduction PNGs in `reports/collaborator_2026-06/figures/`.
Method per figure: (1) read the cropped paper figure in `references/paper_figures/`
and the source PDF in `references/papers/`; (2) read our reproduction; (3) list every
concrete difference; (4) fix what is fixable; (5) re-render, re-read, compare; (6) record
honestly what could **not** be matched and why.

Generators (all re-run, 0 execution errors):
- `notebooks/_build_anarchic_reproduction_vs_papers.py` -> `anarchic_reproduction_vs_papers.ipynb`
  -> saves `fig_epsK_cloud, fig_consistency, fig_CBd_phiBd, fig_CBs_phiBs, fig_D0_funnel, fig_ReIm_M12`.
- `notebooks/_build_sidebyside_vs_papers.py` -> `constraint_sidebyside_vs_papers.ipynb`
  -> saves `fig_STU, fig_Zbb_gLgR` (newly wired to disk; see note below).
- `notebooks/constraint_plots_vs_literature.ipynb` also re-executed clean (does not write report PNGs).

> Provenance note for `fig_STU.png` / `fig_Zbb_gLgR.png`: before this pass, **no** checked-in
> script saved these two names — they had been produced by an ad-hoc command that was lost.
> They are now saved reproducibly by `_build_sidebyside_vs_papers.py`.
>
> **Regression fix (2026-06-23, follow-up):** an interim revision had switched the STU/Zbb
> export to a `save_report_axis()` helper that cropped to the OURS (right) panel only, so
> `fig_STU.png` and `fig_Zbb_gLgR.png` rendered as **single (ours-only)** panels — inconsistent
> with the other six report figures, which are all **two-panel side-by-sides (real cropped
> PAPER figure on the LEFT, ours on the RIGHT)**. Both are now **back to two-panel side-by-side**:
> `_build_sidebyside_vs_papers.py` saves the WHOLE figure via `save_report(fig, name)`
> (`fig.savefig(..., dpi=150, bbox_inches="tight")`) — the same helper/approach the anarchic
> build uses for the other six. The LEFT axis imshow's the real cropped paper PNG
> (`casagrande_0807.4937_fig4_ST.png` for STU, `..._fig8_gLgR.png` for Zbb, titled
> "PAPER — arXiv:0807.4937 (Casagrande et al. 2008) Fig. 4 / Fig. 8"); the RIGHT axis keeps the
> improved OURS panel (STU overlay with paper-2008 + PDG-2025 ellipses + RS wedge; densified Zbb
> stripe + Z-pole ellipse). Notebook re-executed, 0 errors. **All 8 report figures are now
> consistent two-panel (paper left | ours right) side-by-sides.**

---

## 1. fig_CBd_phiBd.png — Bauer arXiv:0912.1625, Fig. 6 (upper-left, C_Bd vs phi_Bd, S1)

**Paper plots:** blue RS scatter clustered at (C_Bd, phi_Bd) ~ (1, 0); a black SM cross;
and yellow (68% CL) / orange (95% CL) **experimentally favored region** centered at the
global-fit value, NOT at the SM. Axes: C_Bd in [0.4, 1.6], phi_Bd in [-15, +15] deg.

**Differences found (before):**
- We OMITTED the experimental 68/95% CL ellipses entirely (the most prominent feature).
- Axis ranges were wrong: x = [0.6, 1.6] (should be [0.4, 1.6]); y = [-12, 15] (should be [-15, 15]).

**Fixed:**
- Added Gaussian 68%/95% CL ellipses (yellow-in-orange, matching the paper's color order)
  centered at **C_Bd = 0.89 ± 0.17, phi_Bd = (-5.8 ± 2.8) deg** — the exact CKMfitter
  global-fit marginals Bauer quotes in the surrounding text (PDF p.57). No correlation is
  quoted, so the ellipses are axis-aligned (documented assumption).
- Set axes to x = [0.4, 1.6], y = [-15, 15], matching the paper.
- Added a legend distinguishing the experimental CL bands from our scatter.

**Result: NOW MATCHES.** Same layout as Fig. 6: our blue cloud sits as a tight knot at
(~1, ~0), the experimental ellipses sit just below/left at (0.89, -5.8 deg), SM cross at (1,0).

**Could-not-match / caveats:** The published contours are CKMfitter *frequentist* regions
and may be slightly non-elliptical; we draw the Gaussian-marginal approximation (centers exact,
shapes nominal). Our scatter has a wider tail than Bauer's (a few draws reach large C_Bd) —
this is a real property of our matched ensemble (|M12^NP|/|M12^SM| ~ 10% with occasional larger
draws), not a styling choice.

---

## 2. fig_CBs_phiBs.png — Bauer arXiv:0912.1625, Fig. 7 (upper-left, C_Bs vs phi_Bs, S1)

**Paper plots:** blue RS scatter clustered near (C_Bs, phi_Bs) ~ (1, 0) — itself **narrow**;
a black SM cross; and **TWO** stacked yellow/orange CL blobs (a twofold ambiguity) at
phi_Bs ~ -19 deg and ~ -70 deg. Axes: C_Bs in [0.6, 2.0], phi_Bs in [-90, +90] deg.

**Differences found (before):**
- Both experimental CL regions OMITTED.
- Axis x = [0.6, 1.6] (should be [0.6, 2.0]).
- (Reported separately:) "our phi_Bs spread looks too narrow vs the paper's ±60+ deg."

**phi_Bs INVESTIGATION (the key item):**
phi_Bs is defined (Bauer Eq. 31) by C_Bq e^{2i phi_Bq} = M12_full / M12_SM, so
phi_Bs = 1/2 * arg(M12_SM + M12_NP) relative to SM. I recomputed it both ways from the
complex-M12 parquet (`anarchic_bauer_s1_m12.parquet`, 3 TeV tile):
- **arg of NP only:** std ~ 41 deg — the NP phase IS uniformly random (no bug there).
- **arg of (SM + NP) [the correct, paper] definition:** phi_Bs p5/50/95 = **-21.5 / 0 / +20.8 deg**.
The full-amplitude phase is **suppressed by |M12^NP|/|M12^SM| ~ 0.10** in our ensemble:
a random NP phase only tilts the *total* by ~arctan(0.1 * sin phi) ~ ±15-20 deg. So our
±20 deg spread is **physically correct and not a bug**. Crucially, the paper's own blue
S1 scatter is ALSO narrow (clustered near phi_Bs ~ 0, within ~±20 deg) — confirmed both in
the figure crop and Bauer's text. The wide ±60-70 deg features in the published figure are
the **experimental allowed CL ellipses** (the bimodal psi-phi fit), which we were missing —
NOT the model scatter. Our scatter width therefore already matches Bauer's scatter width.

**Fixed:**
- Added BOTH experimental CL ellipse pairs: **C_Bs = 0.93 ± 0.19** with
  **phi_Bs = (-19.0 ± 10.8) deg** AND **phi_Bs = (-69.9 ± 10.1) deg** (the twofold
  psi-phi ambiguity, phi_s <-> 90 - phi_s; Bauer PDF p.58).
- Set axes to x = [0.6, 2.0], y = [-90, 90].

**Result: NOW MATCHES.** Two stacked CL blobs at -19 / -70 deg, our narrow blue cloud near
(1, 0) which does NOT reach either ellipse — exactly the paper's message (RS S1 cannot
populate the experimentally preferred non-SM solutions).

**Could-not-match / caveats:** The narrowness of OUR cloud relative to the *experimental
ellipses* is the physics point, not a defect — documented above. Same frequentist-vs-Gaussian
contour caveat as Fig. 6.

---

## 3. fig_D0_funnel.png — Gedalia arXiv:0906.1879, Fig. 1 (D0-D0bar funnel)

**Paper plots:** x = sin(2 sigma_D) in [-1, 1], y = x12^NP/x12 in [0, 1]; a **lavender/light
blue-violet filled funnel** (allowed region) that diverges as ~1/|x| and is capped at y=1; a
**yellow horizontal GMFV band** along the bottom (y <~ 0.12); and a small **red/pink LMFV box**
near x = +1.

**Differences found (before):**
- Our funnel was a plain grey wedge ("0.7" grey), no GMFV band, no LMFV box.
- y-limit was [0, 1.2] (paper is [0, 1]).

**Fixed:**
- Recolored the funnel to lavender (`#c9cdec`) with a violet boundary line, matching the paper.
- Added the yellow GMFV band (`axhspan` to y=0.12) with a "GMFV" label.
- Added the red/pink LMFV box near sin(2 sigma_D) ~ +1 with an "LMFV" label.
- Set y-limit to [0, 1].
- Kept our anarchic D-cloud overlaid at M_KK = 1/3/10 TeV (red/orange/green), which sinks
  into the funnel as M_KK rises (50% inside at 1 TeV -> 88% at 10 TeV; median y 0.38 -> 0.004).

**Result: NOW MATCHES** the paper's styling (funnel + GMFV + LMFV regions and axes), with our
cloud overlaid (an addition the paper does not have, by design).

**Could-not-match / caveats:** The GMFV-band height (0.12) and LMFV-box width are read off the
crop, not from an equation we recompute (Gedalia's GMFV/LMFV bands come from their specific MFV
spurion analysis); they are positioned to match the figure, labeled as such. The funnel boundary
slope uses sin(phi)_exp ~ 0.18 (their D-mixing CP bound), reproduced from the same inputs.

---

## 4. fig_Zbb_gLgR.png — Casagrande arXiv:0807.4937, Fig. 8 LEFT panel (g_L^b, g_R^b)

**Paper plots (left panel):** the (g_L^b, g_R^b) plane, x = b in [-0.44, -0.40], y = g_R^b in
[0.05, 0.12]; nested filled 68/95/99% CL ellipses (yellow/orange/red) with an SM cross; and a
**dense, near-horizontal blue stripe** pinned at g_R^b ~ 0.0774 that slides in g_L^b.
(The RIGHT panel is a Higgs-mass sweep curve — see below.)

**Differences found (before):**
- Our stripe was too **sparse** — it read as scattered discrete dots, not a continuous line.
- The exported PNG's colorbar was clipped.

**Fixed:**
- Diagnosed the sparsity: it is NOT missing data. 4369/5000 points lie in-window, but g_R^b
  varies only over [0.07742, 0.07747] (a razor-thin line ~5e-5 tall), so individual markers
  read as dots. Restyled to plot ALL in-window points (sorted by g_L^b) as a continuous
  navy line PLUS dense small markers colored by M_KK -> reads as a continuous stripe.
- Fixed the export crop to include the colorbar (+ "M_KK [TeV]" label).

**Result: NOW MATCHES** the paper's left panel: nested 68/95/99% ellipses, SM dot, and a
continuous horizontal stripe at g_R^b ~ 0.0774 sliding toward less-negative g_L^b as M_KK drops.

**Could-not-match (documented):** The paper's **RIGHT sub-panel is a Higgs-mass sweep**
(m_h = 60..1000 GeV green curve with a blue custodial point) — a quantity our pipeline does
**not** compute. We reproduce only the LEFT panel. The stripe's right (high-g_L, low-M_KK) end
genuinely thins because M_KK is sampled log-uniform on [1,10] TeV (few very-low-M_KK draws); the
connecting line makes it read as a stripe but the marker density there reflects real sampling.
Note also our points pin g_R^b at the SM value (delta g_R^b ~ 0 at this gauge fixing), as the
side-by-side notebook's caveat already states.

---

## 5. fig_epsK_cloud.png — Bauer arXiv:0912.1625, Fig. 4 (|eps_K| vs M_KK, S1 panel) — VERIFIED

**Paper plots (S1 panel):** |eps_K| (log y, ~1e-5..1e3) vs M_KK [TeV] (~1..10); blue cloud
falling like 1/M_KK^2, orange band = points consistent with the measured |eps_K|, plus
5/50/95% quantile curves.

**Checked:** axes (log-y 1e-7..1e3, x 1..10.5), the orange "consistent" overlay uses the
paper-era window [1.2, 3.2]e-3, navy 5/50/95% quantile curves, SM line. Per-tile consistent
fractions (5.5% at 1 TeV -> ~28% at 3 TeV) and the 50% quantile crossing the window at ~5 TeV
reproduce Bauer's S1 trend. **No fix needed beyond confirming ranges/scales.** Status: MATCHES
(single S1 panel; the paper's 2x2 S1-S4 grid is reproduced separately in the notebook's Fig. 1b).

---

## 6. fig_consistency.png — Bauer arXiv:0912.1625, Fig. 5 (% consistent vs M_KK) — VERIFIED

**Paper plots:** % of points consistent vs M_KK [TeV], rising from ~0 toward a plateau, one
curve per scenario S1-S4, three sub-panels P(eps_K)/P(Zbb)/P(total).

**Checked:** our single-panel analogue plots % consistent vs M_KK for the paper-era |eps_K|
window (orange) and the current strict gate (crimson), x [1, 10.5], y [0, 100]. The orange curve
(paper-era) rises 5% -> 70% over 1-10 TeV, tracking Bauer's S1 P(eps_K) shape. **No fix needed.**
Status: MATCHES as an analogue. (Honest scope: we show our 2 gates on one panel, not the paper's
4-scenario x 3-observable grid; that full grid is not our deliverable here.)

---

## 7. fig_ReIm_M12.png — Blanke arXiv:0809.1073, Fig. 2 (Re/Im M12 planes, K + Bs) — RESTYLED

**Paper plots:** a density-COLORED SCATTER (colorbar = local point count ~1..100). Left (kaon):
|Im(M12^K)_KK/Im_SM| vs |Re(M12^K)_KK/Re_SM|, log-log, x 1e-3..1e3, y 1e-5..1e5; a diagonal
correlated cloud with Im typically ~100x Re (the eps_K problem). Right (Bs): Re,Im /|M12^s_SM|,
log-log, x,y 1e-5..1e3, generically O(1).

**Differences found (before):**
- We used a **hexbin**; the paper is a fuzzy density-colored **scatter**. The hexbin washed out
  the diagonal-cloud character and the bright core.
- Axis extents were slightly off (kaon x 1e-5..1e3; should be 1e-3..1e3).

**Fixed:**
- Switched both panels to a density-colored scatter (color = 2D log point-count, bright/dense
  points drawn last), with a "point count" colorbar on the Bs panel — matching the paper's look.
- Set extents exactly: kaon x [1e-3, 1e3], y [1e-5, 1e5]; Bs x,y [1e-5, 1e3].

**Result: NOW MATCHES** — the diagonal correlated cloud and bright core read like Blanke Fig. 2.
Our kaon Im/Re ratio of medians is **103.6x**, reproducing Blanke's "~100x" eps_K statement.

**Could-not-match / caveats:** Blanke fixes M_KK = 2.45 TeV; we use our nearest tile (3 TeV) —
a small, documented offset. Absolute color counts depend on our ensemble size (25k draws), so the
colorbar range is not numerically identical to Blanke's, but the structure matches.

---

## 8. fig_STU.png — Casagrande arXiv:0807.4937, Fig. 4 LEFT panel (S-T plane) — VERIFIED + wired

**Paper plots (left/minimal panel):** S-T plane, axes [-0.4, 0.6]^2, filled 68/95/99% CL
ellipses (U=0), SM near origin, and the blue RS region sweeping to large positive T as M_KK
drops (the minimal-RS T-problem).

**Checked / fixed:**
- Confirmed axes [-0.4, 0.6]^2, the 68/95/99% CL ellipses (ours = blue/PDG-2025 live EW001;
  paper = grey/CGHNP-2008), SM cross at origin, and our continuous-M_KK (S,T) trajectory rising
  to large T. This matches the paper's left panel direction.
- Wired the figure to disk reproducibly (was previously an unsaved `plt.show()` cell) and fixed
  the export to include the M_KK colorbar.

**Result: MATCHES** the minimal-RS left panel.

**Could-not-match (documented):** The paper's **RIGHT sub-panel is the custodial model**, where
T is protected/flattened. Our (S,T) proxy is the **minimal** case only, so we reproduce the LEFT
panel; the custodial right panel is out of scope for this proxy. Also, we draw the live PDG-2025
EW fit ellipse (tighter, SM-centered) alongside the paper's 2008 ellipse — an intentional update,
clearly labeled, not a mismatch.

---

## Summary table

| Figure | Paper (arXiv, fig) | Status |
|---|---|---|
| fig_CBd_phiBd | 0912.1625 Fig. 6 | now-matches (added 68/95% exp ellipse, fixed axes) |
| fig_CBs_phiBs | 0912.1625 Fig. 7 | now-matches (added 2 exp ellipses; phi_Bs narrowness shown physical) |
| fig_D0_funnel | 0906.1879 Fig. 1 | now-matches (lavender funnel + GMFV band + LMFV box) |
| fig_Zbb_gLgR | 0807.4937 Fig. 8 (left) | now-matches left panel; right (Higgs sweep) not computed |
| fig_epsK_cloud | 0912.1625 Fig. 4 (S1) | verified-matches (no fix needed) |
| fig_consistency | 0912.1625 Fig. 5 | verified-matches as analogue (1 panel, not full grid) |
| fig_ReIm_M12 | 0809.1073 Fig. 2 | now-matches (hexbin -> density scatter, exact extents) |
| fig_STU | 0807.4937 Fig. 4 (left) | verified-matches left/minimal; right (custodial) out of scope |

**phi_Bs verdict (explicit):** our ±20 deg phi_Bs spread is the physically correct
arg(M12_SM + M12_NP) result, suppressed by |M12^NP|/|M12^SM| ~ 10% in the matched ensemble; it
matches Bauer's own (narrow) S1 scatter. The wide structure in the published figure is the
experimental twofold CL fit (now drawn), not the model prediction. No bug.
