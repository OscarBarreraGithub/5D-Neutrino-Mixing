# Plan — Is minimal (non-custodial) RS still better than the SM on FLAVOR?

**Status:** PLAN ONLY (no implementation). Converged from two independent plans
(Claude strategist + Codex gpt-5.4), cross-argued to convergence, with the
adjudications below. **Posture:** the note is designed to be able to land
*negative*; the null hypothesis is "RS reduced to a mass-hierarchy mechanism,"
and the burden of proof is on the surviving-flavor-advantage claim.
**Date:** 2026-06-29.

Supporting material: `/tmp/.../lit_analysis_RS_vs_SM_flavor.md` (literature
analysis, 30+ refs), `docs/FLOOR_SUMMARY.md`, `reports/collaborator_2026-06/CONTENT.md`.

> **LANE NOTE (read first).** This note is about **LANE A — anarchic RS** (the
> literature strawman) specifically, NOT lane B (production AS RUN, ~7 TeV) or
> lane C (FPR ideal, ~2 TeV). So its ~10–30 TeV epsilon_K numbers are the
> **anarchic** floors *by design* — they are the thing the note interrogates, not
> our production model. Canonical three-lane definitions:
> [`docs/MODEL_CONVENTIONS.md`](MODEL_CONVENTIONS.md). Always tag a floor with its
> lane.

---

## 0. The thesis

> **Minimal anarchic RS, as a *flavor theory*, no longer buys a viable TeV-scale
> prediction.** Its one genuine SM-absent virtue is geometric (single-input)
> generation of *both* the fermion-mass hierarchy *and* the RS-GIM FCNC-suppression
> pattern; the SM does not explain the hierarchy this way. But RS does **not**
> reduce the flavor parameter count (Agashe-Perez-Soni); the claim is structured
> anarchy/localization, not fewer knobs. Run forward with anarchy, that same input
> forces LANE-A M_KK >~ 10-30 TeV from epsilon_K and >~ 18-20 TeV *irreducibly* from
> oblique T, while every would-be flavor "prediction" (phi_s, C_Bs, D0, S_psiphi)
> has collapsed onto the SM point, and the S-T ellipse has **recentered on the
> SM**. The SM does better on flavor/CP as a precision null test: after masses and
> tree-level CKM are fixed, loop flavor data are already accounted for with no
> added custodial/alignment/tuning structure. Custodialization rescues the EW sector
> but not epsilon_K.

Answerable *now* (and not in 2008) for three reasons the note foregrounds:
(i) the Higgs is found (RS's raison d'etre moved from "trigger EWSB" to "flavor +
composite Higgs," re-introducing the little hierarchy v/f); (ii) the oblique
ellipse recentered on the SM; (iii) the epsilon_K SM tension sharpened (2025
lattice), sharpening the flavor lever.

## 1. Central claim + space of honest conclusions

**Q:** Does anarchic RS make any flavor statement, at any M_KK reachable
*indirectly* by current/near-term experiment, that (a) the SM does not make and
(b) is not already falsified or collapsed to "~ SM"?

**Fair benchmark:** split observables into a **FIT SET** (masses + tree-level CKM;
the hierarchy RS is allowed to claim) and a **TEST SET** (loop flavor/CP, neutral
meson mixing, EW precision). The SM benchmark uses tree-level CKM inputs only, so
possible NP in mixing is not fed back into the reference fit. "RS does better" then
means predictive surplus on the test set after the fit set is fixed, not just extra
flexibility.

| Conclusion | Evidence bar | Prior |
|---|---|---|
| **(A) RS wins narrowly on flavor** | >=1 generic anarchic-predicted correlation, distinct from SM at >1sigma of current resolution, realized at an M_KK that simultaneously survives epsilon_K AND S,T,U | Unlikely |
| **(B) RS reduced to mass-hierarchy mechanism (+ non-falsifiable FCNC correlation)** | RS-GIM correlation is real and SM-absent, BUT only predicts smallness *given* large M_KK, with the required M_KK above indirect reach for every channel | **Most likely (target)** |
| **(C) RS flavor fully moot / strictly worse** | RS-GIM adds nothing beyond a generic anarchic-TeV model | Overreach; **decline** |

**Stance:** commit to **(B)** by *trying to rescue (A) and failing in a documented,
quantitative way*. Explicitly **decline (C)** (the geometric-hierarchy virtue is
genuine; a referee will correctly defend it). The decisive sub-question (the spine):
*is there any M_KK at which a generic anarchic point both (1) survives epsilon_K and
S,T,U and (2) shows a >1sigma deviation in some still-measurable channel?* If the
deviation region always sits below the survival floor, (A) dies and (B) is
established.

## 2. Section-by-section outline (headline figures)

1. **Framing: the inverted narrative (2008 -> 2025).** Three goalpost-moves.
   *Fig 1 (HEADLINE):* recentered S-T plane — 2008 ellipse vs PDG-2025 ellipse(s) +
   minimal-RS T-wedge vs M_KK + custodial wedge. (Built; needs the exact-anchor +
   2008-overlay polish of STEP 1.)
2. **The one surviving virtue: geometric hierarchy + RS-GIM as a single input.**
   Steelman the pro-RS case: FCNC ~ sqrt(m_i m_j)/v. *Fig 2:* the same f-factor
   localization that gives the 5-decade hierarchy sets the FCNC coefficients. Also
   state the caveat: this is structural motivation, not a reduction in flavor
   parameters relative to the SM.
3. **The flavor wall: anarchy over-saturates epsilon_K.** *Fig 3 (HEADLINE):*
   epsilon_K cloud + 5/50/95 quantiles, two floor sets (Bauer-era 95% -> 10 TeV;
   current-input 50% -> ~30 TeV). *Fig 3b:* Re/Im(M12) 106x kaon asymmetry.
4. **The irreducible wall: oblique T and the existence floor.** S,T,U existence
   floor 18-20 TeV, no Yukawa freedom (min ~ median). *Fig 4 (HEADLINE):*
   existence-vs-typical floors — flavor collapses to <~1 TeV under alignment, S,T,U
   does not move.
5. **Have the predictions collapsed? The SM-overlap test.** The operational
   (A)-vs-(B) discriminator. *Fig 5 (HEADLINE, the figure that decides the paper):*
   per **quark** channel (S_psiphi / phi_s, C_Bd/C_Bs, D0 phase), the M_KK below
   which a >1sigma deviation is generic, with the epsilon_K(50%) and S,T,U(existence)
   survival bands overlaid. **mu->egamma is NOT in this figure** (the lepton sector
   was dropped from the production scan, so it is not in the quark M12 parquet; and
   in a brane-Higgs setup the dipole is UV-sensitive, so Beneke-Moch-Rohrwild /
   Agashe-et-al. style CLFV requires a separately declared bulk/narrow-Higgs setup)
   — it is flagged as a separate lepton-proxy follow-up, not a joint
   channel here. **KI-1 guard:** phi_s/S_psiphi must be computed from the stored
   complex M12 with the *complex SM-box* convention arg(M12_NP) - arg((V_tq V_tb*)^2),
   NOT the B002/B004 real-M12_SM adapter formula (which carries the KI-1 convention
   bug; see docs/KNOWN_ISSUES.md).
6. **Corrected EW side note: Z->bb is not the driver.** *Fig:* (g_L^b, g_R^b) — old
   25-30 TeV was the B1 bug; corrected ~5 TeV. Narrows the real issue to T and
   epsilon_K.
7. **The cures and what they cost.** Alignment / 5D-MFV / shining drive flavor floors
   to ~2-3 TeV but trade "anarchy predicts" for "symmetry imposes." Custodial fixes
   T/Zbb, not epsilon_K. Down-aligned cures must also pay the D0/charm cross-check:
   relieving kaons can move pressure into the up sector. *Fig 7:* S1 (anarchic) vs
   S2 (aligned) endpoints; custodial-vs-minimal epsilon_K overlay.
8. **Dual-language cross-check + |V_cb| caveat.** Map our floor to composite-Higgs
   m_rho >~ O(10) TeV; fold the epsilon_K |V_cb| systematic into the floor as a band.
9. **Verdict.** Land on (B), preserving the (A)-style virtue (geometric hierarchy is
   real; the FCNC correlation is real but non-falsifiable at reachable scales).

## 3. Adjudicated calculation list

### (a) Already done — reuse
epsilon_K cloud + Bauer 10 TeV floor + 106x Im/Re; existence-vs-typical 18-20 TeV
STU floor; B/D/D0 clouds; corrected Z->bb ~5 TeV; minimal + custodial dT proxies
(`oblique_stu.py`); **the 990k-draw M12 parquet** (`anarchic_bauer_s1_m12.parquet`:
re/im M12 for K/Bd/Bs/D + `passes_pdg` + `ratio_eps_K/B_d/B_s/D`, co-registered);
`anarchic_bauer_s1_zbb.parquet` (passes_Zbb, couplings); S1/S2/S4 parquets;
`sidebyside_points.parquet` (S_pred/T_pred/passes_EW001/T010/T011).

### (b) Cheap (post-process / scalar, <=1 day each)
- **C1 Recentered-ellipse floor table** — {minimal dT, custodial dT} x {2008, PDG-2025
  U-fixed, PDG-2025 U-free}. Scalar eval of `evaluate_rs_oblique_proxy` on an M_KK
  grid; root-find ratio=1. Entry: `quarkConstraints/oblique_stu.py`.
- **C2 Current-input epsilon_K quantiles** — re-quantile `anarchic_bauer_S1.parquet`
  (epsilon_K enters linearly via `eps_K_np = ratio*budget`). No re-scan.
- **C3 |V_cb| sensitivity band** — re-evaluate the epsilon_K floor for inclusive vs
  exclusive |V_cb|. Scalar.
- **C4 SM-overlap scoreboard (Fig 5)** — post-process the M12 parquet (M12 ->
  phi_s/S_psiphi, C_Bd/C_Bs, D0 phase) for the **quark** channels only, mask by its
  own `passes_pdg`/`ratio_eps_K`, overlay the **scalar** S,T,U existence floor and
  experimental resolutions (phi_s +-16 mrad; ATLAS/CMS/LHCb). Report four compact
  outputs: best-case/existence floor, generic survival-fraction floor, binding-
  constraint profile, and predictive-residual score (largest surviving deviation in
  sigma units). **No new scan** (see D1).
  **Two guards (Codex sign-off blockers):** (i) do NOT include mu->egamma — lepton
  params were dropped, so it is not in this parquet; it is a separate lepton-proxy
  follow-up. (ii) Compute phi_s/S_psiphi from the stored complex M12 with the complex
  SM-box convention, NOT the B002/B004 real-M12_SM adapter (KI-1 bug, docs/KNOWN_ISSUES.md).
- **C5 S1-vs-S2 alignment endpoints** — existing parquets (`anarchic_bauer_S1/S2`).
  S1 median epsilon_K crossing ~19 TeV vs S2 ~1.8 TeV on the stored grids; include
  D0/charm as the conservative cross-check on down-aligned relief (conservative
  because D0/charm mixing is long-distance-dominated, so it bounds but does not
  hard-exclude the up-sector pressure that down-alignment induces).
- **C6 Partial-compositeness paragraph** — map our floor to m_rho >~ O(10) TeV
  (2507.05332, 1205.5803). Analytic.

### (c) Deferred / out of scope
Full custodial EW implementation (2-4 wk); continuous alignment dial (3-5 d,
downgraded — S1/S2 endpoints suffice); exact-tower vs ZMA (1-3 wk); lepton mu->egamma
scan (cut; leptons cannot rescue the kaon wall, and brane-Higgs dipoles are not a
clean headline constraint without a calculable bulk/narrow-Higgs setup); full PC dual
translation.

## 4. The three adjudicated disagreements

**D1 — decisive calculation: NO NEW SCAN required for Fig 5.** Both reconcilers
inspected the schemas. The M12 parquet (990k draws, 45 M_KK points 1-20 TeV) carries
the M12 phases for K/Bd/Bs/D *co-registered* with `passes_pdg` and the epsilon_K/dF=2
ratios. S,T,U is a deterministic scalar in M_KK (no per-draw column needed -> vertical
overlay). The only mask NOT co-registered is Z->bb (separate parquet) — but Z->bb
floors at ~5 TeV, far below the binding epsilon_K (~30) and S,T,U (~18-20) floors, so
it passes trivially in the regime where the verdict is decided and is immaterial to
Fig 5. **Resolution:** build Fig 5 by post-processing existing data with the
epsilon_K + (scalar) S,T,U masks. **Pre-registered escape hatch:** if a fully-joint
epsilon_K+Zbb+STU+M12 per-draw mask is demanded, OR the deviation/survival crossover
is ambiguous at the 20 TeV grid ceiling, run ONE cheap targeted joint run (extend
`anarchic_bauer_s1_zbb.py` to also emit `re/im_m12_*`, or the m12 path to compute Zbb
on the same draw; ~1.7 h on 48 cores, demonstrated). Not a broad new program.

**D2 — alignment: S1-vs-S2 endpoints, continuous dial DEFERRED.** The note's claim
("anarchy fails typically; alignment rescues flavor by adding structure") is fully
made by the two existing endpoints (S1 anarchic, S2 `common_cd` down-aligned). A
continuous dial requires a non-unique interpolation convention -> an alignment-model
paper, not this reassessment. D0 is the required cross-check on this cure, because
down alignment can shift pressure into up/charm; treat D0 as conservative
(long-distance charm) — a bound, not a hard exclusion. Optional single intermediate
point only if a referee challenges the endpoints as too binary.

**D3 — S-T anchors: KEEP the repo's genuine PDG-2025 anchors (Claude correct on
facts).** Verified against the snapshot `flavor_catalog/references/EW001/
pdg_2025_electroweak_stu.txt` (PDG 2025 Table 10.8, Section 10 revised Nov 2025 by
de Blas/Dittmaier/Kogler): U-fixed **S=0.026+-0.075, T=0.047+-0.066, rho=0.90**;
U-float **S=0.021+-0.096, T=0.04+-0.12**. These are genuinely sourced + SHA-snapshotted
and are *newer* than Codex's proposed "PDG-2024" numbers (S=-0.05, T=0.00) — so the
"stale anchors" premise is factually wrong, and substituting older PDG-2024 values
would be a regression. **Adopt Codex's underlying discipline only:** (i) label the
figure exactly "PDG 2025 Table 10.8, U-fixed" with the numbers, (ii) show the U-free
variant as a second ellipse (defuses the W-mass/U-treatment objection), (iii) add an
explicit, separately-cited **2008 reference ellipse** (S~0.07+-0.10, T~0.16+-0.12) for
the recentering story. Any PDG-2024 overlay is an optional external cross-check, never
the headline.

## 5. Custodial recommendation

**Do NOT build full custodial for this note.** The verdict is custodial-invariant:
epsilon_K is a color-octet KK-gluon LR operator, untouched by the EW SU(2)_R/P_LR
structure. Custodial enters as a **one-table / one-overlay** demonstration via the
existing proxy (`ew_model='custodial_rs_plr'`, `custodial_rs_plr_t_coefficient`):
T/Zbb floor drops toward ~2-3 TeV while the epsilon_K floor stays ~10-30 TeV.
**Pre-registered revisit criterion:** escalate to a full custodial EW fit *only if*
the STEP-3 scoreboard unexpectedly rescues (A) (finds a surviving >1sigma channel),
because only then does "can custodial lower M_KK into reach" become live. This
directly answers the user's "is custodial worth the investment" question: **not for
this note**; custodial is motivated by the Phase-2 RS-EW 100M-scan program, not by
this reassessment.

## 6. Risks / referee objections + defenses

- **ZMA vs exact:** we *reproduce* the published ZMA floors (Bauer 10 TeV, Blanke
  100x), so our ZMA is calibrated; the dual-language m_rho cross-check is independent;
  the verdict is robust to O(1) since the floor/deviation gap is a decade.
- **|V_cb| systematics in the 5sigma epsilon_K tension:** C3 propagates both |V_cb|
  choices as a band; the verdict holds even at Bauer-era inputs (10 TeV >> deviation
  region). The 5sigma is corroborating, NOT load-bearing.
- **Anarchy-prior dependence:** use Bauer's exact prior (modulus-uniform, Y_max=3,
  scanned c) and reproduce his numbers; report S4 (Y_max=12) as robustness; the prior
  only moves the floor *up*.
- **Floor-definition ambiguity:** quote typical (median) AND existence floors always;
  state the (A)/(B) verdict against the typical floor; map them to best-case vs
  generic/survival-fraction language in the scoreboard.
- **"You forgot a channel":** Fig 5 scans all four NP-sensitive channels jointly; the
  lepton sector is flagged future-work (its mu->egamma floor is comparable/lower, so
  it cannot rescue (A); and the brane-Higgs dipole is UV-sensitive, so it is not a
  clean quark-scoreboard input).
- **Fit-locus vs anarchic-scatter:** Fig 5 and all sec-3/5 figures use the
  anarchic-forward ensemble (matching the literature method), not the fit locus.
- **SM benchmark contamination:** use tree-level CKM only for the SM reference fit;
  do not let loop-level NP-sensitive mixing observables define the null-test baseline.

## 7. Critical path (fastest confirm-or-kill, ~1.5 new-compute days)

1. **STEP 0** (no compute): freeze the evaluative standard ("better on flavor" =
   viable + predictive FCNC suppression *without* imposed flavor symmetry). Declare
   the fit/test split: fit masses + tree-level CKM, test loop flavor/CP + EW; score
   predictive surplus, not parameter flexibility.
2. **STEP 1** (~0.5 d, scalar): lock the S-T recentering figure with exact labeled
   PDG-2025 anchors (U-fixed + U-free) + explicit 2008 reference ellipse + RS/custodial
   wedges (C1). *Kill condition:* if recentering does not raise the minimal floor, the
   sec-1 premise is wrong.
3. **STEP 2** (~0.5 d): lock the epsilon_K band — paper-era ~10 TeV, current ~30 TeV,
   |V_cb| band (C2 + C3). No re-scan.
4. **STEP 3 (decisive, ~1 d post-process):** build the SM-overlap scoreboard (C4).
   *Decision gate:* if no epsilon_K+STU-surviving point shows a >1sigma deviation in
   any channel within the <=20 TeV grid -> (A) dead, (B) confirmed. *Escape hatch:*
   extend the grid / run the cheap targeted joint run only if ambiguous at 20 TeV.
   Output: best-case floor, generic floor, binding constraint, predictive residual.
5. **STEP 4** (~0.5 d): custodial-only no-rescue table (proxy config flip).
6. **STEP 5**: S1-vs-S2 alignment endpoints.
7. **WRITE** once four facts are stable: recentered S-T disfavors positive T;
   epsilon_K is the typical wall; custodial doesn't touch epsilon_K; alignment rescues
   only by adding structure.

## 8. Residual disagreements (post-convergence)

- **D1 escape hatch (minor):** the 20 TeV M12-grid ceiling vs the ~30 TeV epsilon_K
  floor — pre-registered as a cheap grid-extension/targeted-joint-run, NOT a blocker.
- **D2 (stylistic):** continuous dial deferred unless a referee asks "how much
  alignment is needed."
- **D3 (resolved on facts):** keep PDG-2025; the older PDG-2024 substitution is
  rejected; Codex's discipline (labels, U-free, 2008 overlay) adopted.

## 9. Changelog

- **2026-06-29 external gap-analysis fold-in:** accepted fit/test + predictive-surplus
  framing, parameter-count caveat, CLFV brane-Higgs dipole demotion, D0-as-alignment-
  cost cross-check, tree-level-CKM SM benchmark discipline, and compact statistical
  output names; rejected broad model-ladder/new-scan/rare-decay expansion as outside
  this LANE-A plan-only kill test. Independent Claude review: APPROVE-WITH-NITS;
  added the D0/charm long-distance-dominated caveat (D0 is a conservative bound, not
  a hard exclusion) per the one flagged nit.
