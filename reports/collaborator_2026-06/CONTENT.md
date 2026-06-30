# Collaborator report — content source (physics facts; one paragraph per figure)

TITLE: Reproducing RS Warped Flavor & Electroweak Constraints — Fixed-Code Validation Against the Literature
SUBTITLE: Minimal Randall–Sundrum, quark sector. Audit fixes + reproduction of published figures.
DATE: June 2026

## Intro (2 short paragraphs)
We work in a minimal Randall–Sundrum (RS) warped extra dimension; fermion mass hierarchies come from geometric (wavefunction) localization rather than hierarchical 5D Yukawas. An adversarial audit of our quark-sector constraint code found six implementation errors (three blocker-level), now fixed and dual-reviewed: **B1** the Z→bb fermion-KK admixture (mistranslated from Casagrande et al. 0807.4937 — wrong sign and ~190× too large), **B2** the ε_K SVD→PDG rephasing (ε_K was convention-dependent), **B3** the ΔF=2 O4/O5 chiral coefficients + state normalization, **M1** the Z→bb (T010) gate, **M2** the EW001 ΔT convention (oblique T was ~6× too weak), **M5** a tagging bug. Net effect: the headline minimal-model floor changed qualitatively — the old "25–30 TeV (108 TeV at 1σ) Z→bb floor" was an artifact of B1 and collapses to ~5 TeV once corrected.

We validate the fixed code two ways. (i) A **constrained-fit scan**: draw bulk masses + a seed, fit the Yukawas to reproduce the measured quark masses + CKM (a near-unique solution → a thin "locus"). (ii) An **anarchic-forward scan**: draw complex O(1) anarchic Yukawas, compute observables forward, keep mass+CKM-reproducing points (the multi-decade "cloud" the literature uses). The papers scatter anarchic Yukawas; our production scan fits — so the two give very different spreads, which matters for how floors compare to the literature.

**Lane framing (important).** The reproductions in this report are the **anarchic baseline (LANE A)** — the literature strawman — used to validate the code against published figures. They are **NOT the production / our model**: production AS RUN (LANE B) is the fit-aligned MFV with a sharp ε_K wall at **~7 TeV**, and the FPR ideal (LANE C, V5KM-aligned) is **~2 TeV**. So the ~30 TeV anarchic ε_K cloud below is the strawman, not "our floor". Canonical three-lane definitions: `docs/MODEL_CONVENTIONS.md`.

## Headline findings (box)
- The fixed code **reproduces the published RS literature** when run in anarchic-forward mode: ε_K 95%-quantile floor = **10 TeV** (= Bauer 0912.1625 headline), D⁰ funnel containment rising to ~88% at 10 TeV (Gedalia), kaon Im/Re(M12) = **106×** (≈ Blanke's ~100× ε_K problem), C_Bd/C_Bs ≈ 1.
- **Minimal-model floors (corrected):** *typical* anarchic floor with current data ≈ **30 TeV**, set by ε_K (flavor); *existence* (fine-tuned) floor ≈ **18–20 TeV**, set irreducibly by **S,T,U**. Z→bb no longer leads (≈5 TeV).
- Flavor constraints are **tunable** (alignment drives the NP → 0); **S,T,U is irreducible** (no Yukawa freedom). This is precisely why custodial RS exists.

## Floor table (two columns: typical/median vs existence/best-tuned)
Note: the eps_K typical floor below is **LANE A (anarchic)**, the literature strawman. Production AS RUN (LANE B) is ~7 TeV; the FPR ideal (LANE C) is ~2 TeV. See `docs/MODEL_CONVENTIONS.md`.
| Constraint | Typical (median) floor | Existence (best/tuned) floor | Tunable? |
| eps_K | ~30 TeV (LANE A anarchic, current; NOT production) | <~1 TeV | yes (align Im M12 -> 0) |
| Delta m_d (B_d) | few TeV | <~1 TeV | yes |
| Delta m_s (B_s) | few TeV | <~1 TeV | yes |
| D0 | low | <~1 TeV | yes |
| Delta m_K | low | <~3 TeV | yes |
| Z->bb (T010) | ~5 TeV | ~4.6 TeV | no (gauge-dominated; c_Q3 pinned by m_t) |
| S,T,U (EW001) | ~20 TeV | ~18-20 TeV | no (no Yukawa dependence) |
| COMBINED | ~30 TeV (eps_K) | ~18-20 TeV (S,T,U) | — |

## Per-figure paragraphs (each ~1 paragraph). Figure files in figures/.

### fig_explorer_floors.png — Constraint census (our 100k minimal scan)
Per-constraint veto fraction vs M_KK from our 100k minimal-model scan with all fixes in. The corrected ranking: **S,T,U (EW001) leads at ~20 TeV**, ε_K (rigorous) at 7 TeV, **Z→bb (T010) collapsed to 5 TeV** (it was the spurious ~25–30 TeV B1 artifact), Δm_s at 2 TeV; D⁰/Δm_d/Z→bb-A_b never bind. Dropping Z→bb leaves the combined floor unchanged at 20 TeV — concrete proof it is no longer the driver.

### fig_STU.png — S,T,U oblique — Casagrande et al. 0807.4937, Fig. 4
The paper plots the S–T plane with 68/95/99% CL ellipses and the minimal-RS "wedge" swept over M_KK (1→10 TeV) and warp volume L. Minimal RS climbs steeply in **T** — the volume-enhanced T ∝ πL/(2cos²θ_W)·v²/M_KK², the well-known "RS T problem." Ours (right) shows the same wedge against the tighter PDG-2025 ellipse (S=0.026±0.075, T=0.047±0.066 vs the paper's S=0.07, T=0.16), giving a higher floor; our **M2** fix corrected the ΔT convention (it had been ~6× too weak). Because S,T,U carries no Yukawa freedom, it is **irreducible** and sets the existence floor at ~18–20 TeV — which is exactly the motivation for custodial RS (custodial symmetry protects T).

### fig_Zbb_gLgR.png — Z→bb couplings — Casagrande et al. 0807.4937, Fig. 8
The (g_L^b, g_R^b) plane with the Z-pole CL ellipses and the SM point. We use essentially the same R_b/A_b as the paper, so the SM point and ellipses overlay identically (the sanity check). With the **B1** fix, the *total* coupling shift R_b sees — gauge-KK (dominant, +8.7×10⁻⁵ at 20 TeV, reproducing the audit on the nose) plus the corrected fermion admixture — moves g_L^b to less-negative values as M_KK drops, the CGHNP direction. The original bug had the fermion piece wrong-sign and ~190× too large, faking a 25–30 TeV exclusion; corrected, Z→bb is gauge-dominated and bites only to ~5 TeV.

### fig_epsK_cloud.png — ε_K — Bauer et al. 0912.1625, Fig. 4
The canonical RS ε_K problem: anarchic |ε_K| vs M_KK spanning ~6 decades, colored by ε_K consistency, with 5/50/95% quantile curves. The key methodological point: the paper **scatters** anarchic O(1) Yukawas (fat cloud) while our production scan **fits** to masses+CKM (thin locus). We now run the fixed code in an anarchic-forward mode that **matches Bauer's S1 draw setup** (`scripts/anarchic_bauer_s1.py`): |Y| uniform in modulus over [0.1, Y_max=3] with uniform phase, and per-draw **scanned** bulk masses (the RH-top c_u3 flat in Bauer's prior; the rest fixed by the leading-order Froggatt-Nielsen relations from the drawn Yukawa minors). This reproduces the **per-tile ~6-decade cloud** (was ~3.3 with the older fixed-c, |Y|<2.1 generator), the **95%-quantile floor at M_KK ≈ 10 TeV** (= Bauer's "M_g(1) ≳ 10 TeV"), and the upper-bulk amplitude (p75 ≈ 118× SM = Bauer's "typically ~100× SM"). The fitted c's land on Bauer's Table 1 central values. The figure is a HYBRID of two layers from the SAME Bauer-S1 prior: the grey backdrop and the 5/50/95% quantiles come from the FAT 990,000-draw `anarchic_bauer_S1` ensemble (112,024 mass+CKM-gated draws on the 45-tile 1–20 TeV grid), which samples the |ε_K| tails far better than the companion run and so restores Bauer's headline per-tile spread (median ~3.7, max ~4.9 decades, 95%-quantile crossing the upper window edge near M_KK ≈ 10 TeV); the dense blue (Z→bb-passing) and orange (|ε_K|-window) overlays come from the 202,500-draw `anarchic_bauer_s1_zbb` companion (22,880 mass+CKM-gated, 5,401 Z→bb-passing blue, ~5,000 in-window orange). Both layers are the same anarchic prior, so blue is a representative Z→bb-passing subset of the grey distribution. One honest, real difference (not a plotting artifact): our blue (Z→bb-consistent) points fill in only ABOVE ~4.5 TeV, whereas Bauer's blue blankets the plane down to 1 TeV — because our corrected Z→bb (PDG-2025 R_b/A_b + exact coupling) is genuinely STRONGER than Bauer's 2009 version and vetoes low-M_KK draws he would have kept. The spread, floor, axes, coloring and scenario ordering all match. A companion 2×2 S1–S4 panel mirrors the paper's full layout (S3 "little" sits an order of magnitude higher from UV dominance, exactly as published). With current (vs paper-era) constraints the ≥50% floor rises to ~30 TeV.

### fig_consistency.png — Consistency fraction — Bauer et al. 0912.1625, Fig. 5
Percentage of anarchic points consistent vs M_KK: a rising S-curve with 1/M_KK² decoupling, where Z→bb consistency turns on before the combined fraction. Ours reproduces the same ordering and decoupling, shown for both paper-era and current inputs.

### fig_CBd_phiBd.png — B_d mixing — Bauer et al. 0912.1625, Fig. 6
The B_d-mixing new-physics plane: amplitude C_Bd vs phase φ_Bd. The paper's anarchic cloud sits near the SM point (1, 0) but with a genuine wide scatter. With the matched ensemble (Y_max=3, scanned c) our cloud now reproduces that spread: C_Bd median 1.02 with p5 = 0.89 (= Bauer's C_Bd = 0.89±0.17), φ_Bd centered at 0 with tens-of-degrees tails — the full anarchic cloud rather than the artificially narrow locus of the old |Y|<2.1 generator.

### fig_CBs_phiBs.png — B_s mixing — Bauer et al. 0912.1625, Fig. 7
The analogous B_s plane (C_Bs, φ_Bs, and S_ψφ). Ours: C_Bs median 1.04, p5 = 0.94 (= Bauer's C_Bs = 0.93±0.19), φ_Bs centered at 0 with broad tails, and S_ψφ spanning [−0.57, 0.04, 0.62] — reproducing Bauer's "large φ_Bs allowed" (the old narrow cloud could not show this).

### fig_D0_funnel.png — D⁰ mixing — Gedalia, Grossman, Nir, Perez 0906.1879, Fig. 1
The allowed "funnel" in (x₁₂^NP/x, sin 2σ_D). Our matched anarchic D-meson cloud (carrying the complex M12 phase) sinks into the grey allowed funnel as M_KK rises — inside-funnel fraction 50% at 1 TeV → 77% at 3 TeV → 88% at 10 TeV — a faithful reproduction (the wider matched ensemble gives a broader, more paper-like cloud than the earlier fixed-c run).

### fig_ReIm_M12.png — ΔF=2 amplitudes — Blanke et al. 0809.1073, Fig. 2
Re/Im of the KK contribution to M12 for the K and B_s systems. For kaons, Blanke note Im(M12) ≫ Re by ~100× (the ε_K problem); **ours = 106×** (ratio of medians at 3 TeV). For B_s the contribution is O(1). Both reproduce the paper.

### fig_existence_vs_typical.png — Existence vs typical floors (our analysis)
For each constraint we plot the **minimum** ratio-to-bound (the best, fine-tuned point → the *existence* floor) and the **median** ratio (the *typical* anarchic point) vs M_KK. The flavor constraints (ε_K, Δm_d, Δm_s, D⁰, Δm_K) are fully **tunable**: even at 1 TeV the best point beats the bound by ~1000× (alignment sends ε_K ∝ Im M12 → 0), so their existence floor is ≲1 TeV. **S,T,U is irreducible** — min ≈ median (spread only ~1.5×) because it has no Yukawa dependence — crossing the bound only at ~18–20 TeV. Hence even with maximally fine-tuned Yukawas, minimal RS cannot exist below ~18–20 TeV, set by oblique S,T,U; Z→bb is the next irreducible one at ~4.6 TeV.

## Conclusion (1 paragraph)
The fixed, dual-reviewed code reproduces the published RS flavor/EWPT literature across seven figures from four papers, matching their quoted numbers (ε_K 10 TeV floor, D⁰ containment, kaon 100× ε_K problem, C_Bd/C_Bs ≈ 1). The corrected minimal-model picture: the *typical* anarchic floor with current data is ~30 TeV, dominated by the (tunable) ε_K constraint; the *existence* floor for a fine-tuned point is ~18–20 TeV, dominated by the (irreducible) S,T,U oblique constraint. The earlier 25–30 TeV Z→bb "floor" was a sign-bug artifact (now 5 TeV). Pushing below ~20 TeV requires custodial protection of T (custodial RS), which is the next branch of this program.
