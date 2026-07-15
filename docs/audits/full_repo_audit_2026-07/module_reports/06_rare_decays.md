# Report 06 — Rare decays / FCNC (quarkConstraints: bsgamma, rare_b_*, rare_kaon_*, rare_charm_*, semileptonic_lfu, leptonic_tree, edm, top_fcnc, higgs_lfv, mu_e_conversion, lfv_three_body)

**Structural summary.** Two-layer design: (a) SM-anchored observable formulas (Buras-convention Hamiltonians, κ-parametrizations, form-factor integrals) and (b) explicitly-flagged NEEDS-HUMAN-PHYSICS RS "proxy" matchings (KK-gluon overlaps → Z-like penguin). SM layer largely solid: reproduced BR(Bs→μμ)=3.648e-9 (time-integrated, y_s=0.0645), BR(Bd→μμ)=1.05e-10, BR(K+→π+νν)_SM=8.47e-11, BR(KL→π0νν)=2.95e-11, BR(KL→μμ)_SD=0.82e-9, BR(B→τν)=8.6e-5; X_t=1.481 checks out; C7/C8 LL running exponents correct; NP added coherently (|SM+NP|²) everywhere. One confirmed hard bug: factor-2 in top-FCNC radiative-dipole widths. Proxy matchings carry a factor-2 and a cross-module sign inconsistency in the lepton axial coupling.

### [MAJOR] t→qγ and t→qg dipole widths are a factor 2 too large
- **File:** quarkConstraints/top_fcnc.py:215-245
- **Category:** factor-of-2 / wrong-formula
- **Claim:** For the module's stated convention `L = e A_μ q̄ iσ^{μν}k_ν/m_t (λ_L P_L + λ_R P_R) t`, the width is Γ = (α/4)·m_t·(|λ_L|²+|λ_R|²), but the code returns `0.5*alpha*m_t*(...)` — 2× too large; same error in the gluon mode (`0.5*C_F*alpha_s` should be `0.25*C_F*alpha_s`).
- **Evidence:** Direct trace computation gives Γ = |c|²Σλ² m³/(16π) with c=e/m_t ⇒ (α/4)m_tΣλ². Cross-anchored against the μ→eγ master formula (reproduces Kuno-Okada BR(μ→eγ)=384π²|A_R|² exactly). Numerical check: ratio = 2.0 exactly. The Z-vector and Higgs-scalar widths in the same file are **correct** (verified against Γ(t→bW) and standard t→qh).
- **Confidence:** high
- **Fix:** 0.5 → 0.25 in `top_photon_dipole_partial_width` and `top_gluon_dipole_partial_width` (t→qγ/qg BRs currently overstated ×2).

### [MINOR] Z'-proxy C9/C10 matching is 2× the exact tree-level result for the stated Lagrangian
- **File:** quarkConstraints/rare_b_dilepton.py:316-351, 385-388 (mirrored in rare_charm_dilepton.py:327-362, 398-401)
- **Category:** factor-of-2 / convention-inconsistency
- **Claim:** Code multiplies the prefactor −π/(√2 G_F α λ_t M²) by `mu_axial = g_R − g_L`, but the coefficient of μ̄γ^μγ₅μ is (g_R−g_L)/2; exact matching gives half the coded value (and analogously C9). Sign of prefactor relative to `L_eff = +J·J/(2M²)` also opposite to tree-level result.
- **Evidence:** matching (4G_F/√2)λ_t(α/4π)C10 = Δ_qb·(ΔA/2)/M² ⇒ C10 = +π Δ_qb ΔA/(2√2 G_F α λ_t M²); code = 2× with opposite sign. Module declares proxy "order-one coefficients," so normalization-hygiene bug, but directly rescales all NP shifts in B005/B015/B016/B019/B021 and C004-C008 by 2 and flips constructive↔destructive interference with C10_SM.
- **Confidence:** medium
- **Fix:** Halve lepton V/A deltas or document the extra 2 and the sign.

### [MINOR] Lepton axial-coupling sign convention flips between kaon and B/charm dilepton modules
- **File:** quarkConstraints/rare_kaon_dilepton.py:748-754 vs rare_b_dilepton.py:316-333, rare_charm_dilepton.py:327-344
- **Category:** sign / convention-inconsistency
- **Claim:** For the same SM-Z lepton, B/charm modules use axial = g_R − g_L = +g_Z/2, while kaon module sets `electron_axial_delta = -neutral_delta` = −g_Z/2 (and `muon_axial_delta = +neutral_delta` in its Y_NP path); identical RS quark couplings interfere with the SM with opposite relative sign in K vs B/charm modes.
- **Confidence:** medium
- **Fix:** One axial-coupling sign convention across all three modules.

### [NOTE] B+→τν uses repo CKM |V_ub|=0.00368 and f_Bd for f_B+
- **File:** quarkConstraints/leptonic_tree.py:77-79
- **Category:** stale-data
- **Claim:** Default |V_ub|=3.68e-3 below PDG inclusive (~4.1e-3), gives BR_SM=0.86e-4 vs measured 1.09e-4. Documented that live B009 constraints override from YAML anchors — acceptable if YAML actually overrides; sidecar not reviewed.
- **Confidence:** medium (silent-default risk)
- **Fix:** Verify B009 YAML overrides both values in production.

### [NOTE] K+→π+νν experimental anchor may be stale vs NA62 2024/25
- **File:** quarkConstraints/rare_kaon_snd.py (constants fine; experimental side in K004.yaml — not reviewed)
- **Category:** stale-data
- **Claim:** SM parametrization correct (reproduces 8.47e-11), but NA62 2024 BR = (13.0+3.3−2.9)e-11 (>5σ observation) should now be two-sided measurement, not upper limit; sidecar unverified.
- **Confidence:** low
- **Fix:** Check K004/K005 sidecars vs NA62 2024 and KOTO 2025.

### [NOTE] Verified-correct items (no action)
- `g_SM² = 4G_F²M_W²/(2π²)` looks like a typo but is an exact identity equal to Buras' g_SM² ≈ 1.78e-7 GeV⁻² — NOT a bug (rare_b_nunu.py:271-274, rare_kaon_snd.py:338-341).
- Bs→μμ: β-power structure, A_ΔΓ, time-integrated factor, C10_SM=−4.103, C10∓C10′ combinations all correct.
- b→sγ: LL running exponents 16/23, 14/23, 8/3 mixing correct; C7_SM^eff=−0.304 coherent interference correct. Inclusive B→Xsℓℓ partonic kernel matches standard LO.
- R_K/R_K* definitions match Buras et al. 1409.4557.
- D0→eμ LFV two-body formula correct.
- D+→π+μμ: γ0 normalization internally consistent with cited Sahoo et al.
- EDM conversion 1 GeV⁻¹ = 1.973269804e-14 cm correct.

### Cross-module inconsistencies
1. Lepton axial sign: kaon vs B/charm dilepton modules disagree.
2. Lepton coupling normalization: full L/R-difference vs (g_V,g_A) with ½ absorbed — same RS point gets different effective NP normalization per sector.
3. τ_Bs vs 1/Γ_s: <0.1%, note only.

### Not reviewed (session cut short)
- higgs_lfv.py, mu_e_conversion.py, lfv_three_body.py (SINDRUM II, μ→3e, τ-LFV limit currency, Au/Al overlap integrals UNAUDITED).
- tests/ for all modules; catalog YAML sidecars (90%-vs-95% CL check unverified).
- Second half of rare_b_baryon_dilepton.py a_9_9/a_10_10 kernels (structurally checked only).
