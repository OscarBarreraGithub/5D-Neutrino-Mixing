# Report 05 — RS electroweak sector (rs_ew_spectrum, rs_ew_couplings, oblique_stu, zpole, zpole_lfv, rs_charged_current, rs_higgs_yukawas, collider_resonance, rs_semileptonic_wilsons)

**Structural summary.** Architecture: `rs_ew_spectrum` builds an exact gauge-NN Bessel tower and the overlap kernel `a(c) = Σ (M_KK²/m_n²)χ_n(1)Ω_n(c)`; `rs_ew_couplings` maps it to light-Z shifts `δg = s_Z·g_SM·(m_Z²/M_KK²)·(a(c)−a_ref)` and 4-fermion contacts; `oblique_stu`/`EW001` implement the correlated (S,T) ellipse; `rs_charged_current` does W⁰+KK diagonalization with a single δG_F subtraction. Numerically/algebraically verified: χ² construction (2-dof 5.991, analytic covariance inverse correct), semileptonic Wilson matching (prefactor `−π/(√2 G_F α λ)` = `1/(2K)`; sign and single-`4G_F/√2` correct, consistent with Buras Z′ primer), CC normalization `C_SM = g₂²/(2m_W²) = 4G_F/√2`, W zero-mode mass entry, η_W = −1 consistent with `s_Z = −1`, and the fixed (B1) Zbb bracket algebra (maps exactly onto CGHNP under `c_CGHNP = −c_repo`, `F² = 2f²`). The x₁² geometric→physical ΔT conversion is correct given the physical M_KK the production scan passes.

### [MAJOR] Higgs-Yukawa B-profile is the un-fixed twin of the B1 Zbb bug (wrong c-sign convention + missing F²=2f²)
- **File:** quarkConstraints/rs_higgs_yukawas.py:292-304
- **Category:** convention-inconsistency / sign / factor-of-2
- **Claim:** `_casagrande_higgs_B_profile` implements the CGHNP bracket **in CGHNP's own convention** but feeds it repo-convention `c` and repo `f_IR`, exactly the error the repo's own B1 audit fixed in the Zbb module.
- **Evidence:** rs_higgs_yukawas uses `denom_left = 1-2c; denom_right = 3+2c; value = (1/denom_left)*(1/F² - 1 + F²/denom_right)` while the audited fix in rs_ew_couplings.py:1875-1901 states: "The previous code used `1/(1 - 2c) * (1/F^2 - 1 + F^2/(3 + 2c))` — the c-sign was wrong in BOTH denominators and the F² = 2f² factor was missing", replaced with `1/(1+2c)·(1/(2F²) − 1 + 2F²/(3−2c))` (dictionary c_CGHNP = −c_repo, F²_CGHNP = 2f²_repo). The Higgs builder passes repo c_L, c_E and repo f_IR. For UV-localized leptons (c > 1/2), `1/(1−2c)` is **negative** → wrong sign of δ_L/δ_E (predicts *enhanced* h→ℓℓ instead of CGHNP brane-Higgs *suppression*), and the dominant 1/F² term is 2× too large. Git log confirms the B1 fix commits touched only rs_ew_couplings.py.
- **Confidence:** high (inconsistency; repo's own audit defines the correct form). Numerical impact currently small (shifts ∝ m_ℓ²/Λ²; off-diagonal LFV exactly zero for diagonal v1 fit), but any h→ff / Higgs-LFV number from this module has wrong sign and magnitude.
- **Fix:** replace the bracket with the converted repo-variable form (reuse `_casagrande_zbb_B_profile`).

### [MAJOR] EW001 (S,T) fit anchors are stale pre-2022 values labelled "PDG 2025" — now floor-driving (known, still live)
- **File:** flavor_catalog/processes/top_higgs_ew/EW001.yaml:88-99 (consumed by EW001.py)
- **Category:** stale-data / statistics
- **Claim:** Active ellipse centred on `S = 0.026 ± 0.075, T = 0.047 ± 0.066, ρ = 0.90` instead of PDG 2024 U-fixed `S = −0.05 ± 0.07, T = 0.00 ± 0.06, ρ = 0.93`; S even has the wrong sign.
- **Evidence:** Logged in docs/audits/oblique_stu_pdg_staleness.md as known, but that doc's June-2026 UPDATE says the anchors now matter because oblique S,T sets the ~18–20 TeV existence floor. Since RS predicts ΔT > 0, the stale positive T central (+0.047) is **anti-conservative**: it weakens the floor relative to current PDG.
- **Confidence:** high
- **Fix:** promote the HEPfit/PDG-2024 anchor block already present in the YAML and re-quote the existence floor.

### [MINOR] `custodial_rs_su2r` is silently relabelled `custodial_rs_plr` in coupling metadata
- **File:** quarkConstraints/rs_ew_couplings.py:905-909
- **Category:** logic-bug / doc-code-mismatch
- **Claim:** For `model_label="custodial_rs_su2r"` the emitted metadata reads `ew_model = "custodial_rs_plr"`, contradicting `custodial_protection_included: False` in the same dict and mislabelling every downstream consumer that reads `metadata["ew_model"]` (EW001 does).
- **Evidence:** three-way enum collapsed to two. Numerically benign today (su2r and plr share the same ΔT coefficient; top-partner loops refused for su2r), but any future su2r/plr split keyed on this metadata, and all provenance records, are wrong.
- **Confidence:** high
- **Fix:** emit the validated `ew_model` string unchanged.

### [MINOR] Oblique proxy hardwires L = 35 while the geometry gives L ≈ 35.9–37.5
- **File:** quarkConstraints/oblique_stu.py:55 (`DEFAULT_RS_VOLUME_LOG = 35.0`), used by EW001.py:437-444 with no override
- **Category:** convention-inconsistency / numerics
- **Claim:** ΔT ∝ L is evaluated with fixed L = 35, but the repo geometry gives warp_log ≈ 35.9 (Λ = 3 TeV) to ≈ 37 (Λ = 1 TeV); the spectrum object carries the exact warp_log, never plumbed into EW001.
- **Evidence:** ln(4.07e15) = 35.94. Underestimating L by 3–7% *weakens* ΔT (anti-conservative), opposite direction of the deliberately "conservative" dropped −1/(2L) term. Custodial coefficient ∝ 1/L is *over*-estimated by same amount.
- **Confidence:** high on mismatch, low severity (few-% floor shift).
- **Fix:** pass spectrum.warp_log into evaluate_rs_oblique_proxy.

### [NOTE] c_S = 30 (physical-M_KK) ⇔ geometric 5.0 v²/Λ², ~20% below the CGHNP brane-Higgs S
- **File:** flavor_catalog/processes/top_higgs_ew/EW001.yaml via oblique_stu.py:291-294
- **Category:** stale-data / wrong-formula (suspected)
- **Claim:** ΔS = 30·(v/M_KK_phys)² ⇔ 5.0·(v/Λ)² geometric vs CGHNP brane-Higgs S ≈ 2π(1−1/L) ≈ 6.1·(v/Λ)²; if the "30" was actually quoted in a Λ-like convention the discrepancy is ×6 not 20%.
- **Confidence:** low
- **Fix:** re-derive c_S from CGHNP/PDG with an explicit convention note.

### [NOTE] M_KK-convention hazard: physical vs geometric carried implicitly through `m_kk_gev`
- **File:** quarkConstraints/scales.py:20 (ξ=1.0) vs scripts/run_full_catalog_scan.py:145 (ξ=2.4487) vs oblique_stu.py (coefficients valid only for physical M_KK)
- **Category:** convention-inconsistency
- **Claim:** `rs_minimal_oblique_proxy` and `EW001._resolve_m_kk_gev` silently assume physical x₁·Λ; low-level quark helpers default geometric. Feeding geometric inflates ΔT/ΔS by ~6 and shifts the floor ~2.45×. Production is correct as wired today; nothing enforces it.
- **Confidence:** high (hazard), no current numeric error found.
- **Fix:** tag the convention on the extra (e.g. `m_kk_physical_gev`) or assert.

### Verified-correct (for the record)
- oblique_stu: χ² 2-dof 5.991 correct; ΔT_minimal = x₁²·πL/(2c_W²)·(v/M_KK_phys)² ≡ πL/(2c_W²)·(v/Λ)² matches Agashe/CGHNP Eq. (147); T-driven floor ≈ 18–20 TeV reproduced. Custodial branch genuinely implements sign-flipped 1/L-suppressed T; P_LR Zbb zeroing actually implemented; su2r correctly unprotected in couplings.
- rs_semileptonic_wilsons: 4G_F/√2 exactly once; signs consistent (reproduces Buras Z′-primer ΔC₁₀ form); ν-blocks no double 1/M².
- rs_charged_current: C_SM = 4G_F/√2 ✓; η_W = −1 consistent; single δG_F subtraction ✓.
- rs_ew_couplings Zbb post-B1: bracket maps exactly onto CGHNP under stated dictionary; prefactor m_b²/(2Λ²) matches; flavour sum uses converted denominator.
- zpole/zpole_lfv: A_f, R_q, A_FB = ¾A_eA_f, radiator calibration, BR inversions check out; charge-state factor 2 correct for i≠j.
- SM chiral couplings t3 − Q·s² correct in both modules.
- collider_resonance: mass-limit proxy only; compares physical M_KK consistently.

**Cross-module inconsistencies:** (1) CGHNP B-bracket in two contradictory forms (fixed Zbb vs unfixed Higgs). (2) su2r/plr enum collapse. (3) L fixed 35.0 vs exact warp_log. (4) M_KK convention enforced nowhere.

**Not reviewed:** tests for these modules; docs/CUSTODIAL_*.md claim-by-claim; collider sidecar YAML staleness; Carena top-partner loop formulas (labelled proxy, not re-derived); CGHNP Eq. (153) custodial coefficient −π/(4c_W²L) independent verification; sign of delta_g_R_b vs CGHNP.
