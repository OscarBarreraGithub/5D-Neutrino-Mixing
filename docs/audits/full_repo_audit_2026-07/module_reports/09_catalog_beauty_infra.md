# Report 09 — Catalog infra + beauty (flavor_catalog_constraints/{base,point_builder,rs_ew_builder,TEMPLATE}, primary/beauty, secondary/beauty)

**Structure.** Each constraint = one self-contained file implementing `Constraint.evaluate(ParameterPoint) -> ConstraintResult`; experimental anchors from YAML sidecars in flavor_catalog/processes/(secondary/)beauty/ loaded at import with loud validation; physics reached only via physics_adapters/* wrapping quarkConstraints/* cores (no reimplementation drift found — B001/B003 validate YAML GeV anchors against deltaf2.py constants bit-for-bit). base.py enforces real-finite result fields (NaN/complex cannot leak into ratio); registry isolates exceptions as *failing* results.
**ID map:** B001 Δm_d | B002 S_ψKs | B003 Δm_s | B004 φ_s(J/ψφ) | B005 BR(Bs→μμ) | B006 BR(B0→μμ) | B009 B+→τν | B011 B→Xsγ | B012 B→K*γ | B015 B→Xsℓℓ | B016 B+→K+ℓℓ | B017/018 R_K | B019 R_K* | B021 Λb→Λμμ | B022 B+→K+νν | B023 B→K*νν | B025 R_D | B026 R_D* | B032 ΔA_CP(Kπ) | B033 S_φKs | B034 φ_s^sss | secondary: B007 Bq→ee, B008 Bq→ττ, B013 Bs→φγ, B014 B→ργ.
**Experimental inputs notably current** (HFLAV 2025 Δm_d=0.5069, Δm_s=17.766, BR(Bs→μμ)=3.34e-9±0.27, R_K=0.949 SM-compatible with 0.846 row explicitly superseded, Belle II B+→K+νν evidence 2.3e-5 recorded, FLAG 2024 f_Bd/f_Bs).

### [MAJOR] B022 HARD veto excludes the SM/decoupling limit — 1σ window centered on the Belle II excess
- **File:** flavor_catalog_constraints/primary/beauty/B022.py (`_build_budget_band`, `hard_budget = max(combined_upper, combined_lower)`; `ratio = |pred − exp|/σ`, `passes = ratio <= 1.0`)
- **Category:** statistics | convention-inconsistency
- **Claim:** B022 vetoes any point whose *total* predicted BR(B+→K+νν) lies more than 1σ from the Belle II evidence value, so the SM limit (and essentially every RS point with small b→sνν NP, i.e. the whole scan) fails a HARD constraint.
- **Evidence:** budget = combined_sigma only — unlike every other beauty HARD constraint (B001/B005/B018/B025 use `|exp − SM| + σ_combined`, containing SM by construction). SM = 5.58e-6 (HPQCD 2023), exp = 2.3e-5 ± ~0.7e-5 → pull ≈ 1.74e-5, ratio ≈ 2.4 → passes=False, severity=HARD. A 1σ two-sided window also falsely vetoes ~32% even for a true model.
- **Confidence:** high
- **Fix:** Demote B022 to SOFT, or use the family-standard `|exp−SM|+nσ` budget on the NP shift (n≥2).

### [MAJOR — KNOWN, KI-1] B002/B004 NP phase measured against a real SM box amplitude
- **File:** B002.py:277, B004.py:254, used at B002.py:405, B004.py:402
- Confirmed present exactly as documented in docs/KNOWN_ISSUES.md KI-1; no drift beyond KI-1's description. Fix per KI-1 (complex SM box).

### [MINOR] Silent-pass on invalid input is inconsistent with registry policy (HARD constraints)
- **File:** e.g. B005.py evaluate() `except Exception:` branch returns passes=True ("NOT EVALUATED — invalid rs_semileptonic_wilsons"); same pattern in every beauty file (TEMPLATE.py:97-109 design).
- **Category:** logic-bug (policy)
- **Claim:** A malformed or missing extra makes every HARD beauty constraint silently pass, while the same exception raised out of evaluate would be converted by the registry into a **failing** result — a mis-wired point builder passes the entire family with only a diagnostics note.
- **Confidence:** high on behavior; medium on severity (documented "graceful degradation")
- **Fix:** Return ratio=None with a scan-level "unevaluated ≠ passed" counter, or fail-closed in production tagging.

### [MINOR] Stale ΔF=2 core Δm anchors, papered over by a non-physical ps⁻¹→GeV conversion
- **File:** quarkConstraints/deltaf2.py:670,683; consumed via B001.py:320 `gev_per_ps_inverse = code_inputs.delta_m_bd_exp_gev / experimental.value`
- **Category:** stale-data | numerics
- **Claim:** Core GeV values correspond to Δm_d=0.5066/Δm_s=17.757 (older PDG) vs sidecars HFLAV 2025 0.5069/17.766; B001/B003 reconcile by *deriving* the unit conversion (6.5773e-13 GeV·ps instead of ħ=6.58212e-13, 0.07% off) — silently absorbing a data-version mismatch. Impact ≲0.1% of budgets.
- **Confidence:** high
- **Fix:** Update core anchors to HFLAV 2025 and pin conversion to ħ with an assert.

### [NOTE] B001 mixes SM-central (CKMfitter core 0.547 ps⁻¹) and SM-σ (HPQCD 0.062) provenance — shifts central budget ~0.008 ps⁻¹.
### [NOTE] One-sided-limit budgets applied direction-blind (B006; likely B007/B008): |np_shift|/(UL−SM) vetoes destructive interference symmetrically — over-conservative. B007 sidecar mixes 90% and 95% CL limits in one list. Confidence medium.
### [NOTE] Legacy core ΔF=2 budget (max(Δm/2, |exp−SM|/2)) allows 100% NP saturation — catalog correctly overrides it; recorded so nobody re-promotes the core default.

## Verified clean
- B009 B+→τν: correct formula, τ_B+=1.638 ps (not τ_B0), f_B=0.190; SM 0.863e-4 vs YAML 0.865e-4 ✓; exp anchor current HFLAV Dec 2025.
- Bq→μμ core: standard formula w/ ΔΓ_s time-integration, scalar operators with β³ and chiral enhancement; F_BS/F_BD no copy-paste; B005 budget family-standard, inputs current (exp 3.34e-9, SM 3.64e-9 Czaja–Misiak 2024).
- b→sγ: C7_SM_eff = −0.304 correct sign; primed dipoles in quadrature; 8/3 C8→C7 mixing present.
- R_K/R_K* on post-2022 SM-compatible values; R_D/R_D* budgets contain SM; correlation loaded (diagnostics only).
- B032/B033/B034 are INFO non-vetoing — do NOT replicate KI-1 into a vetoing path.
- Infra: base.py rejects NaN/Inf/complex in veto-relevant fields; make_point rejects unknown extras; B001/B003 pin YAML to deltaf2.py rel_tol=1e-12.

## Cross-module inconsistencies
1. Budget-convention split: B022 (pure 1σ, excludes SM) vs family standard (|exp−SM|+σ); CL conventions vary (1σ/90%/95%) without normalization to common veto strength.
2. Exception policy split: in-file `try/except → passes=True` vs registry `exception → passes=False`.
3. Data-version split: deltaf2.py hardcoded constants vs HFLAV-2025/HPQCD sidecars, reconciled by derived non-ħ conversion.

## Not reviewed
Full line-by-line: B011/B012/B015/B016/B021/B032 internals, secondary B013/B014; rare_b_nunu.py internals (ν-flavor factor 3, K vs K* C_L±C_R signs); scalar-operator convention in rare_b_dilepton vs matching side; anchors.py/registry.py full read; per-constraint tests beyond existence.
