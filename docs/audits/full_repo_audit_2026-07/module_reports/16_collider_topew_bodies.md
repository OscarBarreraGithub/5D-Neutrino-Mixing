# Report 16 — Gap-closing: collider_rs + top_higgs_ew constraint bodies (CR002-014, T001-T020, T012/T014, tests)

**Structural summary.** All ~25 bodies share one template: ~250 lines of YAML-anchor plumbing plus a short evaluate() delegating every physics formula to an adapter. No body re-implements a width or re-applies a factor already in an adapter. CR bodies are pure benchmark mass-exclusion proxies (ratio = m_limit/m_proxy ≤1 passes); T bodies are BR-vs-limit or 1σ-pull comparisons. Directions and generation indices uniformly correct. The dominant real defect is inherited: the confirmed ×2 dipole-width bug flows unmodified into four HARD constraints, and their tests hardcode the buggy outputs.

### [CRITICAL] ×2 dipole bug propagates into T003-T006 exclusions — independently confirmed from first principles
- **File:** T003.py:50-56 (also T004, T005, T006), consuming quarkConstraints/top_fcnc.py:216-247
- **Category:** factor-of-2
- **Claim:** T003/T004 (t→qγ) and T005/T006 (t→qg) route straight through the buggy dipole widths; bodies add no correction, so every BR and exclusion ratio is ×2 too large (over-exclusion).
- **Evidence:** Independent numerical spinor-trace computation (4×4 gamma matrices, t at rest, transverse polarizations): Γ = **0.25**·α·m_t·(|λ_L|²+|λ_R|²) = 0.3375 GeV at λ=1, exactly half the code's 0.6751. Gluon width same ×2.
- **Confidence:** high
- **Fix:** 0.5 → 0.25 in both dipole widths in top_fcnc.py (single point of truth).

### [MAJOR] Tests bake the ×2 into hardcoded expected values — they will mask the fix
- **File:** tests/.../test_T003.py:208, test_T004.py:217, test_T005.py:214, test_T006.py:215
- **Category:** test-bug
- **Claim:** Tests assert the buggy outputs to 16 digits (6.068421034828811e-6 photon; 1.118111915877719e-4 gluon = photon × (4/3)(α_s/α), confirming both encode the same ×2). Corrected values: 3.0342105174144055e-6 and 5.590559579388595e-5.
- **Confidence:** high
- **Fix:** Halve the four constants when the adapter fix lands.

### [MINOR] Malformed-extra silently passes in the T-series but fails closed in the CR-series (new silent-pass variant)
- **File:** T014.py evaluate (also T015-T020): `except (AttributeError, KeyError, TypeError, ValueError)` → `_unevaluated_result` with passes=True. Contrast CR002/CR009/CR012/CR014: identical exception set → passes=False.
- **Category:** logic-bug
- **Claim:** A *present but malformed* couplings object (wrong shape, non-finite entry — `_coupling_entry` deliberately raises ValueError on non-finite, which is then swallowed) yields passes=True for HARD T-series constraints; a shape/NaN corruption in z_delta_g_* matrices vetoes nothing.
- **Confidence:** high (behavior); medium (policy intent)
- **Fix:** Distinguish invalid-extra (fail closed) from absent-extra (unevaluated pass), matching the CR convention.

### [NOTE] VLQ pair constraints use m_VLQ = M_KK (KK-gluon scale), anti-conservative for custodial partners
- CR002/003/004/008/010 compare exclusions to common M_KK ≈ 2.45·Λ_IR, but custodial partners (T_5/3, c_tR near 0) are characteristically *lighter* — 1.5-1.6 TeV VLQ exclusions bite later than they should (systematically permissive). Documented NEEDS-HUMAN-PHYSICS proxy (`VLQ_PAIR_MASS_PROXY_ASSUMPTION_V1`).

### [NOTE] CR009 contact-interaction mapping: takes min() over all eight PDG chirality/sign entries (weakest 23.9 TeV) per an explicit stated policy; residual soft spots: Λ_RS ≡ kk_ew_mass ignores the CI g²=4π convention (ratio not apples-to-apples, off by √(4π)/g_KK), and destructive "range_endpoint" entries pooled as one-sided lower bounds.

### [NOTE] Γ_t = 1.41 GeV fixed BR denominator across T001-T008 (PDG measured; NNLO SM ~1.32 would raise BRs ~7%) — consistent within catalog, cite provenance.

## Verified correct
- Γ(t→qZ) vector width ≡ Aguilar-Saavedra form; phase space (1−ρ)²(1+2ρ) correct.
- Γ(t→qh) = m_t/(32π)(1−ρ)²(|y_L|²+|y_R|²) standard.
- T001-T008 wiring: generation indices correct; SM (~1e-14) excluded from HARD ratio; min(left,right) documented.
- h→ℓℓ′ (T018-T020): Γ = m_h(|Y_ij|²+|Y_ji|²)/8π, Γ_h = 4.07 MeV, already charge-summed matching charge-summed limits (no missing/extra 2); inversion exact; pairs correct.
- Z→ℓℓ′ (T015-T017): BR = 2(|δg_L|²+|δg_R|²)/(Σ_SM + LFV), charge factor 2 correct; Σ_SM re-verified vs sin²θ_eff = 0.2315; inversion incl. (1−limit) exact.
- T014: Γ_FCNC weight = 2·N_c·radiator·(|δg_L|²+|δg_R|²), single [i,j] entry (no double count), BR denominator correct; test helper independently recomputes (2.52251/3.64956 verified).
- T012: R_c calibration to 0.1721, A_c formula, charm index, 1σ pull style (already catalogued).
- CR pass logic: ratio direction, single GeV→TeV application (no double 2.4487), invalid input fails closed, CR012 intervals diagnostic-only, CR014 width benchmark carried.
- CR011: INFO, non-vetoing, direction correct.

## Not reviewed
CR001/CR005-007/CR013 (prior audit), EW001-003/T010/T011 (prior audit), adapter internals and YAML values (prior audit, trusted), test_contract/conftest, zpole_charm radiator deep internals.
