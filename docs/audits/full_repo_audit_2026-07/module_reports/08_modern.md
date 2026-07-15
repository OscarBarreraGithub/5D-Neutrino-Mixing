# Report 08 — modern/ subpackage (quarkConstraints/modern/)

**Structural summary.** ~90% frozen-schema boilerplate; physics in three thin layers: `couplings.py`/`matching.py` (wrap mass-basis KK-gluon couplings into tree Wilsons C1_VLL/VRR = g²/6M², C4 = −g_Lg_R/M², C5 = +g_Lg_R/3M²), `evaluation.py` (runs Wilsons M_KK→2 GeV via qcd_running then deltaf2._hadronic_eval_for_system), `phenomenology.py` (vendored byte-copy of the kaon/B/D hadronic evaluation gating the 5-system acceptance used by scan.py). Verified: operator basis, GGMS M12-ready normalization, KK-gluon color/Fierz coefficients, and RG direction/transpose numerically — all correct. `modern` does **NOT** share the paper_0710_1869 rg.py bug: qcd_running.py applies the coefficient-basis (transposed) ADM in the correct direction (C1×0.73, C4×3.54 from 3 TeV→2 GeV — sane). Real problems: budgets/bounds bookkeeping, B-meson scale hygiene, units-inconsistent surrogate projection.

### [MAJOR] ε_K NP budget is the bare central-value difference |ε_exp − ε_SM| with no uncertainties
- **File:** quarkConstraints/deltaf2.py:654-655,795 (vendored at modern/phenomenology.py:436-437,473)
- **Category:** statistics
- **Claim:** Gating budget `abs(2.228e-3 − 2.161e-3) = 6.7e-5` ≈ 0.44σ of BGS theory uncertainty (σ_SM ≈ 1.5e-4), exp/SM errors ignored. A 95% CL budget ≈ |Δ| + 2σ ≈ 3.7e-4, ~5.6× looser; since floor ∝ budget^(−1/2), default overstates modern ε_K floor ~×2.3. Override hook exists but scan.py acceptance path has no override plumbing.
- **Confidence:** high
- **Fix:** Default budget |Δ| + 2σ_combined or explicit documented CL; plumb override into phenomenology sidecar.

### [MAJOR] ε_K reported `bound` (4.18e-4) is not the budget actually used (6.7e-5) — factor 6.24 mismatch in artifact
- **File:** modern/inputs.py:951 vs modern/phenomenology.py:472-474,1405-1407
- **Category:** logic-bug / convention-inconsistency
- **Claim:** Artifact stores bound=4.18e-4 while ratio_to_bound computed against 6.7e-5; consumers reconstructing amplitude = ratio×bound off ×6.24. Third value 2.0e-8 in deltaf2.py:216 for same key. Three inequivalent "ε_K bounds" coexist; only 6.7e-5 gates.
- **Confidence:** high
- **Fix:** Store the actually-used budget in the artifact.

### [MAJOR] B_d/B_s Wilsons over-run below m_b while bags/chiral factor are m_b-scale — dominant C4 overestimated ~25%
- **File:** modern/evaluation.py:177,196-198; modern/phenomenology.py:366-392,542-560
- **Category:** wrong-formula (scale mismatch, effective double-running)
- **Claim:** All systems evolved to mu_had=2.0 GeV including B_d/B_s, but B matrix elements use FLAG bags at m_b (B1=0.87, B4=1.02, B5=0.96) and chiral factor with m_b(m_b)=4.18 — the 4.18→2 GeV running double-counted.
- **Evidence:** evolve(4.18→2.0): C4×1.250, C1×0.946. |M12^NP(B)| via C4 inflated ~25%, tightening B-mixing floors ~12% in M_KK. Chiral factors also mix scales.
- **Confidence:** high
- **Fix:** Per-system mu_had matched to each bag-parameter scale.

### [MAJOR] `matching_to_deltaf2_summary` compares dimensionless surrogate amplitude to dimensionful GeV bounds (and skips RG)
- **File:** modern/matching.py:517-543
- **Category:** units / logic-bug
- **Claim:** Projection computes `effective_amplitude = (3000 GeV)²·|w·C(M_KK)|` (dimensionless) and divides by |M12| budgets in GeV (1.667e-13, 5.844e-12, 1.742e-15, 3.125e-15) — unit-inconsistent, orders of magnitude wrong. Identical latent bug in dead fallback branch of evaluation.py:226-242. No QCD running applied here either.
- **Confidence:** high on inconsistency; medium on production impact (scan path uses hadronic evaluators, not this projection).
- **Fix:** Delete the projection or wire it to legacy dimensionless bounds with evolution.

### [MINOR] Verifier "verification" is schema/self-consistency only — can never catch a physics error
- **File:** modern/verifier.py:357-700,850-895
- **Category:** circular-validation
- **Claim:** Verifiers check schema-ID string equality and internal redundancy of stored fields with abs_tol=0.0; nothing recomputed from couplings/Wilsons against an independent reference. scan.py's verifier_failed_point_count certifies bookkeeping, not physics — none of the bugs above would trip it.
- **Confidence:** high
- **Fix:** Add at least one independent recomputation per verifier.

### [MINOR] Stale experimental inputs
- **File:** deltaf2.py:683,697 (mirrored modern/phenomenology.py:559,569; modern/inputs.py:987,999)
- **Category:** stale-data
- **Claim:** Δm_s corresponds to 17.757 ps⁻¹ (old HFLAV) not current 17.765; Δm_D = 6.25e-15 GeV (x≈0.39%) vs HFLAV 2024 x=0.407% (≈6.5e-15). Verified current: Δm_d ✓, ε_K ✓, κ_ε ✓, f_K/f_Bd/f_Bs ✓ FLAG 2024, B_K ✓, CKM target ✓ (% level). Effects ≤4%.
- **Confidence:** high (values), low (impact)
- **Fix:** Bump Δm_s and Δm_D.

### [NOTE] Asymmetric exclusion semantics across observables (no common CL)
- **File:** deltaf2.py:956-968; modern/phenomenology.py:661-691
- **Category:** statistics
- **Claim:** ε_K allows NP ≤ 0.44σ_theory (hyper-tight) while Δm_d/s/K/D allow NP |M12| up to Δm_exp/2 (100% of measurement) via `max(...)` picking the looser bound — Δm_d,s measured to ≤0.5% with ~10% SM predictions (a ~2σ budget would be ~5× tighter). Per-system "ratio=1" surfaces have no common statistical meaning; max-aggregation inherits this.
- **Confidence:** high that it's intentional-but-uncalibrated policy
- **Fix:** One documented CL convention for all five systems.

## Cross-module inconsistencies
1. **g_s\* = 3.0 (modern/inputs.py:802) vs perturbative g_s(M_KK) ≈ 0.98 (legacy deltaf2 default)** — both carry the same `operator_convention_id = "kk_gluon_tree_v1"`; |C| differs ×9.4 between "identically labeled" evaluations — apples-to-oranges hazard for compare_2007_vs_modern.py.
2. Three coexisting ε_K bounds under one backend key.
3. Vendored B/D constant blocks in phenomenology.py not covered by the pin test named in the comment (kaon only) — drift risk.
4. RG code: modern uses qcd_running.py — verified correct direction and transposed ADM (eigenvalues {−16, 2}); does NOT share the paper_0710_1869 inverted-RG bug; nothing in modern/ imports that module.

## Verified-correct
KK-gluon tree matching coefficients incl. 1/(2M²); GGMS normalization byte-identical across matching.py/deltaf2.py/phenomenology.py; Δm bound uses |M12|<Δm/2 (correct 2|M12| convention, conservative); _nf_for_scale threshold segmentation correct.

## Not reviewed
modern/artifacts.py, bridge_artifacts.py (skimmed via callers); scripts/compare_2007_vs_modern.py (g_s* and CL mismatches are the concrete risks to check there); modern/scan.py lines ~850-1550 and post-1760; modern/__init__.py; tests incl. the claimed kaon-constant pin test.
