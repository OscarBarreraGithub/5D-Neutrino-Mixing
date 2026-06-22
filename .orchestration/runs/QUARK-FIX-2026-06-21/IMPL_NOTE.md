# QUARK-FIX-2026-06-21 — implementation note

Implementer lane (Opus). Implements PLAN v3 (dual-APPROVED) to the letter. NO git
commits; all changes left in the working tree. The big 1M scan (Phase II) and the
website/docs restatement (Phase III) are NOT done (out of scope for this pass).

HEAD at implementation: `5328a22` (working tree only).

---

## Code changes (per item)

### B2 — ε_K SVD column rephasing to PDG CKM convention
`quarkConstraints/fit.py`
- NEW `_branch_free_unit_phase(z)` (~line 234): returns `z/|z|` with a documented,
  signed-zero-safe fallback (`+1` when `|z| <= _CKM_PHASE_TOL = 1e-300`). Never routes
  a near-zero pivot through `np.angle` (avoids the `np.angle(-0.0+0j)=π` branch).
- NEW `_real_sign(x)`: `+1.0 if x.real >= 0.0 else -1.0` (so `±0.0 -> +1.0`).
- NEW `_rephase_to_pdg_convention(U_L_u, U_R_u, U_L_d, U_R_d)`: pins FIVE entries real
  ≥0 — `V_ud, V_us, V_cb, V_tb` (the four PDG anchors) PLUS `V_ts` (the 5th condition
  that removes ALL residual rephasing freedom) — leaving the single physical phase in
  `V_ub`. Global anchor `pa_0=1` (arg(V_ud)=0). All phase arithmetic uses branch-free
  unit complex numbers. Applies the SAME phase to L and R columns of each Dirac field
  (masses preserved). Idempotent by construction.
  - NOTE/refinement vs PLAN §1.2: pinning only the 4 named entries real does NOT fully
    fix the convention (a 1-parameter residual remains that shifts arg(V_ub) — verified
    numerically). A 5th real-positivity condition is mathematically required for full
    convention-stability. I pin `V_ts` real as the 5th; this keeps `V_ub` complex
    (carrying the physical phase) exactly as PLAN §1.3 requires, and gives the
    rephasing-invariance the audit demanded (×6 collapse). This is a faithful
    realization of approach A, not a deviation in intent.
- Wired into `mass_matrix_observables` (~line 250): rephase before `ckm = U_L_u† U_L_d`.

Verified at the default benchmark (r=0.25): V_ud=0.974, V_us=0.225, V_cb=0.0415,
V_tb=0.999 (all real ≥0), |V_ub|=0.00368, J=3.12e-5 — physical PDG CKM. Rephasing-
invariance (arbitrary physical column rephasing -> identical CKM, ~2e-16), determinism
(same point twice -> bit-identical), masses/|CKM|/Jarlskog invariant — all confirmed.

### B3 — ΔF=2 O4/O5 un-swap + 1/(2 m_M), ALL SIX ME sites (together)
Corrected coefficients (GGMS Eq. 8, M12-ready): `⟨O1⟩=(1/3)`, `⟨O4⟩=(R/4+1/24)` (singlet,
LARGE), `⟨O5⟩=(R/12+1/8)` (crossed, SMALL).
- `quarkConstraints/deltaf2.py`: `_kaon_matrix_elements` (lines ~720), `_meson_matrix_elements` (~908).
- `quarkConstraints/modern/phenomenology.py`: `_kaon_m12_np_from_bridge_match` (~452),
  `_evaluate_epsilon_k_from_bridge` (~482), `_evaluate_delta_mk_from_bridge` (~516),
  `_meson_matrix_elements_inline` (~578).
Verified: ⟨O4⟩/⟨O5⟩ (per-bag) = 2.8532 at R=25.746 (singlet > crossed); ⟨O1⟩ = (1/3)f²mB1
exact; absolute ⟨O1_K⟩ = 2.2128e-3 GeV³ (PLAN §6.1: ≈2.21e-3).

### B1 — CGHNP Zbb fermion-KK retranslation
`quarkConstraints/rs_ew_couplings.py`
- `_casagrande_zbb_B_profile` (~1835): corrected diagonal bracket
  `1/(1+2c)·(1/(2F²) − 1 + 2F²/(3−2c))` (c-sign fixed in both denominators, F²=2f² factor
  restored).
- `build_rs_zbb_fermion_kk_mixing` (~957): restructured flavour sums to CGHNP form —
  `B_d = B_correct(c_d3, f_d3) + (1/(2 f_d3²))·Σ_{i=0,1} row_ratio_i/(1−2 c_d_i)` and the
  symmetric `B_Q`. The light-gen sum now carries the COMMON 3rd-gen `1/(2 f_3²)` factor,
  NOT a per-light-gen full bracket (fixes the F²(c_d3)/F²(c_d_i) ~ 10²-10⁴ inflation).
Verified: B_d/B_Q now POSITIVE (was negative -> sign flip), δg_L^b > 0 (SM-Zbb sign),
δg_R^b < 0; corrected fermion piece is small (< gauge piece).

### M2 — EW001 ΔT physical-M_KK vs Λ_IR convention
`quarkConstraints/oblique_stu.py`
- NEW `GAUGE_KK_ROOT_NN = 2.448687135269161`, `PHYSICAL_MKK_OVER_LAMBDA_IR_SQ = x1²`.
- `minimal_rs_t_coefficient` / `custodial_rs_plr_t_coefficient`: multiply the ΔT
  coefficient by `x1²` so the geometric-Λ_IR coefficient is expressed in the physical-M_KK
  convention (matching the ΔS=c_S calibration). Equivalent to ΔT using `(v/Λ_IR)²`.
  Module docstring + `OBLIQUE_STU_RS_PROXY_V1` updated.
Effect: minimal EW001 floor moves from ~6 TeV to ~16 TeV (chi2 at 6 TeV: 3.6313 -> 519.55).

### M5 — tag_result "proxies"≠"proxy" -> structured tag-hint
`scripts/run_full_catalog_scan.py`
- NEW `_structured_tag_class(diag)` (reads `diag["tag_class"] ∈ {rigorous,proxy,partial,stub}`)
  and `_mentions_proxy(text)` (plural-robust "proxy"/"proxies").
- `tag_result`: honors the structured `tag_class` (after the stub/missing guards) before
  prose; replaces both `"proxy" in needs_text` and `"proxy" in status_text` with
  `_mentions_proxy(...)`; added "proxies resolved" to the resolved-phrase allowlist.
`flavor_catalog_constraints/primary/top_higgs_ew/T001.py` + `T002.py`: add
`"tag_class": "proxy"` to the evaluated diagnostics.
`flavor_catalog/website/scripts/build_scan_explorer.py` (~303, F7): veto count now
requires `tag ∈ {rigorous, proxy}`, aligning the explorer with the harness/comparison
(a partial/stub HARD failure no longer counts as a veto).
Verified: T001 now tags `proxy` (was `partial`).

### M1 — T010 R_b veto -> T011-style loose-edge NP-shift budget
`flavor_catalog_constraints/primary/top_higgs_ew/T010.py`
- `ObservableBudget`: ADDED `central_residual` and `hard_veto_budget` fields.
- `_build_budget`: computes `central = |exp − sm_fit|`, `hard_veto_budget = central + combined`
  (mirrors T011._build_budget).
- `ObservablePull`: ADDED `np_shift`, `hard_veto_budget`, and a `np_shift_ratio` property
  (`|np_shift| / hard_veto_budget`). `_observable_pull` sets `np_shift = predicted − sm_prediction`.
- `evaluate`: veto scalar is now `max(np_shift_ratio)` (was `max(abs_pull)`); `budget`
  reported is `hard_veto_budget`. `_pull_diagnostics` surfaces np_shift/np_shift_ratio/
  hard_veto_budget/central_residual.

**Recomputed M1 headroom (from code, NOT hardcoded):**
- R_b SM-limit pull vs experiment `(sm_prediction − exp)/combined_sigma = −0.99601...`,
  so the OLD 1σ gate passed the SM by only **0.003985844028832619σ** (≈0.004σ — confirms
  the audit; exact code-derived value).
- combined_sigma(R_b) = 6.726812023536855e-4; central_residual = 6.700000000000039e-4;
  hard_veto_budget = 1.3426812023536894e-3.
- NEW gate: SM-limit NP shift = 0 -> ratio = 0 -> SM passes by construction (no knife-edge).
- Floor behavior confirmed at fixture points: minimal T010 vetoes at 3 TeV (ratio 1.67),
  PASSES at 8 TeV (0.236) and 15 TeV (0.067) — the spurious ~108 TeV artifact has
  collapsed; T010 ceases to be the minimal-floor driver above a few TeV.

---

## Tests

### NEW literature-anchored ABSOLUTE pins (M7 — computed from convention dictionary / paper)
- `tests/test_epsilon_k_physics.py` class `TestB3GGMSMatrixElementCoefficients`:
  - O4/O5 coeff ratio from GGMS rationals (≈2.853 at R≈25.7); production must match.
  - singlet O4 > crossed O5 at R≫1.
  - ⟨O1⟩ = (1/3)f²mB1 (catches the missing 1/(2m_M)).
  - absolute ⟨O1_K⟩ ≈ 2.213e-3 GeV³.
- `tests/test_quark_fit.py` (B2): rephasing-invariance (×6 collapse -> CKM invariant),
  determinism (bit-identical), PDG-convention (4 entries real≥0, V_ub complex),
  invariants (masses/|CKM|/Jarlskog), idempotence.
- `tests/test_rs_ew_phase6a_zbb_fermion_mixing.py` (B1): `_manual_B` REPLACED by
  `_cghnp_diagonal_bracket` computed from the `c=−c_repo, F²=2f²` dictionary (NOT a byte-copy
  of production); new `test_cghnp_diagonal_bracket_matches_convention_dictionary_and_production`
  (production must agree with the dictionary oracle) + sign pin (UV b_R -> δg_L^b>0) +
  magnitude pin (`test_casagrande_zbb_fermion_piece_is_smaller_than_gauge_piece`).
- M1 literature-anchored pin: at the SM-limit (NP shift 0) T010 ratio = 0, SM passes by
  construction — added in `test_T010.py` and `test_rs_ew_phase6a_zbb_fermion_mixing.py`.
- M5: `test_plural_proxies_in_needs_human_is_tagged_proxy_not_partial` and
  `test_structured_tag_class_hint_overrides_prose_matching` in test_full_catalog_scan_harness.py.

### Updated snapshot pins (legitimately moved — re-derived at the benchmark)
- `tests/test_quark_deltaf2.py:164` ε_K 1.9286313761001348 -> 20.667539197050676 (B2+B3;
  removed the now-stale `1.92<=...<2.37` band).
- `tests/test_quark_plot_data.py:33` deltaf2_max_ratio -> 20.667539197050676 (B2+B3).
- `tests/constraints/primary/kaon/test_K001.py`: run/unrun 0.26555/2.84474 ->
  1.4067965606347435 / 4.932284791072436 (B3-only, fixed couplings); a parametrized
  expected_pass flipped True->False (corrected kaon ε_K is larger -> small-coupling point
  now crosses the veto at ratio ~1.09); related loose bounds relaxed.
- `tests/constraints/primary/beauty/test_B003.py`: predicted 4.7604e-14 -> 4.696762557259569e-14,
  ratio 0.18065 -> 0.1782339184264785 (B3; budgets NOT touched — they are ⟨O⟩-invariant anchors).
- `tests/constraints/primary/beauty/test_B001.py`, `test_B002.py`, `tests/constraints/primary/kaon/test_K002.py`:
  predicted/ratio/phase snapshots re-pinned (B3, + B2 phase for B002); loose pass/fail bounds
  relaxed where the corrected ⟨O⟩ legitimately shifted the ratio. (Done by a subagent under
  strict rules: NEVER touch anchor budgets; only predicted/ratio/expected_pass.)
- `tests/test_rs_ew_phase6a_zbb_fermion_mixing.py`: δg_L^b/δg_R^b 2.3519e-5/-5.4919e-4 ->
  3.7365283558008076e-06 / -1.9279722543603958e-05 (B1); T010/T011 predicted+ratio re-pinned
  (B1+M1).
- `tests/constraints/primary/top_higgs_ew/test_T010.py`: `_manual_t010_ratio` helper rewritten
  to the M1 NP-shift budget formula; SM-limit M1 anchor pin added.
- `tests/constraints/primary/top_higgs_ew/test_EW001.py`: chi2 3.6312587955526667 -> 519.5501818374169
  (M2); pass/fail parametrize updated (6 TeV now FAILS; 18 TeV passes — ~16 TeV floor).
- `tests/test_rs_ew_custodial_pr1.py`: `test_15tev_minimal_zbb_veto_survives...` renamed/rewritten
  to assert the CORRECTED behavior (minimal T010 no longer vetoes at 15 TeV; still vetoes at 3 TeV).
- `tests/test_rs_ew_custodial_pr2.py`: T010 oracle ratio 1.4680590546247065 -> 0.2364937629532263 (M1).

### Anchor budgets explicitly NOT moved (⟨O⟩-invariant)
- K001:127 anchor.budget 3.0370868171657726e-4; B003 budgets 5.844e-12 / 2.635167648629676e-13;
  C002 anchor.budget/value/m12_budget_gev (C002 predicted/ratio are in-test recomputations, so
  they did not need re-pinning).

---

## Deviations / notes for the reviewer
1. **B2 needs a 5th real-positivity anchor (V_ts).** PLAN §1.2/§1.3 named only the 4 PDG
   entries; pinning only those leaves a residual rephasing that shifts arg(V_ub) (the
   convention is NOT fully fixed). Mathematically a 5th condition is required. I pin V_ts
   real (keeping V_ub complex as §1.3 requires) — this realizes approach A faithfully and
   yields the rephasing-invariance the audit demanded. Flagging because it is a (minor)
   spec refinement, not a verbatim transcription.
2. **ε_K benchmark moved to 20.67, not the plan's ≈3.32.** The plan's ×1.7227 (-> 3.3225)
   was a B3-ONLY monkeypatch figure (reviewer §1.3). With B2's phase-convention fix ALSO
   applied (as PLAN §1.4 mandates: re-pin after BOTH land), the convention-stable ε_K at the
   default benchmark is 20.67. This is expected and correct (the CKM is now physical PDG form;
   b/d systems still pass; determinism + rephasing-invariance verified).
3. **K001 expected_pass flip + EW001 6 TeV pass->fail + pr1 15 TeV veto collapse** are
   physics-correct consequences of the fixes (larger kaon ε_K; ×6 ΔT; collapsed T010 artifact),
   not regressions. Each was re-derived from code and the test intent updated with a comment.
4. No production-code change outside the six plan items. No git commit. Phase II (scan) and
   Phase III (restatement) NOT performed.
