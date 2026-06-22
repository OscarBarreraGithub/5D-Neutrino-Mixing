# Slice 4 — Catalog adapters & experimental anchors (HARD quark-sector ΔF=2 / radiative)

Auditor: Claude (slice-4, QUARK-AUDIT-2026-06-10). READ-ONLY audit; python executed for
numeric verification. Scope: K001, B002, B003, B004, C001, C002, B011, B012, their YAML
sidecars, `flavor_catalog_constraints/physics_adapters/{deltaf2,bsgamma,ckm_extraction}.py`,
and `tests/constraints/primary/{kaon,beauty,charm}`. Upstream ΔF=2 matrix-element bugs
(O4/O5 swap, missing 1/(2m_M)) and the Zbb mistranslation are CONFIRMED FINDINGS OF OTHER
SLICES; here I audit the adapter/anchor layer and quantify the interactions only.

All numeric claims below were re-derived by running the repo code (constraint constructors,
`evolve_deltaf2_wilsons`, `_ll_bsgamma_running_coefficients`, `repo_default_ckm_phases`).
All 77 tests in the eight relevant test files pass on HEAD.

---

## Summary verdict

The adapter/anchor layer itself is in good shape: anchors are current and provenance-pinned,
unit conversions are exact, NP-room constructions are documented and arithmetically correct,
phase logic for B002/B004 uses the complex amplitude correctly, and the bsγ leading-log
running machinery is numerically correct. I found **0 BLOCKER and 2 MAJOR issues that
originate in this slice**, plus one BLOCKER-level **interaction** with the upstream
matrix-element bugs that must be stated clearly: for LR-dominated kaon points the combined
upstream bugs make K001 **anti-conservative** (ε_K^NP understated ≈ 1.7×), and this layer's
choice of the loose-edge budget (4.5× the central residual) compounds that leniency.
Budgets themselves were NOT calibrated against the buggy ⟨O⟩ layer (they are exp/SM-anchor
constructions), so fixing the core does not require re-deriving any budget.

Counts: BLOCKER 0 (1 upstream interaction flagged), MAJOR 2, MINOR 8.

---

## INTERACTION FINDINGS (upstream bugs confirmed by other slices — impact on this slice)

### I-1. ΔF=2 ⟨O⟩ bugs propagate to all six mixing gates; direction is NOT uniformly conservative — BLOCKER-level interaction (confidence: high)

Context given: `quarkConstraints/deltaf2.py` `_kaon_matrix_elements()` /
`_meson_matrix_elements()` have the O4/O5 basis swap and are missing the 1/(2m_M)
normalization (all ⟨O⟩ are 2× the M12-ready values).

I verified numerically how this distorts the quantities the catalog gates actually consume,
using realistic RG-evolved KK-gluon coefficients (matching pattern C1=+1/6, C4=−1, C5=+1/3
per g_L g_R/M_KK², evolved 3 TeV → 2 GeV, where C4 → −3.24, C5 → +0.28):

| System | VLL/VRR M12: code/correct | LR M12: code/correct |
|---|---|---|
| K (ε_K, Δm_K) | 2.00 | **0.58** (understated) |
| B_d (B001/B002) | 2.00 | 1.06 |
| B_s (B003/B004) | 2.00 | 1.06 |
| D0 (C001/C002) | 2.00 | 0.92 |

(code/correct per-operator factors for the kaon: O1 2.0×, O4 0.70×, O5 5.7×; the
post-running LR combination is C4-dominated, hence 0.58 net. "Correct" = standard
SUSY/BMU-basis M12-ready MEs: ⟨O1⟩→(1/3)f²m B1, ⟨O4⟩→(R/2+1/12)f²m B4/2,
⟨O5⟩→(R/6+1/4)f²m B5/2, R=(m_M/(m_q1+m_q2))².)

Consequences for this slice:

- **K001 (ε_K)**: RS KK-gluon ε_K is typically LR(C4)-dominated (chiral+RG enhancement,
  R_K ≈ 25.7). Code understates ε_K^NP by ≈ 1.7× for such points → K001 HARD veto is
  **anti-conservative**; ε_K-driven M_KK floors are understated by ≈ √1.7 ≈ 1.3× at fixed
  budget. Combined with the loose-edge budget policy (see M-2), the deployed K001 gate can
  be up to ≈ 7.8× looser in ε_K (≈ 2.8× in M_KK) than a central-budget, corrected-ME gate.
- **B003/C001 (|M12| gates)** and the VLL-dominated regions of all systems: predictions 2×
  overstated → conservative (vetoes too strict by ≤ 2× in M12, ≤ √2 in M_KK floor).
- **B002/B004/C002 (phase/Im gates)**: φ_NP = arg(1 + M12^NP/M12^SM) and Im M12^NP inherit
  both the magnitude distortion and a *complex-direction* distortion (the swap changes the
  relative weight of the C4 vs C5 terms, which carry the same phase here, so the dominant
  effect is the magnitude factor in the table). C002's |Im M12^NP| for D0 is ≈ 0.92–2.0×
  depending on operator mix.
- **Budgets are clean**: K001's band (exp/BGS2020), B003's room (HFLAV/CKMfitter/FLAG
  f_Bs√B̂), C001/C002's Δm_D/2 and sin φ_M rooms, and B002/B004's σ-combinations are all
  built purely from YAML/doc anchors — none was tuned against the buggy ⟨O⟩ layer. Fixing
  `deltaf2.py` MEs requires **no budget re-derivation**, only re-running floors.

### I-2. Zbb / rs_ew_couplings bug: no interaction with this slice (confidence: high)

Grep over all eight constraints + the deltaf2/bsgamma adapters + their cores: no import or
use of `rs_ew_couplings` / Zbb quantities. The confirmed Zbb mistranslation does not touch
the ΔF=2 or radiative gates audited here.

### I-3. ε_K Im-part SVD phase-convention dependence (known from another slice)

Not re-derived. Noted: it equally affects B002 (arg M12_Bd^NP), B004 (arg M12_Bs^NP), and
C002 (Im M12_D^NP), since all take the phase/Im of the same convention-dependent mass-basis
couplings. Any upstream phase-convention fix must be revalidated against all four
phase-sensitive gates, not just K001.

---

## MAJOR findings (this slice)

### M-1. Bag-parameter scale pairing inconsistent with the 2 GeV Wilson evaluation point
- Files: `quarkConstraints/deltaf2.py:636-701` (constants), consumed via
  `physics_adapters/deltaf2.py` by K001/B001/B002/B003/B004/C001/C002 (all use
  `mu_had = 2.0`).
- Code values: kaon B4=0.903, B5=0.691 are FLAG 2024 MS-bar(**3 GeV**); B_d/B_s B1=0.87,
  B4=1.02, B5=0.96 are MS-bar/NDR at **μ=m_b** (per `docs/audits/bag_param_inventory.md:31-45`);
  D0 B_i are at 3 GeV (ETM 2015). Wilsons are LL-evolved to **2 GeV** before contraction.
- Correct: C_i(μ)·⟨O_i(μ)⟩ must share μ. Mismatch size (LL, repo's own `run_alpha_s`):
  kaon O4 ≈ 19% overstatement (C4 grows ~(αs(2)/αs(3))^{16/(2β0)} between 3 and 2 GeV);
  B-system VLL ≈ 7% (2 GeV vs m_b), LR larger (O4 anomalous dimension −16).
- Impact on floors: ≤ ~10% on M_KK floors; conservative direction for the kaon LR piece.
- Status: kaon part is explicitly documented in-code (TODO at `deltaf2.py:643-650`);
  the B-meson μ=m_b pairing is only flagged "Yellow" in the inventory and is NOT mentioned
  in the B003 constraint docstring, which calls the path "audited".
- Severity: MAJOR (systematic, affects every mixing gate; partially documented).
  Confidence: high.

### M-2. Test suite pins adapter ≡ core, never pins absolute physics — the confirmed core bugs were structurally undetectable
- Files: `tests/constraints/primary/kaon/test_K001.py:208-226` (and the analogous
  `test_pass_fail_and_numbers_match_*` in test_B003/test_C002, plus
  `test_constraint_matches_shared_bsgamma_module_recomputation` in test_B011).
- All numeric assertions compare the constraint output to a re-invocation of the same
  `quarkConstraints` evaluators (`evaluate_epsilon_k_with_running`, etc.). There is no
  literature-anchored absolute pin (e.g., ⟨O1_K⟩(2 GeV) ≈ 0.00221 GeV³ M12-ready, or an
  M12 benchmark against a published RS/NP point). Consequently the O4/O5 swap and the
  missing 1/(2m_M) pass 77/77 tests, and any future ME regression will too.
- Impact: explains how a BLOCKER-class physics error survived the "1054 passed" gate.
- Recommendation: add 2-3 absolute matrix-element pins (kaon O1/O4/O5 at fixed μ with
  tolerances from FLAG) and one end-to-end |M12| pin against an independently computed
  benchmark.
- Severity: MAJOR (process/test-design). Confidence: high.

---

## MINOR findings

### m-1. K001 hard veto = loose band edge, applied sign-blind
- `K001.py:212-231`: hard budget = |ε_exp − (ε_SM − σ_comb)| = 3.04e-4, i.e. 4.5× the
  central residual 6.7e-5; veto compares |ε_K^NP| (abs). The loose edge is one-sided
  (assumes NP fills the exp−SM gap direction); a *wrong-sign* NP of 3.0e-4 would actually
  sit ≈1.6σ_comb below experiment yet still passes. Documented policy
  (`docs/audits/epsilon_k_sm_decision.md:71-100` — the doc's "uses the smaller of" sentence
  at the upper-edge paragraph contradicts its own worked number 3.0e-4; the code implements
  the worked number). Conservative-for-exclusion by design, but combined with I-1 it makes
  the deployed ε_K gate very permissive. Confidence: high.

### m-2. `_SM_CHOICE_SENSITIVITY = 0.15e-3` hardcoded in K001.py, not in the YAML sidecar
- `K001.py:57`. All other budget inputs are YAML-loaded; this one is a constant whose
  provenance lives only in the audit doc (UTfit/CKMfitter span). Provenance hygiene only;
  value matches the doc. Confidence: high.

### m-3. B003 SM Δm_s provenance mix and stale core exp constant
- SM central is the core constant `DELTA_M_BS_SM = 1.17e-11` GeV (= 17.775 ps⁻¹,
  "CKMfitter") while σ_SM is built from FLAG 2024 f_Bs√B̂_Bs (256.1±5.7 MeV → 2×2.23% =
  4.45% on Δm; arithmetic verified: σ_SM(Δm) = 5.21e-13 GeV, hard M12 budget = 2.64e-13
  GeV). Mixing a CKMfitter central with a FLAG-only σ omits the |V_ts V_tb|² parametric
  uncertainty (~3-5% on Δm) → budget modestly tighter than a full treatment (conservative).
  Separately, the core `DELTA_M_BS_EXP = 1.1688e-11` GeV ↔ 17.757 ps⁻¹ disagrees with the
  YAML exp (17.766 → 1.16938e-11) by 0.05%; it only feeds the *legacy diagnostic* budget
  (`core_legacy_m12_budget` = Δm_exp/2 = 5.844e-12), not the catalog veto. Same class:
  `DELTA_M_BD_EXP = 3.334e-13` ↔ 0.5066 ps⁻¹ vs PDG 0.5069 (0.07%). Confidence: high.

### m-4. B012 "SM" prediction ≡ experimental normalization; budget = σ_exp only
- `B012.py:104` (`sm_value = neutral.value` = 4.163e-5) and `:143-157` (budget = 0.092e-5).
  By construction the no-NP point has zero residual and only the experimental σ is NP room;
  no exclusive form-factor theory σ (which is several × larger) is granted. Documented
  circularity ("normalization at the catalog boundary"); strictly conservative. With the
  C7-proxy issue (m-5) this gate's floors should not be quoted as physical. Confidence: high.

### m-5. B011/B012 C7 "proxy" matching has no loop factor or RS-calibrated normalization
- `quarkConstraints/bsgamma.py:317-358`: C7_NP(M_KK) = (g_bs/g_s)·(3 TeV/M_KK)²·1.0.
  A genuine RS dipole is loop-suppressed (×αs/4π or 1/16π² and a chirality flip), so the
  proxy can overstate C7_NP by O(10²-10³) → B011/B012 vetoes are enormously conservative
  and their M_KK reach is not interpretable as physics. This is prominently documented
  (NEEDS-HUMAN-PHYSICS, `BSGAMMA_RS_MATCHING_ASSUMPTION_V1`, repeated in both constraint
  notes/diagnostics), hence MINOR as a *documentation-consistent* limitation — but any
  floor table mixing B011/B012 with the ΔF=2 floors must caveat this. Confidence: high.

### m-6. B011 SM C7 convention: representative C7_eff(μ_b=4.8 GeV) = −0.304 in a
  rate-normalized formula
- `bsgamma.py:67-77`. Because BR is normalized to the YAML NNLO anchor and only the *ratio*
  |C7_SM+C7_NP|²/|C7_SM|² is used, the absolute C7_SM choice enters only the NP
  interference weight; using the LO-ish −0.304 instead of the NNLO effective ≈ −0.29..−0.38
  (μ_b-dependent) shifts the linearized NP sensitivity by O(10-25%). Dwarfed by m-5.
  Documented as "representative". Confidence: medium.

### m-7. B002 SM-σ anchor uses HFLAV *measured* β (22.63° ± 0.45°) around the in-core CKM
  central (22.55°)
- `B002.py` budget: σ_sin2β = |cos 2β|·2σ_β(HFLAV) = 0.0111. Using the measurement's own
  uncertainty as the SM-prediction uncertainty is mildly circular (a true CKM-fit-without-
  sin2β σ would be larger, e.g. ~0.7° on β). Conservative; documented in the docstring and
  `_SM_PHASE_SOURCE_POLICY`. Central consistency verified: in-core sin2β = 0.7083 vs exp
  0.710 (residual 0.0017 ≪ budget 0.0177), φ_s = −0.0379 vs exp −0.041 (residual 0.0031 ≪
  budget 0.0160) — no hidden budget consumption at the SM point. Confidence: high.

### m-8. Prompt-vs-YAML anchor deltas are updates, not errors (for the record)
- S_ψKs: catalog 0.710±0.011 (HFLAV Summer 2025, all-charmonium; mode-specific 0.712)
  vs prompt's 0.699 (older PDG/HFLAV). φ_s: catalog −0.041±0.016 (HFLAV PDG-2025 inputs;
  J/ψKK-only −0.050) vs prompt's −0.049. Δm_s: 17.766±0.006 ps⁻¹ (HFLAV Fall 2024
  recommended; PDG-2025 17.765 recorded as cross-check). x_D=0.405%, y_D=0.636% vs prompt
  0.4/0.62. All snapshot-pinned and fact-check-stamped in the sidecars. Confidence: high.

---

## Verified correct (explicit checks per the slice prompt)

1. **Anchors vs 2025-26 values** — ε_K^exp = 2.228(11)e-3 (PDG 2026) ✓; ε_K^SM = 2.161e-3,
   BGS2020, grouped σ {0.153, 0.076, 0.065}e-3 → 0.183e-3 total = published 0.18e-3 ✓
   (κ_ε = 0.94(2) folded into the non-perturbative group, stated and consistent ✓);
   Δm_K = 3.484e-15 GeV ✓; Δm_s: 17.766 ps⁻¹ × 6.582119569e-13 = 1.169379e-11 GeV —
   conversion exact ✓; Δm_d = 3.334e-13 GeV ↔ 0.5066 ps⁻¹ (0.07% from PDG, m-3) ✓;
   f_K = 155.7 MeV paired with the (2/3)f²m_K B ⟨O1⟩ convention (not the f/√2 convention)
   ✓; f_Bs = 230.3 MeV, f_Bd = 190.0 MeV, f_D = 212.0 MeV ✓;
   **B_K: MS-bar(2 GeV) = 0.5503 paired with Wilsons evolved to 2 GeV — the classic
   RGI-vs-μ ~30% mistake is NOT present** (B̂_K = 0.7533 is recorded in the sidecar as
   provenance only) ✓; Δm_D anchor 0.997e10 ħ s⁻¹ → 6.562e-15 GeV ✓ and C001 explicitly
   overrides the stale core 6.25e-15 default ✓.
2. **NP-room logic** — K001: documented central/tight/loose band, hard veto = loose edge,
   all arithmetic reproduced (central 6.70e-5, σ_comb 2.367e-4, loose 3.037e-4) ✓;
   B003: (residual + σ_comb)/2 in M12 space, reproduced (2.635e-13 GeV) ✓; B011:
   |exp−SM| + √(σ_exp²+σ_SM²) = 0.9e-5 + 2.55e-5 ✓ with exp 3.49(19)e-4 / SM 3.40(17)e-4
   exactly the prompt's values ✓; circularities (B002 β-σ, B012 normalization) are
   documented in docstrings/diagnostics ✓.
3. **Gating** — all eight gates: ratio = |NP|/room, passes = ratio ≤ 1, Severity.HARD ✓.
   No double-counting: Δm rooms are halved into M12 space exactly once (B003 budget /2,
   C001 Δm/2, C002 Δm/2 normalization; core Δm^NP = 2|M12^NP| identity respected) ✓.
   B002/B004 use the full complex M12 (diagnostics assert `phase_uses_complex_m12_not_abs`)
   ✓; B004 budget is direction-aware in the asymmetric SM σ ✓; C002 places |Im| (not Im,
   not |full|) against a symmetric max-|sin φ_M| room — consistent one-sided treatment ✓.
4. **Units end-to-end (K001, B003)** — C(GeV⁻²)·⟨O⟩(code GeV³ = f²m_M·#) → M12 in GeV;
   ε_K = κ_ε Im M12/(√2 Δm_K) dimensionless; B003 compares GeV to GeV. Chain is
   dimensionally consistent; the *only* normalization defect is the upstream missing
   1/(2m_M) (other slice; folded as the global 2× in I-1). ħ conversion constant
   6.582119569e-13 GeV·ps is the CODATA value, used once ✓.
5. **B011/B012 running** — LL exponents γ77⁰ = 32/3 → η^{16/23}, γ88⁰ = 28/3 → η^{14/23}
   (nf = 5), C8→C7 mixing (8/3)(η^{14/23}−η^{16/23}); threshold segmentation and the
   upper-triangular composition law verified; reproduces the textbook M_W→m_b numbers
   (u77 = 0.690 vs 0.692, u78 = 0.087 vs ≈0.09) ✓. nf selection at segment tops correct
   (6 above m_t, 5 below; m_b = 4.18 < μ_b = 4.8 so no spurious nf=4 segment) ✓.
   No C7-C7′ interference in the inclusive rate (correct to m_s/m_b) ✓.
6. **C002 factor-2 audit** — budget = sin(max|φ_M^95|) = sin(1.48°) = 0.02583 from the
   parsed HFLAV CKM2025 snapshot (all-CPV row, CI [−1.48, 1.35]°) ✓; constraint
   |Im M12^NP| ≤ (Δm_D/2)·sin φ_M uses |M12| = Δm/2 exactly once — no factor-2 error ✓.
   The "total |M12| saturated by Δm_exp/2, no SM LD phase" proxy is loudly documented
   (NEEDS-HUMAN-PHYSICS) ✓.
7. **CKM phase layer** — β = arg(−V_cd V_cb*/V_td V_tb*), β_s = arg(−V_ts V_tb*/V_cs V_cb*),
   φ_s = −2β_s: standard rephasing-invariant definitions, correct index mapping verified
   against the matrix layout; repo values sin2β = 0.7083, β_s = 1.087°, φ_s = −0.0379 rad
   are PDG-fit-consistent ✓.
8. **Adapter boundary** — `physics_adapters/deltaf2.py` and `bsgamma.py` are pure
   re-export/override shims; budget overrides preserve core |M12^NP| and only replace the
   room (checked all five wrappers); `bd_source_as_bs_source` slot-mapping for b→d keeps
   hermiticity ((1,2) and (2,1) conjugate) ✓.
9. **Tests** — 77/77 pass; YAML-anchor pinning and degrade-gracefully paths covered
   (but see M-2 for what they cannot catch).

---

## Recommended actions (for the orchestrator, no code changed by this audit)

1. After the slice-1 ME fix lands, re-run K001/B001/B003/C001/C002 floors; expect ε_K
   floors to RISE ≈1.3× (LR points) and |M12| floors to drop ≤√2 (VLL points). No budget
   changes needed (I-1).
2. Add absolute literature-pinned ME/M12 tests (M-2).
3. Unify bag-parameter scales with μ_had = 2 GeV or evolve Wilsons to each system's bag
   scale (M-1); update `bag_param_inventory.md` and the B003 docstring caveat.
4. Refresh `DELTA_M_BS_EXP`/`DELTA_M_BD_EXP` core constants to the YAML values or label
   them legacy-diagnostic-only (m-3).
5. Move `_SM_CHOICE_SENSITIVITY` into K001.yaml (m-2); fix the self-contradictory
   "smaller of" sentence in `epsilon_k_sm_decision.md` (m-1).
6. Keep B011/B012 out of any quoted physical M_KK floor table until a real C7 matching
   exists (m-5).
