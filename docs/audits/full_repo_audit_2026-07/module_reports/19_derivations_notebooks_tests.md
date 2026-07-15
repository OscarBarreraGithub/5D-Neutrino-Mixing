# Report 19 — Gap-closing: derivations/, notebooks, remaining scripts, test expected-values spot audit

**Structural summary.** derivations/ contains only two real .tex files (conventions.tex, rs_ew_gauge_kk_coupling.tex, plus kk_modes/bessel_equation.tex); the mass-formula derivations are explicitly marked NEEDED in derivations/masses/README.md. The RG-evolution test pins were reproduced independently to 12 digits — physically correct, not self-pins. The GGMS matrix-element oracle test is genuinely independent. The dominant test-suite weakness is ~105 long-float default-tolerance pytest.approx snapshot pins across ~40 files (self-oracle pattern). One notebook still headlines the known-buggy pre-B1 "25–30 TeV Z→bb floor" with no superseded banner.

## PRIORITY 1 VERDICT — Seesaw prefactor: **UNDERIVED, but code is internally self-consistent (not a bug)**

- No seesaw derivation exists: derivations/masses/README.md:13-14 marks charged_leptons.tex and seesaw.tex NEEDED, and lines 44-45/73-74 list "Where does the factor of 2 come from?" and "Why (f_N^UV)⁻²?" as open questions. conventions.tex:146,154 restates both formulas as boxed *conventions* matching the code verbatim — derives nothing.
- **However the code's 2 is NOT a naive error.** The one real derivation, derivations/rs_ew_gauge_kk_coupling.tex:95-96, derives brane values of the normalized fermion zero mode: **g₀(c,1)² = 2·f_IR²(c)** and **g₀(c,ε)² = 2·f_UV²(c)**. The brane factor 2 applies to *both* brane bilinears: m_D = 2vk·f_L·Y·f_N (IR Yukawa) AND M_eff = 2·M_N·(f_N^UV)² (UV Majorana). Then m_ν = m_D²/M_eff = 4k²v²f²f²Y²/(2M_N f_UV²) = **2k²v²f_L²f_N²Y²/((f_N^UV)²M_N)** — exactly the code. The "naive 4" forgets the same brane-2 in the Majorana term. Numerical sanity check lands m_ν ≈ 0.008 eV at the benchmark — correct ballpark.
- No validation exists either: masses/README.md:86 references neutrinos/seesawValidation.py (full 6×6 Takagi cross-check) which does not exist; no test diagonalizes the full neutral mass matrix.
- Charged-lepton 2vk and f_IR/f_UV: conventions.tex matches warpConfig/wavefuncs.py and yukawa/charged_lepton.py exactly; supported by the same g₀²=2f² result; assembled derivation likewise NEEDED. No derivation-vs-code contradiction found anywhere.

### [MAJOR] Notebook headlines the known-buggy pre-B1 "25–30 TeV Z→bb floor" with no superseded banner
- **File:** notebooks/wq_quarkonly_explore.ipynb (cells 0, 2, 34) + notebooks/wq_quarkonly_explore.py
- **Category:** doc-code-mismatch / stale-data
- **Claim:** Headline finding — "Z→bb (T010/T011) rigorous and vetoing... pushing the rigorous M_KK floor to ~25–30 TeV (strict survival 0% below 20 TeV)" — is exactly the B1 bug the repo corrected; no SUPERSEDED banner (unlike the .tex notes).
- **Confidence:** high
- **Fix:** Banner or delete/rerun.

### [MINOR] Derivation READMEs claim .tex files DONE that do not exist
- **File:** derivations/zero_modes/README.md:13, flavor/README.md:13, masses/README.md:15
- **Claim:** f_factors.tex, pmns_convention.tex, neutrino_spectrum.tex all tagged DONE but the directories contain only README.md; geometry/ marked DONE with no .tex at all.
- **Fix:** Commit the missing .tex or retag NEEDED.

### [MINOR] ~105 long-float default-tolerance snapshot pins ("self-oracle" tests) beyond the known offenders
- **File:** tests/ repo-wide grep
- **Category:** circular-validation / test-bug
- **Claim:** 105 pytest.approx(<12+-digit float>) assertions pin adapter outputs rather than independent oracles. Worst offenders beyond the known list: test_rs_ew_custodial_pr2.py (8), test_rs_ew_phase6a_zbb_fermion_mixing.py (6), test_mu_to_e_gamma.py:141-146 (5 — pure pipeline snapshots), test_B002.py (5), test_ckm_extraction.py (4), test_B003.py (4); test_quark_deltaf2.py:170 is an explicitly acknowledged re-pinned snapshot (20.667539197050676).
- **Confidence:** high (pattern), medium (severity)
- **Fix:** Tag snapshot pins regression-only; ensure each constraint also has an independent-oracle or scaling test.

### [NOTE] Cross-lane quark mass-formula prefactor differs by 2 (v·F·Y·F anarchic vs 2v·F·Y·F fitter)
- **File:** scripts/run_rs_anarchy.py:446 (and anarchic_bauer_s1_zbb.py:214) vs quarkConstraints/fit.py:367
- **Claim:** Lane A builds M_q = v·F·Y·F while Lane B fits with 2v·F·Y·F; each internally consistent, but Yukawa-prior scales/perturbativity statements not comparable across lanes without the 2. The lepton convention (m = v·f·Ȳ·f in dimensionless Ȳ) matches the *anarchic* form, not the fitter's. Likely deliberate (Bauer matching) but undocumented — and directly relevant to the Report 12 Bauer-bridge factor-2 finding.
- **Confidence:** medium
- **Fix:** Document per-lane mass prefactor in MODEL_CONVENTIONS.md.

### [NOTE] Seesaw factor-2 and (f_UV)⁻² are unvalidated as well as underived — the promised Takagi cross-check (neutrinos/seesawValidation.py) is absent; add a one-generation 2×2 Takagi test pinning the prefactor.

## Verified correct
- test_wilson_rg_audit.py pins reproduced independently to 12 digits (C1=0.729130912, C4=3.538163975, C5→C4=0.894757449, C5=0.853891628).
- test_epsilon_k_physics.py GGMS oracle genuinely independent (rationals in-test); ratio and M12 normalization structurally correct.
- test_diagnostics.py:45 prefactor exactly matches the fitter — no test-vs-code mismatch (the factor-2 issue is cross-lane).
- test_quark_model.py structural-sound; test_contract.py meaningful engineering contract (no physics, by design).
- _audit_followup_crossings.py, c02c_post_process.py, rs_anarchy_mkk_min_hist.py: M_KK^min = M_KK·√(max_ratio) and g_s*=3 rescale both correct; rs_anarchy_gate_sensitivity.py log-interpolated crossings correct.
- wq_quarkonly_1m_plan.py seed-disjointness correct; analyze_wq_quarkonly.py reservoir sampling (Algorithm R) and veto accounting correct; anarchic_bauer_s1_zbb.py internally consistent (Z99=2.5758 correct).
- qcd/alphaS.ipynb and quark_benchmark_validation.ipynb clean (correct β coefficients).
- conventions.tex f_IR/f_UV/BC formulas match code exactly.

## Not reviewed
build_wq_quarkonly_comparison.py beyond structure; ~140 of 181 test files (targeted set + pin grep only); notebook output cells and remaining large notebooks beyond stale-claim grep (mkk_vs_yukawa C=0.02 hits are deliberate 2008-comparisons); rs_ew_gauge_kk_coupling.tex line-by-line beyond the boxed results; kk_modes/bessel_equation.tex; plot_rs_anarchy_summary.py + hist variants beyond shared-pattern check.
