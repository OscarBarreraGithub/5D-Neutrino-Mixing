# Minor Findings Ledger

Scope: every minor finding and note in `AUDIT_COMPENDIUM.md` section 4 and section 8.4 is listed below. Items marked DOC are documented here as known limitations or caveats for a later physics or validation pass. The additional section 8.4 paragraph headed "Additionally verified correct" is not counted as a defect list.

Status totals: FIXED 8, DOC 50, WONT-FIX 2, ALREADY-COVERED 0, ALREADY-CORRECT 1, ESCALATE 2.

ESCALATE items (RESOLVED in escalation cycle, commit 422dc18, Claude audit APPROVE):

- RD-01: "Z-prime proxy C9/C10 matching is off by x2 and sign." FIXED: corrected by -0.5 (1/2 chiral
  projector + tree-level exchange sign; WET/Buras-Girrbach) in rare_b_dilepton, rare_charm_dilepton,
  rs_semileptonic_wilsons. Symmetric LFU ratios unchanged; absolute NP shifts + interference move.
- NMF-03: "T-series silently pass malformed coupling objects." RESOLVED: verified T010/11/12/14/15/16/17
  already fail-close on malformed z_delta_g_* via M-18 (_coupling_entry -> invalid_extra), executed NaN
  cases directly; added regression coverage. (Also fixed a related M-18 over-reach straggler where L-series
  legitimate unevaluated cases were wrongly fail-closed: L002/L004/L005/L007/L008/L010 now follow the L009
  pattern; malformed still fails closed.)

## Section 4: Lepton Core, Report 01

- LC-01, flagship benchmark nonperturbative: "flagship benchmark point `c_E=[0.75,0.60,0.50]` fails its own `|Ybar|<4` bound." DISPOSITION: DOC. Rationale and files: documented as a known nonperturbative example needing benchmark refresh, not a safe numerical edit. Files touched: this ledger.
- LC-02, `_F_exact` tiny-positive sign: "`_F_exact` fake-positive sign at `x<=1e-12`." DISPOSITION: DOC. Rationale and files: real low-level numerical cleanup, but all report evidence says it is below current physics relevance and needs targeted numerical tests. Files touched: this ledger.
- LC-03, fallback scan direction: "fallback root scan only searches upward of last seed." DISPOSITION: DOC. Rationale and files: fallback robustness issue, not an unambiguous one-line fix without defining a new scan contract. Files touched: this ledger.
- LC-04, tau mass label: "`M_TAU=1.77686` predates the 2023 PDG fit, but is labelled PDG 2024." DISPOSITION: FIXED. Rationale and files: label corrected to a legacy PDG 2022 retained value while keeping the valid numeric value. Files touched: `yukawa/constants.py`, `yukawa/charged_lepton.py`.
- LC-05, NuFIT pin: "NuFIT 6.1 inputs are pinned before public validation." DISPOSITION: DOC. Rationale and files: needs a data refresh after the public release is verified, not a safe local correction. Files touched: this ledger.
- LC-06, naturalness massless-neutrino gate: "naturalness filter `min|Ybar|>=0.1` spans both sectors, so exactly massless-neutrino options cannot pass." DISPOSITION: DOC. Rationale and files: physics policy choice for scan selection, needs a deliberate model-scope decision. Files touched: this ledger.
- LC-07, Takagi branch-cut accuracy: "Takagi `sqrtm` branch-cut accuracy loss near singular spectra." DISPOSITION: WONT-FIX. Rationale and files: report marks the loss bounded and benign for current scans, so no code change is justified in this minor pass. Files touched: none.
- LC-08, seesaw prefactor derivation: "seesaw prefactor 2 vs 4 internally consistent but underived." DISPOSITION: DOC. Rationale and files: derivation note is needed before changing formulas because the code and tests are internally consistent. Files touched: this ledger.

## Section 4: Quark Model, Report 02

- QM-01, CKM target vintage: "CKM fit target not PDG 2024 central." DISPOSITION: DOC. Rationale and files: retuning CKM targets changes fit baselines and needs a dedicated calibration pass. Files touched: this ledger.
- QM-02, z equals 1.92 label: "`wilson_upper_limit` default `z=1.92` is mislabeled as 95 percent one-sided or conflated." DISPOSITION: FIXED. Rationale and files: comments and figure labels now call this the z=1.92 audit convention, about 97.3 percent one-sided coverage, rather than standard 95 percent. Files touched: `quarkConstraints/finite_stats.py`, `scripts/rs_anarchy_mkk_min_hist_by_cvals.py`, `scripts/rs_anarchy_cfw_comparison.py`, `docs/audits/zero_pass_inventory.md`.
- QM-03, `np.isclose` c plateau: "`np.isclose` `c=0.5` plateau creates finite coupling at exactly flat profiles." DISPOSITION: DOC. Rationale and files: boundary convention affects profile numerics and needs analytic expectation plus regression tests. Files touched: this ledger.
- QM-04, legacy hadronic fallback: "legacy `use_hadronic=False` fallback uses kaon LR weights for all systems." DISPOSITION: DOC. Rationale and files: legacy path should either be retired or rebuilt with per-meson weights, which is not a safe minor edit. Files touched: this ledger.
- QM-05, diagnostics scale layout: "`diagnostics.py` fitter scale is hard-coded to `m_t(m_t)` despite legacy 3 TeV layout." DISPOSITION: DOC. Rationale and files: changing the reported scale affects fit comparability and should be handled with a dedicated diagnostics migration. Files touched: this ledger.

## Section 4: QCD, Report 03

- QCD-01, alpha_s package split: "`alpha_s(M_Z)` mismatch: standalone `qcd` package uses 0.1180, DeltaF2 wrapper tests use 0.1179." DISPOSITION: FIXED. Rationale and files: DeltaF2-facing defaults and tests now import the shared PDG 2024 `ALPHA_S_MZ=0.1180`; stale docs and RG snapshots were aligned. Files touched: `quarkConstraints/qcd_running.py`, `quarkConstraints/bsgamma.py`, `tests/test_qcd_running.py`, `tests/test_wilson_rg_audit.py`, `tests/constraints/primary/kaon/test_K001.py`, `tests/constraints/primary/kaon/test_K002.py`, `docs/STATE_OF_PROJECT.md`, `docs/quark_scan_assumptions_compact.tex`, `review_local/deltaF2_framework_review.tex`.
- QCD-02, mass-running scheme anchor: "`run_msbar_mass` anchors alpha_s scheme ignoring `n_f_ref`." DISPOSITION: DOC. Rationale and files: scheme contract needs a real API decision and validation across thresholds. Files touched: this ledger.
- QCD-03, top mass stale constant: "`M_TOP_POLE=172.69` is stale." DISPOSITION: FIXED. Rationale and files: reference constant and notebook label updated to PDG 2024 value 172.57 GeV, with the code comment clarifying it is the direct-reconstruction top mass retained for reference. Files touched: `qcd/constants.py`, `qcd/alphaS.ipynb`.
- QCD-04, light-quark sigma semantics: "light-quark sigma symmetrization inconsistent with docstring." DISPOSITION: DOC. Rationale and files: uncertainty semantics should be corrected with consumers checked, not changed silently. Files touched: this ledger.

## Section 4: DeltaF2, Report 04

- DF2-01, VIA additive terms: "VIA additive terms likely double-count." DISPOSITION: DOC. Rationale and files: needs a physics derivation and a coordinated M12 recalibration, not a quick arithmetic edit. Files touched: this ledger.
- DF2-02, stale D-mixing constant: "`Delta m_D` experimental constant is stale." DISPOSITION: DOC. Rationale and files: catalog override exists, but core constants need a dedicated synchronization pass. Files touched: this ledger.
- DF2-03, epsilon_K precision: "`epsilon_K^SM` has an extra unpublished digit." DISPOSITION: DOC. Rationale and files: precision display should be harmonized with source provenance before changing tests or snapshots. Files touched: this ledger.
- DF2-04, Delta m_K default bundle: "`Delta m_K` evaluator is not in the default bundle." DISPOSITION: DOC. Rationale and files: adding the evaluator changes production constraint coverage and needs a deliberate scan-policy change. Files touched: this ledger.
- DF2-05, kaon bound surrogate: "`DeltaF2ObservableSummary.bound` is legacy surrogate for kaons while `ratio_to_bound` uses hadronic budget." DISPOSITION: DOC. Rationale and files: object semantics need an API cleanup so callers do not mix display and budget bounds. Files touched: this ledger.

## Section 4: RS-EW, Report 05

- RSEW-01, custodial model alias: "`custodial_rs_su2r` silently relabelled `custodial_rs_plr`." DISPOSITION: DOC. Rationale and files: alias provenance should be made explicit in a compatibility note before altering model names. Files touched: this ledger.
- RSEW-02, oblique L hardwire: "oblique proxy hardwires L=35 while geometry gives 35.9 to 37.5." DISPOSITION: DOC. Rationale and files: weak-effect proxy drift needs a geometry-aware recalibration. Files touched: this ledger.
- RSEW-03, c_S coefficient: "`c_S=30` geometric coefficient may be about 20 percent low." DISPOSITION: DOC. Rationale and files: report marks confidence low, so this is a derivation item rather than a safe fix. Files touched: this ledger.

## Section 4: Rare Decays, Report 06

- RD-01, Z-prime C9/C10 convention: "Z-prime proxy C9/C10 matching is off by x2 and sign relative to standard weak-effective-theory conventions." DISPOSITION: ESCALATE. Rationale and files: this can rescale many rare-decay constraints and may affect interference signs, so it needs its own fix cycle with convention tests. Files touched: this ledger.
- RD-02, lepton axial sign split: "kaon vs B/charm lepton axial sign opposite." DISPOSITION: DOC. Rationale and files: could be convention-canceling or physical depending on operator definitions, so it needs a sign-convention audit. Files touched: this ledger.
- RD-03, B to tau nu silent defaults: "B to tau nu silent defaults for Vub, decay constants rely on YAML overrides." DISPOSITION: DOC. Rationale and files: default provenance should be refactored with adapter tests, but current production path relies on explicit YAML values. Files touched: this ledger.

## Section 4: Catalog Kaon and Charm, Report 10

- KCH-01, K008 ISU coefficients: "K008 hardcoded ISU coefficients mismatch YAML SM limit." DISPOSITION: DOC. Rationale and files: resolving this requires source-level coefficient reconciliation, not a safe label edit. Files touched: this ledger.
- KCH-02, K005 Grossman-Nir limit: "K005 uses KOTO limit while NA62-implied Grossman-Nir stronger." DISPOSITION: DOC. Rationale and files: this is a constraint-policy choice that changes exclusions and should be handled as a catalog update. Files touched: this ledger.
- KCH-03, C001 stale core override: "C001 catalog overrides stale core Delta m_D but core users do not." DISPOSITION: DOC. Rationale and files: core and catalog constants need synchronization in a dedicated constants pass. Files touched: this ledger.
- KCH-04, C003 direct CP non-veto: "C003 direct CP observation non-vetoing." DISPOSITION: DOC. Rationale and files: currently an honest INFO-style catalog limitation, but turning it into a veto needs a defined statistical model. Files touched: this ledger.

## Section 4: Catalog EW, Collider, and Lepton, Report 11

- ECL-01, point builder M_KK name: "point_builder exports bookkeeping M_KK as `kk_gluon_mass_gev`." DISPOSITION: DOC. Rationale and files: naming cleanup can break downstream consumers and needs a compatibility field plan. Files touched: this ledger.
- ECL-02, Gfitter year label: "Gfitter year 2026 anchors are 2018." DISPOSITION: FIXED. Rationale and files: remaining EW001 source manifest year corrected to 2018; the main EW001 YAML already carried the 2018 source text. Files touched: `flavor_catalog/references/EW001/source_manifest.yaml`.
- ECL-03, T010 and T011 gate semantics: "T010/T011 gate direction-blind undefined CL." DISPOSITION: DOC. Rationale and files: needs a catalog CL and direction convention before tightening veto behavior. Files touched: this ledger.
- ECL-04, collider sigma times BR path: "collider sigma x BR path exists but is never populated; all CR reduce to mass cuts." DISPOSITION: DOC. Rationale and files: populating production cross sections changes the collider model and needs data sources plus tests. Files touched: this ledger.
- ECL-05, missing global constraints: "missing Sigma mnu, 0nubb, g-2 constraints, hadronic EDMs INFO-only." DISPOSITION: DOC. Rationale and files: documented as scope gaps requiring new physics adapters and source choices. Files touched: this ledger.
- ECL-06, L over E anchors: "L/E anchors verified current." DISPOSITION: ALREADY-CORRECT. Rationale and files: reports rechecked the anchors and found no defect to fix. Files touched: none.
- ECL-07, muN to eN material mismatch: "muN->eN Al budget borrows Au limit." DISPOSITION: DOC. Rationale and files: replacing this needs a nucleus-specific conversion budget and source validation. Files touched: this ledger.

## Section 4: Scripts and Pipeline, Report 12

- SP-01, empty catalog index stub: "`flavor_catalog/catalog_index.yaml` is an empty stub while the registry has 103 constraints." DISPOSITION: FIXED. Rationale and files: file now clearly states it is a deprecated scaffold and that the registry is the authoritative source. Files touched: `flavor_catalog/catalog_index.yaml`.
- SP-02, stale temp double count: "`build_constraint_matrix.py` can double-count stale `.tmp` if a prior run crashes." DISPOSITION: DOC. Rationale and files: cleanup policy should be implemented with atomic output handling and tests, not ad hoc deletion. Files touched: this ledger.
- SP-03, collider cut mismatch: "5.5 TeV `pass_COLLIDER` cut contradicts documented 4 TeV." DISPOSITION: DOC. Rationale and files: threshold choice changes scan conclusions and needs a policy update rather than a silent constant change. Files touched: this ledger.
- SP-04, little-RS label: "S3 little RS `L=7.0` labelled `ln(1e3)=6.908`." DISPOSITION: FIXED. Rationale and files: label now says `L=7.0`, keeping the existing numerical benchmark unchanged. Files touched: `scripts/anarchic_bauer_s1.py`.
- SP-05, EW003 allowlist: "EW003 appears in quark-only allowlist, but requires extra forbidden." DISPOSITION: DOC. Rationale and files: allowlist semantics need cleanup with scan expectations and coverage tests. Files touched: this ledger.
- SP-06, fit_success default: "`fit_success` defaults True when diagnostics missing." DISPOSITION: DOC. Rationale and files: changing the default can invalidate historical scans and needs a migration policy. Files touched: this ledger.
- SP-07, survival naming inversion: "`survives_all_HARD_strict` and inclusive names are inverted." DISPOSITION: FIXED. Rationale and files: code and project docs now clarify that strict is the historical rigorous-HARD-only field and inclusive also counts proxy or partial HARD vetoes. Files touched: `scripts/run_full_catalog_scan.py`, `docs/STATE_OF_PROJECT.md`.

## Section 8.4: New Minor and Note Findings

- NMF-01, LFV phase policy split: "Opposite unknown-phase envelope policies between muN->eN destructive and mu->3e/tau constructive." DISPOSITION: DOC. Rationale and files: phase treatment needs a common policy choice and validation, not a safe formula edit. Files touched: this ledger.
- NMF-02, K009 vector and axial compression: "K009 collapses ISU vector/axial direct-CP weights to one coefficient." DISPOSITION: DOC. Rationale and files: requires source reconstruction and adapter tests for the separated weights. Files touched: this ledger.
- NMF-03, malformed T-series extras: "T-series constraints silently pass malformed coupling objects while CR-series constraints fail closed." DISPOSITION: ESCALATE. Rationale and files: malformed inputs can make HARD constraints pass silently, so this needs fail-closed behavior and targeted tests in its own fix cycle. Files touched: this ledger.
- NMF-04, VLQ mass proxy: "VLQ pair searches use `m_VLQ = M_KK`." DISPOSITION: DOC. Rationale and files: lighter custodial partner modeling requires model-dependent mass maps and cannot be safely inferred here. Files touched: this ledger.
- NMF-05, CR009 CI convention: "CR009 pools destructive contact-interaction endpoints as one-sided and ignores CI `g^2=4pi`." DISPOSITION: DOC. Rationale and files: needs an explicit contact-interaction convention and per-endpoint treatment. Files touched: this ledger.
- NMF-06, K_L LFV projection: "K_L LFV modes omit K_L CP projection and K0->pi0 isospin Clebsch." DISPOSITION: DOC. Rationale and files: current impact is zero for shipped branches, but any live K_L LFV use needs the projection and Clebsch implemented first. Files touched: this ledger.
- NMF-07, exact xi_KK rescale: "`exact` xi_KK rescale drops RG start-scale drift." DISPOSITION: DOC. Rationale and files: correcting this requires coordinated RG start-scale treatment and validation. Files touched: this ledger.
- NMF-08, CFW fixed g_s comparison: "CFW comparison fixed `g_s=1.05` while run varies." DISPOSITION: DOC. Rationale and files: figure comparison caveat, needs a plotting update tied to scan metadata. Files touched: this ledger.
- NMF-09, figure filters and non-finite ratios: "three figure scripts disagree fit-quality filters and map non-finite ratios to passing." DISPOSITION: DOC. Rationale and files: visualization pipeline needs a shared fit-quality helper and non-finite policy, but it is outside safe physics constants edits. Files touched: this ledger.
- NMF-10, phase0 tolerances: "`calibrate_phase0` tolerances from achieved residuals not PDG errors." DISPOSITION: DOC. Rationale and files: calibration tolerance definition should be changed only with a benchmark recalibration plan. Files touched: this ledger.
- NMF-11, Lane A and B mass prefactors: "Lane A quark masses prefactor `v` while Lane B `2v`, undocumented and relevant to Bauer bridge M-11." DISPOSITION: DOC. Rationale and files: physics convention needs derivation and bridge tests before any formula change. Files touched: this ledger.
- NMF-12, derivation README status: "derivation READMEs tag never-written tex as DONE." DISPOSITION: DOC. Rationale and files: documentation bookkeeping should be corrected with the missing derivation artifacts identified. Files touched: this ledger.
- NMF-13, long-float snapshot pins: "about 105 long-float snapshot pins in tests." DISPOSITION: DOC. Rationale and files: snapshot looseness needs a test-maintenance pass so precision changes remain detectable. Files touched: this ledger.
- NMF-14, K lifetime field name: "`kplus_lifetime_s` field holds K_L lifetime." DISPOSITION: DOC. Rationale and files: existing adapter diagnostics already document the legacy `kplus_*` field names carrying K_L values for K_L modes; no rename was made because it would break K020 compatibility. Files touched: this ledger.
- NMF-15, B012 one-sigma budget: "B012 budget bare `1 sigma_exp` with SM=measurement." DISPOSITION: DOC. Rationale and files: needs a statistical treatment of SM uncertainty and measurement central value before changing veto strength. Files touched: this ledger.
- NMF-16, anchor loader monkeypatching: "anchor loaders monkeypatch module-global scaffold, not thread safe." DISPOSITION: WONT-FIX. Rationale and files: current loaders are serial and test-scoped; refactor only becomes necessary if concurrent loaders are introduced. Files touched: none.
- NMF-17, MEG II comparison factor: "`meg2_constraint_comparison` text says 100x but should be 80x." DISPOSITION: DOC. Rationale and files: documentation text correction is low risk, but the exact comparison context should be checked against the current script output before editing. Files touched: this ledger.

## Verification

- Requested command: `python -m pytest tests/ -k "pdg_quark or qcd or catalog_index or lepton or constants or alpha_s" -q 2>&1 | tail -20`. Result content: 12 failed, 244 passed, 1547 deselected. Remaining failures are charged-lepton flavor-override and legacy-overlap tests where `ConstraintResult` fail-closed HARD-extra policy sets `diagnostics["evaluated"] = True` for `invalid_extra`; this was not changed because it is outside the safe minor-label and constants pass.
- Direct constants and QCD selection: `python -m pytest tests/ -k "pdg_quark or qcd or catalog_index or constants or alpha_s" -q`. Result: 68 passed, 1735 deselected.
- Direct alpha_s and RG tests: `python -m pytest tests/test_qcd_running.py tests/test_wilson_rg_audit.py -q`. Result: 30 passed.
- Direct kaon RG snapshot tests touched by the alpha_s unification: `python -m pytest tests/constraints/primary/kaon/test_K001.py::test_reference_couplings_show_qcd_running_enhancement tests/constraints/primary/kaon/test_K002.py::test_reference_couplings_show_qcd_running_enhancement -q`. Result: 2 passed.
