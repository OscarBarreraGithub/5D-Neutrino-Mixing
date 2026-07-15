# Report 02 — Quark model core (quarkConstraints: model, fit, couplings, scan, scales, proxies, benchmarks, validation, finite_stats, diagnostics)

**Structural summary.** The Lane B pipeline is: `model.py` (spurions → C_Q/C_u/C_d → BulkMassMap → f_IR overlaps) → `fit.py` (M = 2v F_Q Y F_x, SVD, PDG rephasing, least-squares in log-SV + left-rotation chart) → `couplings.py` (mass-basis g·F² matrices) → `deltaf2.py` (KK-gluon tree Wilsons, LL RG running, hadronic MEs) → `scan.py`/`validation.py` gates. Core physics checked out well: the spurion chirality assignments (Y Y† for the doublet, Y† Y for singlets) are correct; CKM ordering V = U_Lu† U_Ld is correct; the SVD/rephasing machinery is verified numerically (five anchors pinned real ≥ 0, masses preserved, full-matrix rephasing invariance to 2e-16); the KK-gluon color/Fierz structure (C1 = g²Δ²/6M², C4 = −g²Δ_LΔ_R/M², C5 = +g²Δ_LΔ_R/3M²) is internally consistent with a direct octet-exchange derivation; hadronic inputs (Δm_K, Δm_d, Δm_s, Δm_D, GGMS matrix elements) check against PDG/HFLAV. The genuine problems found are: two logic bugs around the fit-seed canonicalization (confirmed numerically), the KK-gluon coupling normalization (documented "legacy" but a ~13× underestimate of ΔF=2 ratios vs. the physically-normalized choice), and a stale/mislabeled CKM target.

### MAJOR — Reported fit seed does not reproduce the fitted result (gauge "canonicalization" changes the physics)
- **File:** quarkConstraints/fit.py:627-655, 843
- **Category:** logic-bug
- **Claim:** `fit_quark_sector` returns `QuarkFitSolution.seed` produced by `_canonicalize_reported_seed`, which replaces `up_left → identity` and rebuilds `down_left` from the **fitted CKM observables** — but the physical invariant is the relative spurion rotation U_u†U_d, which is *not* the CKM matrix (F-factor distortions intervene), so the "canonical representative" is not gauge-equivalent to the optimum (5-decimal rounding compounds this).
- **Evidence:** `best_seed = _canonicalize_reported_seed(raw_best_seed, ...)`; re-evaluating the reported seed: fitted score `5.5e-9` vs reported-seed re-eval score `0.699`; masses_up drift from `[1.19e-3, 0.602, 162.4]` to `[5.5e-4, 0.588, 167.6]`. Common left rotation of *both* Y_u and Y_d is a gauge freedom; setting only `up_left = I` while taking `down_left` from CKM observables is not.
- **Confidence:** high (numerically confirmed)
- **Fix:** Report the raw optimized seed (or canonicalize by the true gauge transformation `Y_x → W†Y_x` with W = fitted `up_left` unitary), and never round the reported representative.

### MAJOR — Seed chaining in scans double-applies `overall_scale`, destroying the warm start
- **File:** quarkConstraints/scan.py:360-373 (also validation.py:242,266,271)
- **Category:** logic-bug
- **Claim:** `run_quark_scan` chains `current_seed = solution.seed` (whose singular values already have the scale absorbed, `overall_scale=1.0`) but then calls `fit_quark_sector(..., overall_scale=float(overall_scale))` again, so `_canonicalize_fit_seed` multiplies the already-physical singular values by 3.0 a second time; `validation.r_sweep_plot_data` does the same with `scale=2.8` computed once from the default seed.
- **Evidence:** chained template up-SVs `[1.24, 2.96, 12.40]` = 3× the reported `[0.414, 0.986, 4.13]`; initial score of the chained seed at r=0.25 is **5.43** vs 0.30 for the clean default seed — the "warm start" is ~18× worse than a cold start (compounded by the MAJOR above). With `max_nfev=120` capped, this can silently degrade or fail fits mid-sweep.
- **Confidence:** high (numerically confirmed)
- **Fix:** Pass `overall_scale=None` (or `current_seed.overall_scale`) when chaining a canonicalized seed, as `proxies.sweep_r_proxy_summary` correctly does.

### MAJOR — KK-gluon coupling uses perturbative g_s(M_KK) ≈ 1.0, omitting the √(2πkr_c) ≈ 8.5 IR enhancement, combined with M_KK ≡ Λ_IR
- **File:** quarkConstraints/couplings.py:145-147; quarkConstraints/scan.py:378 (`g_s_star=None` hardcoded, ditto validation.py:112,163,277)
- **Category:** wrong-formula / units (partially documented as "legacy repo_v1 behavior")
- **Claim:** The first KK gluon couples to IR-localized zero modes with strength ≈ g_s√(2πkr_c)·f² (≈ 8.5·g_s here, per the module's own docstring citing CFW), but every production-path call passes `g_s_star=None`, so Wilson coefficients are suppressed ×(1/8.48)² ≈ 0.014; simultaneously M_KK = Λ_IR instead of the physical 2.449·Λ_IR inflates 1/M² by 6.0. Net: ΔF=2 ratios in this module are ~13× too small at fixed physical spectrum.
- **Evidence:** benchmark at Λ_IR = 3 TeV: ε_K ratio-to-budget = **20.7** (legacy) vs **269.8** (g_s* = 8.48·g_s, M_KK = 2.449Λ) — floor understated ~3.6× in scale. `scales.py` defines `GAUGE_KK_ROOT_NN = 2.4487` but `DEFAULT_QUARK_XI_KK = 1.0`; FLOOR_SUMMARY.md states all quoted floors are *physical* M_KK = 2.4487·Λ_IR, so any floor derived from this module's `M_KK` column conflates the two conventions.
- **Confidence:** high on the physics; medium on downstream impact (the quoted ~7 TeV Lane-B floor comes from the `flavor_catalog` K001 path, not audited here — but this module's caveat is documented only as "bookkeeping", not as a ~13× amplitude effect).
- **Fix:** Default `g_s_star = g_s·√(2πkr_c)` and `xi_KK = GAUGE_KK_ROOT_NN`, or prominently tag every output of this lane as non-physical bookkeeping units.

### MINOR — CKM target is not "PDG 2024 central" as documented; scan tolerances derived from different values than the target
- **File:** quarkConstraints/benchmarks.py:32-34 (target), quarkConstraints/scan.py:40-51 (tolerances)
- **Category:** stale-data / convention-inconsistency
- **Claim:** `_TARGET_CKM` (θ12=0.2274, θ13=0.00368, θ23=0.0415, δ=1.196) yields |V_us|=0.22544, |V_cb|=0.04149, |V_ub|=0.00368, J=3.118e-5, while `default_quark_targets()`'s docstring calls it "the PDG 2024 central unitary" and scan.py derives 2σ tolerances from PDG 2024 values (0.22501, 0.04183, 0.00382, 3.08e-5); δ=1.196 is a PDG ~2018-era value (PDG 2024: ≈1.147).
- **Evidence:** |V_us| target is +0.9σ_PDG off center; |V_ub| −3.7% (~1.4σ); J +1.2%. Acceptance bands are centered on the stale target, not on the PDG 2024 central used to size them.
- **Confidence:** high
- **Fix:** Regenerate `_TARGET_CKM` angles from the same PDG 2024 global-fit numbers used in the scan.py tolerance comments.

### MINOR — Wilson upper limit default z=1.92 conflates ΔlnL=1.92 with a Gaussian z-score
- **File:** quarkConstraints/finite_stats.py:8
- **Category:** statistics
- **Claim:** `wilson_upper_limit(..., z=1.92)` gives a one-sided ~97.3% (two-sided 94.5%) interval and a k=0 limit of 3.69/n, whereas the standard 95% one-sided upper limit is z=1.645 (2.71/n Wilson) or the exact "rule of three" 3.0/n; 1.92 is the half-χ² log-likelihood threshold (3.84/2), not a z-score. The Wilson formula itself is algebraically correct; tests (tests/test_finite_stats.py:36-56) pin the convention rather than a stated CL.
- **Confidence:** high on the numbers, medium on intent
- **Fix:** Use z=1.645 for a 95% one-sided limit (or document the actual coverage, 97.3% one-sided, wherever "95%" is quoted).

### MINOR — Yukawa-normalization: mass prefactor 2v with O(3) spurions implies Ȳ = 2·(spurion) ≈ 6–8, outside the repo's |Ȳ| < 4 perturbativity convention
- **File:** quarkConstraints/fit.py:367 (`prefactor = 2.0 * v`), benchmarks.py:104/218 (`overall_scale` 2.8–3.0)
- **Category:** convention-inconsistency (factor-of-2)
- **Claim:** With m = 2v·F·Y·F the spurion is Y = k·Y_5D = Ȳ/2 in the repo convention (CLAUDE.md: m = 2vk f Y_5D f, Ȳ = 2kY_5D); the fitted top spurion entry is ~3.5–4.1 (overall_scale ≈ 3 × sv₃ ≈ 1.24), i.e. Ȳ_t ≈ 7–8, well past the |Ȳ| < 4 perturbativity gate used elsewhere in the repo — no perturbativity check exists anywhere in the quark fit/scan path.
- **Evidence:** fitted physical up SVs `[0.414, 0.986, 4.13]`; either the prefactor should be v (if the spurion is Ȳ) or the scan needs a |2Y| < 4 gate.
- **Confidence:** medium (depends on which normalization the quark lane intends; nothing in the package pins it).
- **Fix:** Document the quark-spurion normalization and add a perturbativity gate consistent with it.

### NOTE — f_IR "c = 0.5" plateau is ~5e-6 wide, flattening fitted c values near 1/2
- **File:** warpConfig/wavefuncs.py:35 (`~np.isclose(c_arr, 0.5)`)
- **Category:** numerics
- **Claim:** `np.isclose` default (rtol=1e-5) treats |c−0.5| < ~5.01e-6 as exactly 0.5, so f_IR is piecewise-constant over a 1e-5-wide plateau (rel. step ~9e-5 at the edges) — a flat spot/derivative kink inside the `least_squares` search region (BulkMassMap range [0.30, 0.72] crosses 0.5; benchmark c_u2 = 0.525, c_Q2 = 0.561 are nearby).
- **Evidence:** `f_IR` constant = 0.11794563 for c ∈ [0.5−5e-6, 0.5+5e-6]; exact values 0.11795623 / 0.11793503 at the edges.
- **Confidence:** high (behavior), low (practical impact)
- **Fix:** Use `np.abs(c-0.5) < 1e-12` (or a series expansion) instead of `np.isclose` defaults.

### NOTE — Legacy-fallback D-mixing operator weights understate the chiral LR enhancement (dead path in production)
- **File:** quarkConstraints/deltaf2.py:44,251-264 (`lr1_weight=7.0` universal)
- **Category:** wrong-formula (legacy fallback only)
- **Claim:** The `use_hadronic=False` surrogate applies the *kaon* LR enhancement weight (7.0) to all systems including D and B, whereas the true chiral ratios differ by system; production uses the hadronic path so this only affects anyone reviving the fallback.
- **Confidence:** medium
- **Fix:** Either delete the fallback or give per-system weights.

### NOTE — `diagnostics.extract_msbar_masses_from_yukawa_row` hardcodes the fitter scale to m_t(m_t) despite documenting a legacy 3 TeV layout
- **File:** quarkConstraints/diagnostics.py:128-131
- **Category:** logic-bug (latent)
- **Claim:** The docstring says CSV Yukawas live "at m_t(m_t) … or mu = 3 TeV in the legacy layout", but `fitter_scale = M_TOP_MS` unconditionally — auditing a legacy-layout CSV would silently mis-run all six masses from the wrong reference scale.
- **Confidence:** high (code), medium (whether legacy CSVs are still audited)
- **Fix:** Accept a `fitter_scale_GeV` argument keyed to the CSV layout.

## Cross-module convention inconsistencies
1. **M_KK vs Λ_IR:** `scales.py` defines the physical gauge root 2.4487 but defaults `xi_KK=1.0`; `couplings.py`/`deltaf2.py`/`proxies.py` all emit columns literally named `M_KK` that are actually Λ_IR, while `docs/FLOOR_SUMMARY.md` declares all quoted floors to be physical M_KK = 2.4487·Λ_IR. Any consumer joining scan CSVs with the docs mixes conventions by ×2.45.
2. **g_s vs g_s\*:** `couplings.py` documents the √(2πkr_c) enhancement, yet all three call sites hardcode `g_s_star=None`; combined with item 1 this makes the module's ΔF=2 amplitudes ~13× below the physically-normalized choice.
3. **Yukawa normalization:** lepton stack Ȳ=2kY with |Ȳ|<4 checks; quark stack `prefactor = 2v` plus O(3) spurions implies Ȳ up to ~8 with no gate.
4. **Seed "canonical chart" vs physics:** fit.py's quotient chart treats (up_left, down_left) → (I, CKM-derived) as gauge, but only a *common* left rotation is gauge; leaks into scan.py and validation.py chaining (the two MAJOR logic bugs).
5. **BulkMassMap surrogate:** MODEL_CONVENTIONS §3 decides for the affine eigenvalue map, while model.py/production still route through the saturating squash — documented TODO (not re-flagged).

**Known issues confirmed, not rediscovered:** KI-1 (B002/B004), FLAG B4/B5 scale caveat, B2 rephasing (independently re-verified correct), BulkMassMap surrogate caveat.

**Verified correct (spot-checked numerically):** spurion chirality structure (model.py:235-237), CKM ordering and PDG rephasing (fit.py:273-359), Jarlskog quartet sign, KK-gluon color/Fierz coefficients {1/6, −1, +1/3}, mass-basis rotation of overlaps U†F²U, hadronic constants, GGMS M12-normalization, Wilson-score algebra in finite_stats.py.

**Not reviewed (out of budget):** qcd_running.py internals, deltaf2.py RG basis map beyond structure, pdg_quark_masses.py values, most tests, paper_0710_1869/, modern/, flavor_catalog K001 path.
