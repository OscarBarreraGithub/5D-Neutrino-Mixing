# Slice 5 audit — scan harness & conventions

Auditor: Claude (skeptical code+physics review), 2026-06-10. READ-ONLY on repo code; python used for data verification against raw scan artifacts.

Scope: `scripts/run_full_catalog_scan.py`, `scripts/wq_quarkonly_1m_plan.py`, `scripts/build_wq_quarkonly_comparison.py`, `flavor_catalog/website/scripts/build_scan_explorer.py`, `scripts/analyze_wq_quarkonly.py`, the two 1M analysis reports, `tests/test_full_catalog_scan_harness.py`, sbatch launchers, and prose claims in `docs/STATE_OF_PROJECT.md` / `.orchestration/PHASE2_PROGRAM_LEDGER.md` / `docs/quark_scan_constraint_update_2026-06.md`.

Cross-slice context taken as given (not re-derived): (a) T010's minimal-model Zbb fermion-mixing term is wrong (sign + ~190x); (b) ΔF=2 matrix elements carry a basis swap + factor-2 making ε_K ~1.7x too weak; (c) EW001 ΔT uses a Λ_IR-convention coefficient with physical M_KK (~6x weak).

---

## Headline data facts established from raw artifacts (used throughout)

Verified by direct aggregation of all 500+500 tile summaries and the comparison CSVs:

- Minimal run `wq_quarkonly_1M_20128400`: 1,000,000 rows, 878,707 evaluated, 121,293 skipped = 121,041 `nonperturbative_quark_yukawa` + 252 `quark_fit_failed`. Custodial twin `20675555`: 878,705 evaluated (121,041 + 254).
- Skips cluster by r, not by low M_KK: r=0.05 → 103,211 skips (~52% of draws), r=0.1 → 17,888 (~9%), r≥0.25 → ≤79 per r-slice. Within r=0.05 the skip count mildly INCREASES with M_KK (10,028 @ 1 TeV → 10,631 @ 50 TeV, ~6% relative).
- Survival fractions per (r, M_KK) cell are almost perfectly step functions: across ~9,200 (model, r, M_KK, constraint) cells in `comparison/constraint_veto_by_r_mkk.csv`, only 10 have veto fraction strictly inside (0,1): B003 0.639 @ (r=0.1, 1 TeV), B004 0.571 @ (r=1, 2 TeV), B013 ~0.517 @ (r=0.25, 1 TeV), B001 0.429 @ (r=1, 1 TeV), B011 0.333 @ (r=0.1, 1 TeV) (each mirrored minimal/custodial). Everything else is exactly 0.0 or 1.0.
- T010 (minimal): 694,123 vetoes total, 0 in custodial — matches ledger claim exactly. Its veto fraction is a pure step: 1.0 for M_KK ≤ 20 TeV (r ≥ 0.1) dropping to 0.0 at 30 TeV; for r=0.05 it is 1.0 at 15 TeV and 0.0 already at 20 TeV.
- EW001 vetoes are IDENTICAL in minimal and custodial (352,491 each), as are CR001 (352,491), CR012/CR013 (264,521), B012 (186,476): the custodial branch changes nothing in the proxy set that defines the inclusive floor.
- Exceptions in production: ZERO (`constraint_exceptions=0` in every tile summary, both runs).
- `hard_not_evaluated` in production: T001 and T002 on EVERY evaluated row (878,707 / 878,705 rows each), both runs.

---

## Findings

### F1. BLOCKER (interaction) — Both headline floors are produced, end to end, by constraints other slices have confirmed buggy; the harness aggregates them faithfully, so the headline numbers are not physics-defensible

- Files: `scan_outputs/wq_quarkonly_1M_20128400/analysis/analysis_report.md`; `docs/STATE_OF_PROJECT.md:5`; `.orchestration/PHASE2_PROGRAM_LEDGER.md:81,98`; data verified from `comparison/constraint_veto_by_r_mkk.csv`.
- Issue: The minimal "rigorous floor ~25-30 TeV" is set 100% by T010: it is the only rigorous constraint vetoing anywhere above ~2-3 TeV in the minimal run (694,123 of the rigorous vetoes; T011 only adds 119,859 at ≤2 TeV), and its veto fraction is a Yukawa-independent 1.0/0.0 step. Slice context (a) establishes T010's dominant term is wrong by sign and ~190x, gating at 1σ where the SM's own R_b pull is −0.996. A ~190x-too-large δg_bL inflates the M_KK at which the step turns off by roughly sqrt(190) ≈ 14x (δg ∝ 1/M²). The "custodial strict floor 2-3 TeV" survives only because P_LR zeroes the same buggy term; the headline minimal-vs-custodial contrast (T010-dominated) is therefore a contrast between a buggy constraint and its removal. The custodial INCLUSIVE 7 TeV floor is set by EW001 + CR proxies; EW001 carries the ×6-weak Λ_IR-vs-physical-M_KK coefficient bug (context (c)), so 7 TeV is understated by up to ~×2.45 in scale if ΔT dominates EW001's verdict.
- Why it matters: every harness mechanism I audited below is working as designed — the problem is strictly upstream constraint physics. But the deliverables of THIS slice (analysis reports, comparison artifacts, Scan Explorer JSON, STATE/ledger prose) present these numbers as the headline result with "rigorous" labels.
- Impact on headline floors: minimal 25-30 TeV invalid (would collapse to the ΔF=2/EW001 regime, itself affected by (b) and (c)); custodial strict 2-3 TeV weakened ~×(1.7)^(1/3)–ish by the ε_K bug direction (too weak → true floor higher); custodial inclusive 7 TeV understated via EW001.
- Confidence: high on the dependency chain (verified from raw veto counts); the underlying constraint bugs are taken from the other slices' findings.

### F2. MAJOR — "~25-30 TeV" is more precise than the scan supports, omits the measured r=0.05 exception, and its precision actually comes from a backwards literature conversion, not the grid

- Files: `.orchestration/PHASE2_PROGRAM_LEDGER.md:81,98,106`; `docs/STATE_OF_PROJECT.md:5`; `docs/quark_scan_constraint_update_2026-06.md:79-92`; data: minimal `analysis_report.md` survival table.
- Issue: The grid {1,2,3,5,7,10,15,20,30,50} TeV has no point at 25. Measured: strict survival is exactly 0.0 at 20 TeV and exactly 1.0 at 30 TeV for r ≥ 0.1 → the data supports only "floor in (20, 30] TeV". For r = 0.05 the minimal run has FULL strict survival already at 20 TeV (floor in (15, 20]); no prose I found states this r-dependence of the minimal floor (the ledger gives r-dependence only for the custodial ΔF=2 floor). The "25-30" range numerically equals 10-12 TeV (literature Λ_IR bound) × x₁ ≈ 2.45 (`docs/quark_scan_constraint_update_2026-06.md:85-87`) — i.e., it is the literature bound converted to the physical convention, presented alongside the scan as if the scan resolved 25 vs 21.
- Prose-overclaim check (slice question 5): I found NO instance of "excluded below X" language; STATE_OF_PROJECT:13/113 correctly describes the floor as a discrete display convention and ledger/docs use "~". The defect is the unsupported precision and the missing r=0.05 caveat, not exclusion language.
- Impact on headline floors: presentation-level; correct statement is "minimal strict floor between 20 and 30 TeV for r ≥ 0.1, between 15 and 20 TeV for r = 0.05 (subject to F1)".
- Confidence: high.

### F3. MAJOR — `hard_not_evaluated` conflates "could not evaluate" with "evaluated, passed, but tag=partial"; T001/T002 are flagged as HARD coverage gaps on 100% of rows because of a literal "proxies" vs "proxy" substring miss; a FAILING partial HARD constraint would silently never veto while being mislabeled "not_evaluated"

- Files: `scripts/run_full_catalog_scan.py:908-916` (`_classify_results`: the `else` branch appends every HARD result whose tag is not rigorous/proxy to `hard_not_evaluated`, even when `evaluated=True` and `passes=True`); `scripts/run_full_catalog_scan.py:984-987` (`tag_result`: `"proxy" in needs_text` fails on the word "proxies"); `flavor_catalog_constraints/primary/top_higgs_ew/T001.py:351` (needs_human text "...are used as neutral-current flavor-overlap proxies..."); `scripts/build_wq_quarkonly_comparison.py:546-563` (`_draw_veto_rows` writes these rows with `veto_class="not_evaluated"`, `evaluated=False` — factually wrong for T001/T002, which were evaluated and passed with ratios ~1e-8).
- Verified: in production both T001 and T002 are `evaluated: True, passes: True, tag: 'partial'` on every row, yet appear in `hard_not_evaluated` 878,707 times each, forcing `coverage_complete=False` and an advisory `hard_coverage_gap:T001,T002` on EVERY row of both 1M runs. Neither `analysis_report.md` nor the comparison README discloses this; the analysis report has no coverage section at all.
- Why it matters: (i) the coverage metric is dead weight — it can never be clean, so a real coverage regression would not stand out; (ii) the brittleness mode the slice prompt asked about is live: one English plural flipped two HARD constraints from proxy to partial; had they been failing instead of passing (ratios are 1e-8, so they never would here), they would silently not veto under BOTH strict and inclusive while the paired-vetoes artifact recorded them as not evaluated.
- Impact on headline floors: none numerically today (T001/T002 pass by 7-9 orders of magnitude). Integrity/reporting defect.
- Confidence: high (reproduced from raw rows).

### F4. MAJOR (interaction) — The custodial INCLUSIVE 7 TeV floor is model-independent by construction in this scan: every constraint that defines it (EW001, CR001, CR012, CR013, B012) has bit-identical veto counts in minimal and custodial

- Files: `comparison/constraint_veto_by_r_mkk.csv` (totals above); `.orchestration/PHASE2_PROGRAM_LEDGER.md:81` ("custodial INCLUSIVE floor ~7 TeV (proxy EW001 oblique-S now dominant...)").
- Issue: the ledger phrasing "EW001 ... now dominant" suggests a custodial-specific result; in fact EW001's verdict is unchanged between models (352,491 vetoes each, a pure step at the same grid edge), as are all CR proxies. The custodial run changes ONLY the T010/T011 column. The inclusive floor is the same step-proxy wall in both models; in the minimal model it is merely hidden under T010. Combined with context (c) (EW001 ΔT ×6 weak), the 7 TeV inclusive number inherits that bug identically in both models.
- Per the slice-5 scale-leak sweep instruction, I checked the OTHER proxies for the EW001 pattern (literature coefficient in Λ_IR convention fed a physical M_KK): EW001 is the only allowlisted constraint using a `v²/M_KK²`-coefficient anchor (grep over top_higgs_ew/, collider_rs/, B011-B014: only `EW001.py` matches). CR001/CR012 compare the physical KK mass against experimental resonance-mass limits (correct convention, `CR001.py:340`, `CR012.py:22`); CR007/CR013 map gauge mass → Λ_IR via `GAUGE_KK_ROOT_NN` and then m_G = x_G1·Λ_IR (`physics_adapters/kk_graviton_resonance.py:78-116` — correct, double-divides nothing); B011-B014 use `kk_ew_mass_gev` (physical) only as the RG matching scale (log sensitivity only). No second EW001-style leak found.
- Impact on headline floors: the "custodial inclusive 7 TeV" should not be narrated as a custodial-physics outcome; it is the shared proxy wall, and it moves with the EW001 fix.
- Confidence: high on the identical-counts fact and on the absence of further v²/M² leaks in the quark allowlist; medium on how far 7 TeV moves (depends on whether ΔS or ΔT controls EW001's edge).

### F5. MINOR — ~12% of draws (52% at r=0.05) are dropped pre-evaluation; all quoted fractions are conditional on the fit+perturbativity filter; the analysis report does not flag the halved/conditioned r=0.05 column

- Files: `scripts/run_full_catalog_scan.py:835-851` (`_require_valid_quark_fit`, `_require_perturbative_leptons` analog `|Y| ≤ 4`), `:1793-1803` (`_skip_reason`); data above.
- The denominator question (slice check 2): every published fraction consistently uses EVALUATED rows — `analyze_wq_quarkonly.py:303,623-624` (survives/evaluated, vetoed/evaluated), `build_wq_quarkonly_comparison.py:591-601,720` (skipped excluded; `veto_fraction = vetoed/points` where points only counts rows carrying a constraints dict, i.e. non-skipped), `build_scan_explorer.py:287-291,404`. No denominator mixes in skipped rows; skipped rows carry `survives_*=False` in raw rows but are never counted as survivors or as denominator.
- Drops do NOT cluster at low M_KK — they cluster at low r and rise slightly with M_KK (data above), so they cannot fake an M_KK floor. But at r=0.05 the anarchic measure is conditioned (only draws whose fitted |Y| ≤ 4 survive), which both halves the statistics and tilts the surviving Yukawa distribution; given the step-function behavior of the binding constraints this has no visible effect on floors, but the r=0.05 column should carry the caveat in the report. Ledger:125 already calls the skip "honest"; the per-r concentration is undisclosed in the analysis report.
- Confidence: high.

### F6. MINOR — Minimal/custodial pairing is real but not bit-perfect: 2 of 1,000,000 draws flipped fit success between runs (252 vs 254 `quark_fit_failed`; 878,707 vs 878,705 evaluated) despite byte-identical inputs

- Files: tile summaries (aggregated above); visible in `comparison/survival_by_r_mkk.csv` rows (0.25, 15.0) and (0.25, 30.0) as minimal_evaluated 19996/19990 vs custodial 19995/19989.
- Verified the draw path (slice check 4): in quark-only mode the per-draw RNG (`default_rng(tile.seed + draw_idx)`, `run_full_catalog_scan.py:317-318`) feeds ONLY `_draw_quark_seed` (`:440`); `ew_model` never touches the RNG stream or `fit_quark_sector`'s inputs (`:448-456`). Empirically compared 150 paired rows across three tiles: 0 Yukawa-seed mismatches and identical fitted `bulk_c_Q`. The sbatch launchers differ only by `--ew-model` (same `--base-seed/--tile-seed-stride/--quark-fit-r`). The 2-row delta is therefore environment-level FP nondeterminism in the SciPy least-squares fit (different nodes/BLAS), flipping marginal convergence. `_validate_pairing` still passes because keys pair 1:1; only the evaluated/skipped status differs on those 2 rows.
- Impact: negligible (2 ppm); but "paired draw-for-draw" claims should say "identical inputs, 2 ppm convergence-status exceptions".
- Confidence: high on mechanism location (fit), medium that BLAS/node variation is the specific cause.

### F7. MINOR — Scan Explorer builder/page nits

- `flavor_catalog/website/src/pages/explore.astro:162`: "first grid point where the median veto fraction is <= 0.5" — the quantity is a per-cell MEAN fraction of vetoed draws (`build_scan_explorer.py:404`), not a median. Wording bug.
- `build_scan_explorer.py:299-309`: the explorer counts a veto for ANY active+evaluated HARD failure regardless of tag — unlike the harness/comparison, which require tag ∈ {rigorous, proxy}. Today this is invisible (the only partial HARD constraints, T001/T002, always pass — cross-checks agree with the CSV to ~1e-7, `.orchestration/runs/WEB-EXPLORER/data_impl_summary.md`), but a future failing partial HARD constraint would appear in the website's veto fractions while being excluded from strict AND inclusive floors. Latent inconsistency between the public display and the published semantics.
- `build_scan_explorer.py:550-588` `_cross_checks`: computes `abs_diff` against the comparison CSV but never asserts a tolerance — a regression would be printed, not fail the build.
- Envelope floor = max over per-constraint floors (`explore.astro:243-253`), not the floor of joint survival; with this dataset's step behavior they coincide, and STATE_OF_PROJECT:113/162 documents the convention. Acceptable as display convention; would diverge if several constraints had correlated ~0.4 fractions.
- Threshold sensitivity (slice check 5): empirically a non-issue for THIS dataset — only 10 of ~9,200 cells have fractional veto fractions, none of which sit between 0.05 and 0.5 at a floor-defining transition. Moving FLOOR_THRESHOLD from 0.5 to 0.05 changes no floor; it would only pull B001 (0.43 @ r=1, 1 TeV) and B011 (0.33 @ r=0.1, 1 TeV) into the "biting" display set at cells already excluded by step constraints.

### F8. MINOR — `tag_result` string matching remains brittle beyond the "proxies" case

- `scripts/run_full_catalog_scan.py:990-999`: a `matching_status` containing the word "rigorous" anywhere (e.g. "proxy pending rigorous treatment") escapes the proxy tag via `"rigorous" not in status_text` and, if it also lacks "recast", lands on `rigorous` at `:1002`. The resolved-phrase allowlist (`:991-997`) covers only five exact negations. The tests (`tests/test_full_catalog_scan_harness.py:94,115,139`) lock in current phrasings but cannot protect against new constraint wording. Recommend a structured boolean diagnostic (e.g. `tag_hint`) instead of prose matching.
- Also note `_proxy_flags` (`:1058-1070`) marks any diag key containing "proxy"/"recast" with a truthy value — CR001's float `sigma_times_br_proxy_pb` correctly trips it, but a hypothetical `no_proxy_used: True` would misclassify a rigorous constraint as proxy (conservative direction).

### F9. MINOR — `run_universal_c_sanity` is skipped in quark-only mode, so the production lanes ran with no decoupling sanity gate

- `scripts/run_full_catalog_scan.py:1750-1756`: quark-only mode replaces the sanity with `{"skipped": true}`. Mitigation observed in data: at 30-50 TeV every veto fraction is exactly 0 in both runs, which is an empirical decoupling check; and `_is_sm_tension_only` (`:1025`) — the guard built to catch SM-tension-gated constraints like T010 (R_b pull −0.996) — is only ever invoked from the sanity path, so it provided zero protection in production. Given F1, this is precisely the gate that might have raised a flag earlier (a 50 TeV universal-c minimal point passes T010, so it would NOT have caught the 190x bug at the default sanity scale — the protection gap is structural, not counterfactual).

---

## Verified correct

1. **Scale-convention chain (slice check 1)** — single conversion, used once, each family on its intended convention:
   - M_KK → Λ_IR happens exactly once per tile: `_build_tiles`, `run_full_catalog_scan.py:1569` (`lambda_ir = mkk / cfg.xi_kk`).
   - The RS-EW spectrum is built from Λ_IR and returns the PHYSICAL first gauge KK mass: numerically `kk_ew_mass_gev = 2998.65` for M_KK=3000 input (ratio 0.99955) — the 0.05% offset is DEFAULT_XI_KK=2.4487 vs the actual Bessel root ≈2.4476 at this ε; cosmetic.
   - ΔF=2 / KK-gluon: `compute_quark_kk_gluon_couplings` is called with explicit `M_KK=tile.mkk_gev, xi_KK=cfg.xi_kk` (`:487-492, 597-602`) — the dangerous `xi_KK=1.0` default in `quarkConstraints/couplings.py:100` (the README's LFV ξ=1 warning) is never exercised in the scan path; `kk_gluon_mass_gev = tile.mkk_gev` (physical) for the 1/M² propagator. Correct convention (the ΔF=2 problems are the slice-(b) matrix-element issues, not scale leaks).
   - Zbb/T-family: consume `rs_ew_couplings` built from the actual spectrum (first-principles, no literature-coefficient convention to leak; T010's defect is the separate slice-(a) bug).
   - Collider: CR001/CR012 compare physical mass to limit; CR007/CR013 graviton mapping divides by the gauge root before multiplying by the graviton root (`kk_graviton_resonance.py:78-116`), with couplings-carried ξ=2.4487 taking precedence. Rows confirm `params: M_KK/Lambda_IR = 2.4487`, `xi_KK = 2.4487`.
   - Dipoles B011-B014: physical `kk_ew_mass_gev` as matching scale (log-only sensitivity).
   - Conclusion: ONE convention defect total in the scan's constraint surface (EW001 coefficient, already slice-(c)); no harness-side leak; no quark constraint inherited ξ=1.
2. **878,707 / 1,000,000 accounting (check 2)** — fully explained (121,041 perturbativity + 252 fit failures), reproducible from summaries, consistent across analysis report, comparison CSV, and explorer; no low-M_KK clustering; denominators uniformly exclude skips.
3. **Strict vs inclusive semantics (check 3)** — code (`:932-933`) matches docs (STATE:5): strict = no rigorous HARD veto; inclusive = no rigorous AND no proxy HARD veto (inclusive ⊆ strict, correctly described as the smaller set). `partial`/`stub` can never veto (`:908-916`, `_accumulate_constraint_tally:1140-1145`, comparison `:539-545`, all require tag ∈ {rigorous, proxy}). Exceptions ARE counted and surfaced at every level: per-row `exception_type` + advisory flag (`:902-926`), tile `constraint_exceptions`/`exception_ids` (`:362-370`), run-summary `exception_ids` top-20 (`:1421`); production had zero. The isolation wrapper `_evaluate_constraint_ids` (`:942-965`) mirrors `registry.evaluate_all` and converts exceptions to failing stub results (which land in hard_not_evaluated, not silent passes).
4. **Seeding & disjointness (checks 4/6)** — plan arithmetic correct and conservative: max per-shard span = 9*1,000,003 + 2,000 = 9,002,027 < SHARD_SEED_BLOCK 20,000,000 (the plan's printed requirement 10,002,030 is stricter than needed and still passes); 50 tasks × 10 tiles × 2,000 draws = 1M; per-draw `default_rng(seed)` with sequential integers is safe under PCG64 SeedSequence. Minimal and custodial intentionally share seeds for pairing; comparison builder enforces 1:1 pairing with duplicate detection (PRIMARY KEY + explicit checks) and scan-plan equality minus run-identity fields, which keeps `base_seed` in the compared surface.
5. **Yukawa draw measure (check 6)** — uniform Re, Im ∈ [−1.5, 1.5] per element with |y_ij| ≥ 0.1 enforced by rejection (≤256 tries, then phase-preserving magnitude clamp — probability of reaching the clamp is negligible); identical config in both runs (config-hash surfaces differ only by `ew_model`); perturbativity cut |Y|_max ≤ 4 applied identically (same skip counts to within the 2 ppm of F6).
6. **Numbers spot-check (check 7)** — analysis-report row/evaluated/survival values match `survival_by_r_mkk.csv` cell-for-cell on all sampled cells; ledger's "T010 694123 minimal / 0 custodial" reproduced exactly; explorer cross-checks vs CSV agree to ≤5e-7; survival columns in the reports are fractions (%.6g), so "1"/"0" = 1.0/0.0 exactly.
7. **Resume/atomicity** — tiles write tmp files + `os.replace`; `_completed_summary` requires `complete` flag, matching config hash, and matching draw counts before skipping (`:1485-1500`); summaries carry git SHA + dirty flag.
8. **Tests** — `tests/test_full_catalog_scan_harness.py`: 18/18 pass in 7s in the current tree; cover tag policy phrasings, seed/hash determinism, resume gating, allowlist vs deferred-lepton split, quark-only draw determinism, custodial model threading, tally semantics, and ΔF=2 ratio relaxation with M_KK. Gap: no test pins the `hard_not_evaluated` semantics for an evaluated-but-partial HARD constraint (F3), and no test asserts explorer/comparison veto-definition equivalence (F7).

---

## Severity summary

| Severity | Count | IDs |
|---|---|---|
| BLOCKER | 1 | F1 (headline floors ride entirely on constraints confirmed buggy by other slices; harness itself faithful) |
| MAJOR | 3 | F2 (floor precision/derivation + r=0.05 omission), F3 (hard_not_evaluated conflation + "proxies" tag miss on T001/T002, undisclosed 100% coverage-gap), F4 (inclusive 7 TeV floor identical across models; EW001 unchanged by custodial branch) |
| MINOR | 6 | F5 (r-clustered 12% skip conditioning), F6 (2 ppm pairing nondeterminism), F7 (explorer wording/definition/assert nits), F8 (tag string-matching brittleness), F9 (no sanity gate in quark-only lanes), xi 2.4487-vs-2.4476 cosmetic offset |

Overall: the harness, seeding, pairing, classification plumbing, and bookkeeping are solid and the published numbers are faithfully computed from the raw rows — I found no harness-side scale leak and no denominator games. The defects that matter are (i) what flows INTO the harness (T010/EW001/ΔF=2 bugs make both headline floors non-defensible as physics), and (ii) how the outputs are NARRATED (unsupported 25-30 precision, undisclosed permanent T001/T002 coverage gap, custodial-inclusive floor framed as custodial physics when it is model-independent in this scan).
