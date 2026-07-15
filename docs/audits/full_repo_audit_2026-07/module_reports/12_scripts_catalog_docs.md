# Report 12 — Scripts, catalog assembly, docs consistency (scripts/, flavor_catalog/, docs floor claims)

## Pipeline structure
1. scripts/run_full_catalog_scan.py (lane B) draws anarchic Yukawa seeds → fit_quark_sector (r=0.25, BulkMassMap) → builds RS-EW point → evaluates the 103-constraint registry (or 46-ID quark-only allowlist) → per-tile JSONL + summary.
2. scripts/build_constraint_matrix.py flattens JSONL tiles into constraint_matrix.parquet.
3. docs/FLOOR_SUMMARY.md lane-B floor cites scan_outputs/fix100k_minimal_20260622T080053/constraint_matrix.parquet (K001).
4. scripts/anarchic_bauer_s1.py (lane A) draws Bauer-prior Yukawas, FN-inverts f_IR for the c's, reuses run_rs_anarchy.py forward physics → per-scenario parquet feeding the ~10/~30 TeV anarchic floor quotes.
5. run_rs_anarchy.py is the fixed-c ACPS ensemble behind rs_anarchy_* plotting.
6. flavor_catalog_constraints/primary (95) + secondary (8) is the physical registry.
7. scan_outputs/wq_quarkonly_1M_*/analysis holds a legacy 1M-draw analysis with its own headline figures.
8. reclassify_scan.py recomputes pass/fail (any-fail aggregation).
9. Doc chain: FLOOR_SUMMARY ← MODEL_CONVENTIONS ← collaborator CONTENT.md.
10. No parquet scan artifact is actually present in the tree.

### [MAJOR] Bauer Y_max=3 prior is off by a factor 2 in the repo convention (effective Y_max ≈ 1.5)
- **File:** scripts/anarchic_bauer_s1.py:24-33, 191-201; run_rs_anarchy.py:88 (DEFAULT_V_GEV=174)
- **Category:** factor-of-2 / convention-inconsistency
- **Claim:** The convention-bridge docstring cancels only ONE of the TWO √2's between repo f_IR and Bauer's F, so the FN inversion localizes fermions a factor 2 too strongly (in F-pair product) for the same drawn |Y| ≤ 3 — effectively Bauer's prior with Y_max = 1.5.
- **Evidence:** Verified numerically: f_IR(c_repo, ε) = F_Bauer(−c_repo)/√2 exactly (ratio 1.4142 at several c). Bauer m = (246/√2)·F_L·Y·F_R = 174·F_L·Y·F_R; repo m = 174·f_L·Y·f_R = 87·F_L·Y·F_R — ratio exactly 2.0 for identical (Y, c). Since ε_K's dominant C4 ∝ (f_Q1 f_d1)(f_Q2 f_d2) ∝ 1/Y², this inflates ε_K^NP ~4× and the anarchic M_KK floor ~2× — biasing the "~30 TeV current median / ~10 TeV paper-era" lane-A numbers high.
- **Confidence:** high on the factor-2 math; medium on net floor impact.
- **Fix:** Use m = 2·v·f_L·Y·f_R (repo Ȳ = 2·Y_Bauer) in the FN inversion, or draw |Y| ∈ [0.2, 6] repo-side; re-derive lane-A floors.

### [MAJOR] Evaluated HARD constraints tagged "partial" never veto and are mislabeled "hard_not_evaluated"
- **File:** scripts/run_full_catalog_scan.py:913-921 (_classify_results)
- **Category:** logic-bug
- **Claim:** A HARD constraint that WAS evaluated but tagged "partial" falls through to hard_not_evaluated; if it fails, it excludes in neither survives_all_HARD_strict nor survives_all_HARD_inclusive. The "inclusive" floor is biased low wherever partial-tagged HARD constraints bind; passing partials pollute coverage counts.
- **Confidence:** medium-high
- **Fix:** Explicit hard_partial bucket; decide whether partial failures veto inclusively.

### [MAJOR] Un-bannered stale artifact still asserts the retracted pre-B1 "Z→bb 25–30 TeV floor"
- **File:** scan_outputs/wq_quarkonly_1M_20128400/analysis/headline/_make_headline_figs.py:117,158,176; analysis_report.md
- **Category:** stale-data / doc-code-mismatch
- **Claim:** Shipped headline figures/report state "sharp Z→bb floor near 25–30 TeV" and "Z→bb dominates by ~4×", which FLOOR_SUMMARY explicitly retracts as the B1 bug. No SUPERSEDED banner.
- **Confidence:** high
- **Fix:** Banner or delete.

### [MAJOR] Headline floor artifacts are absent from the tree — lane-A/B floor quotes unreproducible as cited
- **File:** docs/FLOOR_SUMMARY.md:22-24,87-89 vs scan_outputs/
- **Category:** stale-data / pipeline-bug
- **Claim:** The cited lane-B parquet does not exist; zero .parquet / constraint_matrix* files anywhere in the repo; no script computes the documented lane-A statistics ("95% quantile floor", "current median ~30 TeV") — anarchic_reproduction_extract.py has no percentile/median/floor logic. The "sharp wall (100% veto ≤5 TeV, 0% ≥7 TeV)" cannot be checked or regenerated from the repo.
- **Confidence:** high (absence verified; numbers may live on the cluster)
- **Fix:** Commit the parquet (or downsampled) plus the exact floor-extraction snippet.

### [MAJOR] Non-reproducible RNG seeds in the Bauer ensemble (hash(scenario))
- **File:** scripts/anarchic_bauer_s1.py:408
- **Category:** pipeline-bug / statistics
- **Claim:** `seed = base_seed + 1009*idx + 100003*(hash(scenario) % 7919)` uses Python's salted string hash (changes per interpreter invocation) — "same seed" rerun produces a different ensemble; published parquet cannot be regenerated.
- **Confidence:** high
- **Fix:** Deterministic map, e.g. sorted(SCENARIOS).index(scenario).

### [MINOR] flavor_catalog/catalog_index.yaml is an empty stub (processes: [])
- File counts verified 95+8=103, no duplicates/overlap; the yaml index has zero entries — anything consuming it silently renders nothing.
- **Fix:** Populate or delete.

### [MINOR] FN CKM hierarchy drops the Wolfenstein A in f_Q1
- **File:** scripts/anarchic_bauer_s1.py:194
- **Claim:** `fQ = [fQ3·λ³, fQ3·Aλ², fQ3]` uses λ³ instead of Aλ³; with A = 0.826 overshoots f_Q1 by 21% → ε_K^NP ∝ f_Q1² ~+46% at fixed rest.
- **Confidence:** medium-high (O(1) minors partly absorb; PDG gate hides)
- **Fix:** A·λ³.

### [MINOR] S3 "little RS" uses L = 7.0, labeled ln(1e3) (= 6.9078) — 9.7% shift in ε. Fix: math.log(1e3).

### [MINOR] build_constraint_matrix can double-count draws via stale .tmp tiles
- **File:** scripts/build_constraint_matrix.py:83-93 — globs both tile-*.jsonl and tile-*.jsonl.tmp.*; a crashed shard rerun leaves both, duplicating draws in floor statistics.
- **Confidence:** medium (needs crash+rerun)
- **Fix:** Skip .tmp.* when final tile exists.

### [MINOR] Derived collider cut (5.5 TeV) in build_constraint_matrix.py:72 disagrees with FLOOR_SUMMARY "~4 TeV" collider floor — undocumented knob.

### [NOTE] EW003 sits in the quark-only allowlist but its required extra (rs_charged_current) is forbidden — evaluated as stub every row, permanently flips coverage_complete=False. Fix: move to deferred list.

### [NOTE] fit_success defaults to True for rows with missing diagnostics (build_constraint_matrix.py:141). Fix: default False/NA.

## Cross-module inconsistencies
- ε_K budget split-brain: core 6.7e-5 vs catalog K001 3.04e-4 (×4.5) — same-K001 floor differs by entry point (consistent with _EPS_BUDGET at anarchic_bauer_s1.py:97).
- EW001 floor: shipped anchors solve to ~15.96 TeV, not documented 18–20 (FLOOR_SUMMARY:25, CLAUDE.md).
- KK-gluon coupling normalization: perturbative g_s(M_KK) with repo f², no √(2L) IR enhancement; the anarchic "validation" and lane-A faithfulness rest on asserted (not derived) cancellations — the one derivable piece (mass bridge) is off by 2.
- Naming hazard: survives_all_HARD_strict is the *looser* criterion (rigorous-only), inclusive the stricter — inverted names invite misquotes.

## Not reviewed
benchmark_perez_randall.py, benchmark_quark_0710_1869.py, benchmark_quark_mfv.py, audit_perez_randall_consistency.py, audit_wilson_rg.py, calibrate_phase0.py, meg2_constraint_comparison.py, anarchic_complex_m12.py, compute_per_quark_residuals.py, analyze_quark_scan.py internals, export_*/plot_* scripts, sbatch shard-merge, flavor_catalog python internals + website/scripts/ingest_catalog.py, tools/, scanParams interplay, run_full_catalog_scan.py:1178-1933.
