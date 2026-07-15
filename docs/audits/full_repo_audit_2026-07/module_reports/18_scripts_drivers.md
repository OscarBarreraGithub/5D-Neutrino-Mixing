# Report 18 — Gap-closing: comparison drivers, benchmark scripts, exports/plots, website builders, tools

**Structural summary.** The never-read drivers split into (a) post-hoc rescaling/comparison scripts reinterpreting stored scan ratios, (b) self-referential benchmarks, (c) aggregation/website builders. The rescaling family is where the real problems live: the 2007-vs-modern comparison rescales modern hadronic-budget ratios by legacy operator-weight bound ratios (different objects, 6 orders of magnitude apart in units), and both "validation" benchmarks in scope validate the code against itself. The anarchic forward stack (run_rs_anarchy.py, anarchic_complex_m12.py) is the cleanest physics code reviewed, with a genuine selftest. The website explorer's envelope floor is statistically optimistic (max of marginals ignores joint vetoes).

### [CRITICAL] 2007-vs-modern rescale factors for B_d/B_s/D0 are not the bounds the scan used
- **File:** scripts/compare_2007_vs_modern.py:35-41
- **Category:** convention-inconsistency | units
- **Claim:** `ratio_2007 = ratio_modern × (bound_modern/bound_2007)` uses legacy operator-weight bounds ('B_d': 4.0e-7/6.0e-7, 'B_s': 5.5e-6/1.5e-5, 'D0': 8.5e-9/8.0e-8 — matching LEGACY_OPERATOR_WEIGHT_BOUNDS gated behind use_hadronic=False) while the scan's ratio_to_bound_by_system used hadronic |M12| budgets in GeV (1.667e-13, 5.844e-12, 3.125e-15 — verified numerically). The claimed 9.4× D0 tightening (0.106 factor) comes from a normalization the scan never used; "tightening driver" and acceptance-delta figures inherit it. Only ε_K (4.18e-4 — a third value alongside the 6.7e-5/3.04e-4 split-brain) and K numerators match. Suspicious: ε_K and K rescale factors both exactly 0.697 (one denominator likely back-derived).
- **Confidence:** high (numerator mismatch), medium (net figure impact)
- **Fix:** Rescale by same-normalization bound pairs; document 2007-value provenance.

### [MAJOR] 2007-vs-modern acceptance comparison asymmetric (fit/verifier gates only on modern side)
- **File:** compare_2007_vs_modern.py:72,86
- **Claim:** accepted_modern includes fit success + verifier + ratio gates; accepted_2007 is ratios-only — fit-failed rows count as "accepted in 2007, rejected by modern", inflating claimed tightening.
- **Confidence:** high (asymmetry), medium (magnitude)
- **Fix:** Same gates both sides.

### [MAJOR] Lane C "benchmark" PASS is a structural self-check, not physics validation
- **File:** scripts/benchmark_quark_0710_1869.py (7086 lines; checks delegated at :886-897)
- **Category:** circular-validation
- **Claim:** Checks import isolation, schema/convention IDs, artifact hashes, and runs the package's own verifier (expectations are strings like "kaon.q1_vll_vrr.8over3_fk2_mk2_bk_mu.v1") — a PASS cannot detect the inverted RG, missing volume factor, or ×4 ME. No paper target numbers exist anywhere in the 283 KB script. export_quark_0710_1869_artifacts.py exports these self-certified artifacts as "canonical PR6 paper artifacts".
- **Confidence:** high
- **Fix:** At least one hard numeric acceptance check against independent FPR/BMU references.

### [MAJOR] audit_wilson_rg shares every physics input with the evolver it audits
- **File:** scripts/audit_wilson_rg.py:26,79,100-117
- **Category:** circular-validation
- **Claim:** "Independent LO textbook calculation" imports run_alpha_s from the code under test and hardcodes the identical ADM [[-16,-6],[0,2]] — validates only propagator linear algebra; no pass/fail gate. Mitigating: the agent independently re-derived the BMU LR ADM conjugation (T γ T⁻¹ with T=[[0,1],[-2,0]]) — the shared constant is correct.
- **Confidence:** high
- **Fix:** Compare against published Buras/UTfit NLO magic numbers with a tolerance band, and assert.

### [MAJOR] Perez-Randall audit conclusion "O(3-4)×, not a convention drift" overstates
- **File:** scripts/audit_perez_randall_consistency.py:44,95-97
- **Category:** factor-of-2 | statistics
- **Claim:** The 3.5× kY_N tension is largely absorbable by convention freedom the script does not scan: seesaw prefactor 2-vs-4, v=174 vs 246, and the hardwired M_N=M̄_Pl/10. Numerically: repo (2, v=174) → 3.52× tension; (4, v=174) → 2.49×; (4, v=246) → 1.76×; ×3 uncertainty in M_N closes it entirely. The script's arithmetic is faithful to the repo's Eq. (6) reduction — the bug is in the inference. Also resolves benchmark_perez_randall.py's factor 5-7: same convention stack + different geometry, not a coding error.
- **Confidence:** medium (paper text not accessible)
- **Fix:** Print the tension under all four (prefactor, v) conventions and vs M_N.

### [MAJOR] Website explorer envelope floor ignores joint vetoes (max of marginals ≠ union)
- **File:** flavor_catalog/website/src/pages/explore.astro:231-253; build_scan_explorer.py:47,351-413
- **Category:** statistics
- **Claim:** Per-constraint floor = first M_KK with veto ≤ 0.5, envelope = max over constraints; three constraints each vetoing 40% of disjoint draws → no individual crosses 0.5 yet ≥50% jointly vetoed — envelope understates the combined floor. Constraints never exceeding 0.5 anywhere are dropped entirely (biting_ids). Secondary: veto-fraction denominator is all evaluated rows, not per-constraint-evaluated rows; "median veto fraction" caption describes a mean.
- **Confidence:** high (mechanism), medium (size — one dominant constraint per cell masks it)
- **Fix:** Stream a per-cell "vetoed by ANY" fraction and drive the envelope from it.

### [MINOR] "Exact" ξ_KK rescale (ratios ÷ ξ²) ignores RG start-scale and α_s(M) drift
- **Files:** export_accepted_quark_scan.py:10-11,59-90; plot_publication_figures.py:181-236
- **Claim:** Exact only for 1/M²; Wilsons were matched/run from Λ_IR not m_gkk: C4 evolution 3→7.35 TeV adds ≈1.09 dropped by the rescale (~9% on ε_K/K ratios, partially offset by g_s² ≈ −7%) → few-% bias on published floors.
- **Fix:** Document or re-evolve.

### [MINOR] CFW comparison uses fixed g_s=1.05 while the run's g_s(M_KK) spans 1.06→0.955 per tile; M_KK_min = M_KK·√max_ratio freezes g_s and RG along the extrapolation — up to ~10% per-tile systematic on the CFW-marker comparisons (10.5/21/17/33 TeV).
- **File:** rs_anarchy_cfw_comparison.py:41,301,453
- **Fix:** Rescale per draw by gs_star/g_s(M_KK_tile).

### [MINOR] Inconsistent silent fit-quality filters and NaN handling across the figure pipeline
- **Files:** analyze_quark_scan.py:100,115,419; plot_publication_figures.py:154-159; plot_allowed_region.py:120
- **Claim:** Three scripts consuming the same JSONL disagree on fit_score filtering (drop >0.1 + NaN vs default-to-0.0-keep vs no filter); analyze_quark_scan maps non-finite system ratios to 0 (missing system silently passes in max-ratio figures); both fall back M_KK → Lambda_IR (latent ×2.45 axis conflation); exclusion contours drawn on gaussian_filter(σ=1.5)-smoothed grids can shift the ratio=1 boundary.
- **Fix:** Centralize row loading; NaN ratios = excluded; assert M_KK present.

### [NOTE] benchmark_perez_randall is a repo-regression test, honestly labeled (expected values repo-generated; PASS = unchanged since snapshot).
### [NOTE] Legacy-g_s benchmarks/replays (benchmark_quark_mfv, gen_sidebyside_points, extract_plot_quantities, validation.py gates) are internally consistent with Lane B (g_s≈1) but must never be quoted next to modern g_s*=3 numbers (×9.4). anarchic_complex_m12 --selftest verifies direct M12 = public path to 1e-9 (good); its docstring promises a ratio_dm_K column never emitted (cosmetic).
### [NOTE] calibrate_phase0 tolerance floors calibrated to the model's own achievable residuals (p95), not PDG uncertainty; pdg_2sigma recorded but not enforced.
### [NOTE] meg2_constraint_comparison: "100×" should be 80×; C_PAPER=0.02 is 15% inconsistent with its own BR-scaling rule (0.0173); MEG II labeled 2024 vs 2025 elsewhere.

## Verified correct
anarchic_complex_m12 selftest + complex-M12 emission (κ_ε=0.94, correct Blanke/Gedalia-plane inputs); anarchic_reproduction_extract amplitude recovery (budgets numerically equal the run's denominators; ΔM_K/2 confirmed); run_rs_anarchy core (CKM ordering, ascending SVD, mass-basis f²-overlaps, per-tile deterministic default_rng — no hash()); BMU LR ADM conjugation independently re-derived; compute_per_quark_residuals/reclassify_scan (delegation); sharded sbatch seeding (base_seed + shard_id); export_collaborator_5tev_points convention hygiene; build_scan_explorer hard-partial veto exclusion (rigorous/proxy only); tools/aggregate_factchecks prefix→family mapping.

## Not reviewed
analyze_wq_quarkonly.py, build_wq_quarkonly_comparison.py, wq_quarkonly_1m_plan.py, c02c_post_process.py, anarchic_bauer_s1_zbb.py, yukawa_envelope/anatomy scripts, rs_anarchy_mkk_min_hist* variants, plot_rs_anarchy_summary.py, plot_quark_deltaf2_r_sweep.py, compare_quark_bulk_mass_maps.py, export_accepted_quark_scan_with_yukawas.py, export_collaborator_direct_affine_points.py, _audit_followup_crossings.py, rs_anarchy_gate_sensitivity.py, notebooks, citation-anchor scripts, most run_*.sbatch wrappers, ingest_catalog.py rendering paths.
