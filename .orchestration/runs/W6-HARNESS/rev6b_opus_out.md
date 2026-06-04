APPROVE

Findings (concise)

Structure — matches the plan
- Driver `scripts/run_full_catalog_scan.py` is a new harness outside `flavor_catalog_constraints`. Per-tile `(Lambda_IR=mkk/xi_kk, M_KK)` is fixed in `_build_tiles` (1157–1172). `RSEWSpectrum.build` + `RSEWOverlapSplineCache.build(...include_omega=True...)` are built once per tile in `_run_tile` (187–204) and injected via `build_from_rs_ew_inputs(..., spectrum=spectrum, rs_ew_cache=overlap_cache, include_charged_current=True, include_fermion_kk_mixing=True, include_higgs_yukawas=True, lepton_yukawa_result=...)` (371–387). The W6a builder enforces `rs_ew_cache.spectrum is spectrum` and rejects callable+cache mixing (`rs_ew_builder.py` 215–373), so injection is validated end-to-end.
- Quark fit uses a real `QuarkFitResult` from `fit_quark_sector` against `default_quark_targets()` and a deterministic `QuarkFitSeed` derived from SVDs of anarchic complex 3×3 (460–555); lepton inputs draw c_L, c_E[3], c_N (uniform), log10 M_N, ordering, Majorana phases (517–534), then `compute_all_yukawas` (356–368). c=0.5 singularity guard runs on both quark-fit c-values and lepton c-values (344–354, 623–626).

Determinism & registry
- `worker_init` asserts `registry.discover()` → exactly `cfg.expected_registry_count=103` and `import_failures()` empty, else raises (132–146).
- Tile seeds = `base_seed + stride*tile_id`, draw seed = `tile.seed + draw_idx`, per-draw RNG via `np.random.default_rng(draw_seed)` (216–218, 1163–1170); `tests/test_full_catalog_scan_harness.py::test_default_tile_seed_stride_covers_large_smoke_tiles` proves stride > n_draws (no seed overlap), and the explicit-stride test pins seeds `[11, 28]`. `_config_hash` is over the full `ScanConfig`.

Atomic resume
- Per-tile JSONL+summary written to `*.tmp.<pid>` then `os.replace` to final (215–297). `_completed_summary` requires `complete=True`, matching `config_hash`, and `n_requested == n_rows == expected_draws` (1091–1106), exercised by `test_completed_summary_resume_requires_complete_matching_hash_and_draws`.

Honest tagging / veto policy
- `_classify_results` (629–691): only `Severity.HARD` can populate exclusion lists; rigorous-failing→`excluded_by_rigorous`; proxy-failing→`excluded_by_proxy`; HARD partial/stub/exception/missing→`hard_not_evaluated`; SOFT/INFO non-pass→`advisory_flags`. `tag_result` (694–728) routes `missing_extra`/`exception_type`/`evaluated=False` to `stub`, proxy flags or `needs_human` containing "proxy"/"recast" to `proxy`, generic `needs_human`/`matching_status` partial/needs-human/deferred to `partial`, with a "proxy resolved" allow-list. Severity invariant covered by `test_hard_veto_policy_separates_rigorous_proxy_partial_stub_and_advisory`; resolved-proxy wording test pins `B022`-style status as rigorous. `survives_all_HARD_strict` requires no rigorous exclusion; inclusive additionally requires no proxy.

Universal-c sanity
- `run_universal_c_sanity` (834–967) builds its own spectrum/cache at `sanity_mkk_gev/xi_kk`, replaces `bulk_state` with universal c, sets `U_*=I`, `ckm=I`, force-degenerates `Y_N` before re-applying PMNS, and filters SM-tension-only exclusions via `_is_sm_tension_only`. Smoke confirms `passes_no_spurious_hard_exclusions=True`, `rigorous=[]`, `proxy=[]`.

Smoke metrics (smoke_w6b_fresh)
- 4 tiles × 2,500 = 10,000 rows; 533 evaluated; skip dominated by `nonperturbative_lepton_yukawa` (~94.7% — sampling artifact of the wide `c`/M_N priors, not a harness defect); spline `max_a_rel_err ≈ 4.76e-8` (≪ rel_tol 1e-3); 103 constraints, 98 evaluated per point, exception rate 0.0485 (E001/E002/L006/L010/L023 each 138/tile = systematic per-point); tag distribution per tile rigorous 2,760 / proxy 6,624 / partial 4,140 / stub 690. Top rigorous vetoes B022/K004/L001/B003/B004/B023/K001; top proxy CR009/CR006/B016/CR001/CR005/CR012/CR013/EW001/K010/B015.

Tests
- `pytest tests/ -q → 1706 passed, 1 skipped`. Harness unit tests cover policy, resume, tile seeds, and the "proxy resolved" carve-out.

Observations (non-blocking, for the perf gate decision — not code defects)
- `post_cache_seconds_per_evaluated_point ≈ 7.24s` (vs plan §18 target 1.2–3 s/point) extrapolates to ~2×10⁵ core-h for 1e8 draws. Drivers of cost on the evaluated path are full 103-constraint evaluation + ΔF=2/semileptonic adapters; this is the gate the smoke is supposed to surface, and the harness reports it honestly.
- High skip rate driven by lepton-perturbativity rejections; tightening priors (or moving the perturbativity guard upstream of the quark fit) is a follow-up tuning task, not a correctness fix.
- `cache_metrics.rs_ew_cache_injected` in the skipped-row branch is hard-coded True (449); the skip path can short-circuit before injection. Cosmetic — does not affect physics or exclusion lists.
