# W9 plan REVISION — address dual-review fixes (codex + Opus, both NEEDS-FIXES)

Revise `.orchestration/runs/W9-CUSTODIAL-SCAN/plan_codex.md` IN PLACE to resolve ALL items below.
Both reviewers verified grid/seed reuse is byte-exact vs `scan_plan.json` (keep that). Re-verify
each fix against the REAL repo. End the revised plan with `PLAN-READY`.

## From codex review
1. Replace the generic W8 loop-flag language with the EXACT post-W8 contract (see W8 plan
   `.orchestration/runs/W8-CUSTODIAL-PR2/plan_codex.md` loop section), or explicitly mark loops
   deferred for this scan. Do not assume a single `include_top_partner_loops=True` boolean.
2. Comparison pairing must match REAL rows: the raw field is `seed` NOT `draw_seed`
   (`scripts/run_full_catalog_scan.py:495-498,602-605`). State the builder normalizes
   `row["seed"] -> draw_seed`, `params.M_KK/1000 -> mkk_tev`, and r from top-level `quark_fit_r`
   (`_build_cache.py:42-43`).
3. Add explicit byte-identity tests for the minimal QUARK-ONLY config (not only the full/default
   hash). Baseline WQ rows use config_hash `c6939cc65d71f86a`; resume depends on
   `_completed_summary(... config_hash ...)` (`scripts/run_full_catalog_scan.py:1468-1483`).
4. Expand harness tests to cover BOTH `_evaluate_draw` paths — quark-only AND full-catalog — since
   both call sites change (`scripts/run_full_catalog_scan.py:465-479,570-586`).
5. Fix the isolation section to the REAL W7/W8 overlap: W8 touches `point_builder.py`,
   `rs_ew_builder.py`, `EW001.py`, `T010.py`, `T011.py`, `T014.py`; W7 touches `point_builder.py`,
   `rs_ew_builder.py`. Sequence W9's harness edits AFTER W7+W8 commit to avoid clobber.
6. Fix the `paired_vetoes.parquet` enum contradiction (`veto_class` says rigorous|proxy then adds
   `not_evaluated`): make the enum + nullability explicit.

## From Opus review
7. **Existing PR1 hash test IS the byte-identity guard.** `tests/test_rs_ew_custodial_pr1.py:137-138`
   asserts `"ew_model" not in harness._config_payload(cfg)` AND pinned `_config_hash(cfg) ==
   "45e21a07585f7489"`. Adding `ew_model: str = MINIMAL_RS_EW_MODEL` to `ScanConfig` makes
   `asdict` emit it unconditionally → this test BREAKS unless `_config_payload` POPS `ew_model`
   when `== "minimal_rs"` (`scripts/run_full_catalog_scan.py:1486-1491`). State explicitly: this
   test is the guard, the pinned hash `45e21a07585f7489` must stay unchanged, the pop logic keeps
   it green.
8. **Worker round-trip test.** Workers reconstruct config via `_GLOBAL_CFG =
   _config_from_payload(cfg_payload)` (`:226`), which does `ScanConfig(**values)` (`:1494-1497`).
   Add a test that `_config_from_payload(_config_payload(cfg))` preserves `ew_model` for BOTH modes
   across the multiprocessing boundary.
9. **Spectrum/ew_model hard co-edit.** `_validate_injected_spectrum`
   (`flavor_catalog_constraints/rs_ew_builder.py:304-305`) RAISES if `spectrum.model_label !=
   ew_model`. State that `model_label=cfg.ew_model` and `ew_model=cfg.ew_model` MUST change
   together at EVERY build site (incl. universal-c sanity `:1200-1206,1261-1277`), and add a test
   that a mismatched pair raises.
10. **Universal-c sanity model policy (DECIDE + STATE).** Decide explicitly whether the once-per-
    process universal-c sanity guard runs under the scan's `ew_model` or stays pinned to
    `minimal_rs`, with a one-line justification (custodial proxy tags could change what "no
    spurious HARD exclusions" means / affect `--skip-sanity`). Do not leave it as "same for
    universal-c sanity."

End the revised plan with `PLAN-READY`.
