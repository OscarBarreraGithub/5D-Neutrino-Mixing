Independent review completed. The plan is close on grid/seed reuse and comparison intent, but it is not ready to approve.

The baseline grid and seed contract are correctly identified: `scripts/wq_quarkonly_1m_plan.py:13-20`, `:49-70`, and `scan_outputs/wq_quarkonly_1M_20128400/scan_plan.json` match the requested 5 x 10 x 20k draw contract. Reusing that helper unchanged in the custodial sbatch is the right apples-to-apples approach.

The blockers are in implementation concreteness: PR2 loop integration is still conditional and underspecified, the comparison schema names a key that real rows do not contain, full-path harness coverage is missing, and isolation misses real W7/W8 overlap.

VERDICT: NEEDS-FIXES

1. In `.orchestration/runs/W9-CUSTODIAL-SCAN/plan_codex.md:67-71`, replace the generic W8 loop-flag language with the exact post-W8 contract. W8’s plan requires more than `include_top_partner_loops=True`: see `.orchestration/runs/W8-CUSTODIAL-PR2/plan_codex.md:109-116` and `:154-176`. Specify the final CLI/config fields or explicitly mark loops deferred.

2. In `.orchestration/runs/W9-CUSTODIAL-SCAN/plan_codex.md:148-157`, fix comparison pairing to match real rows: the raw field is `seed`, not `draw_seed` (`scripts/run_full_catalog_scan.py:495-498`, `:602-605`). State that comparison normalizes `row["seed"] -> draw_seed` and `params.M_KK / 1000 -> mkk_tev`.

3. Add explicit byte-identity tests for minimal quark-only configs, not only the full/default hash. Existing hashes are `45e21a07585f7489` for the two-MKK full config and `d96cb734f724aedb` for the analogous quark-only config; baseline WQ rows use `c6939cc65d71f86a`. Resume depends on `_completed_summary(... config_hash ...)` at `scripts/run_full_catalog_scan.py:1468-1483`.

4. Expand harness tests in `.orchestration/runs/W9-CUSTODIAL-SCAN/plan_codex.md:51-57` to cover both quark-only and full-catalog `_evaluate_draw` paths. The plan changes both call sites (`scripts/run_full_catalog_scan.py:465-479`, `:570-586`), but only specifies a tiny quark-only plumbing test.

5. Fix isolation against the real W7/W8 plans. W8 also touches `flavor_catalog_constraints/point_builder.py`, `rs_ew_builder.py`, `EW001.py`, `T010.py`, `T011.py`, and `T014.py` (`.orchestration/runs/W8-CUSTODIAL-PR2/plan_codex.md:251-260`). W7 touches `point_builder.py` and `rs_ew_builder.py` (`.orchestration/runs/W7-MUEGAMMA/plan_codex.md:127-145`). W9’s overlap section currently understates this.

6. Fix the `paired_vetoes.parquet` enum contradiction: `.orchestration/runs/W9-CUSTODIAL-SCAN/plan_codex.md:164-166` says `veto_class` is only `rigorous|proxy`, then proposes `not_evaluated`. Make the enum and nullability explicit.