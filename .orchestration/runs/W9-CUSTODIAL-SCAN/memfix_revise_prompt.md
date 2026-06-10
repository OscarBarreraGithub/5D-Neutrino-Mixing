# W9b memfix plan REVISION — address review fixes (codex APPROVE w/ guardrails; Opus NEEDS-FIXES)

Revise `.orchestration/runs/W9-CUSTODIAL-SCAN/memfix_plan.md` IN PLACE to resolve ALL items below.
Keep option (b) (per-draw paired_vetoes, chunked ParquetWriter, schema-identical). Re-verify against
the real `scripts/build_wq_quarkonly_comparison.py`. End with `PLAN-READY`.

## Opus NEEDS-FIXES (blocking)
1. **Dedup-by-draw_seed semantics.** The current `_validate_pairing` seed-group dicts
   (`build_wq_quarkonly_comparison.py:211-218`) collapse to ONE entry per `draw_seed` (last row
   wins), including when the same `draw_seed` appears under different `(r,mkk_tev)` in one run. The
   SQLite reimplementation MUST reproduce this exact collapse so the `draw_seed_group_mismatches`
   count is byte-identical. State it explicitly.
2. **Manifest grid fallback.** Lines 482-484 build the manifest `grid` via
   `_sorted_unique(row["r"] ...)` / `row["mkk_tev"]` over `minimal.rows.values()` — which iterates
   the now-removed retained rows. The plan must re-source this (SQLite distinct values, or the
   `scan_plan`) so it does not break when `plan.get("r_grid")`/`plan.get("mkk_tev")` is absent.
   (config_hashes/git_shas at 141-156 are already covered; the grid fallback is NOT.)

## codex implementation guardrails (fold into the plan)
3. Make the SQLite join key LOSSLESS for floats — use a typed composite key (typed columns), not a
   stringified float, so `(r, mkk_tev, draw_seed)` pairing cannot drift.
4. Keep all manifest behavior that currently relies on `minimal.rows.values()` working after
   `RunRows.rows` is removed (ties to #2).

End the revised plan with `PLAN-READY`.
