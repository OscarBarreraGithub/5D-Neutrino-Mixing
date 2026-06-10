# W9b Memory-Fix Plan: `scripts/build_wq_quarkonly_comparison.py`

## Scope

This is a plan only. Do not change production code in this run. The implementation target is
`scripts/build_wq_quarkonly_comparison.py`, whose committed comparison builder OOMs on the real
1M minimal + 1M custodial inputs after writing `paired_draws.parquet`.

## Verified Hotspots

- `build_comparison` loads both runs fully before any artifact work
  (`scripts/build_wq_quarkonly_comparison.py:62-63`).
- `_load_run` creates `rows: dict[...]` (`scripts/build_wq_quarkonly_comparison.py:117`) and stores
  every normalized row under its pairing key (`scripts/build_wq_quarkonly_comparison.py:130-135`).
  `_normalize_row` copies the raw `constraints` payload
  (`scripts/build_wq_quarkonly_comparison.py:181`) and retains it in each stored row
  (`scripts/build_wq_quarkonly_comparison.py:197`). This keeps roughly all heavy per-row constraint
  dictionaries for both runs resident.
- `_build_artifact_rows` materializes both `paired_rows` and `veto_rows`
  (`scripts/build_wq_quarkonly_comparison.py:250-251`), then appends per-draw veto records with
  `veto_rows.extend(_draw_veto_rows(...))` (`scripts/build_wq_quarkonly_comparison.py:290-295`).
- `_draw_veto_rows` expands each row's constraint payload into long-form per-draw veto dictionaries
  (`scripts/build_wq_quarkonly_comparison.py:305-337`), including sorted HARD rigorous/proxy failures
  (`scripts/build_wq_quarkonly_comparison.py:308-318`) and explicit `hard_not_evaluated` gaps
  (`scripts/build_wq_quarkonly_comparison.py:319-336`).
- `_write_paired_draws` and `_write_paired_vetoes` both build a full PyArrow table from Python lists
  (`scripts/build_wq_quarkonly_comparison.py:660` and
  `scripts/build_wq_quarkonly_comparison.py:684`). The `paired_vetoes` path is the dominant blowup
  because it duplicates a tens-of-millions-row Python list into Arrow memory.
- The existing UI-facing aggregate artifact already exists as `constraint_veto_by_r_mkk.csv`
  (`scripts/build_wq_quarkonly_comparison.py:545-548`) with columns defined at
  `scripts/build_wq_quarkonly_comparison.py:879-907`.
- The current synthetic test reads `paired_vetoes.parquet` as per-draw long-form rows and asserts
  `rigorous`, `proxy`, and `not_evaluated` records
  (`tests/test_wq_quarkonly_comparison.py:79-98`).

## `paired_vetoes` Decision

Choose option (b): keep `paired_vetoes.parquet` at per-draw granularity and stream it with a chunked
`pyarrow.parquet.ParquetWriter`.

Justification:

- This is schema-identical to the committed design. `PAIRED_VETOES_COLUMNS`
  (`scripts/build_wq_quarkonly_comparison.py:838-851`) and the Parquet schema in
  `_write_paired_vetoes` (`scripts/build_wq_quarkonly_comparison.py:668-682`) stay unchanged.
- There is no artifact/schema delta for reviewers to absorb. The UI can continue using
  `constraint_veto_by_r_mkk.csv` for the natural per-`(r, M_KK, constraint, model)` "which
  constraint vetoes" view, while `paired_vetoes.parquet` remains the exact per-draw audit/drilldown
  table.
- Memory becomes bounded by the configured Parquet chunk size rather than by total veto row count.
  The on-disk `paired_vetoes.parquet` may still be large; that is an expected disk-size cost, not a
  resident-memory cost.

Explicit schema delta: none for `paired_vetoes.parquet`, `paired_draws.parquet`,
`survival_by_r_mkk.csv`, `constraint_veto_by_r_mkk.csv`, `manifest.json`, `schema.json`, README, or
`run_index.json`.

## Implementation Plan

1. Replace full in-memory `RunRows.rows` with streaming metadata plus a disk-backed compact index.
   Keep run-level metadata (`root`, `ew_model`, `config_hashes`, `git_shas`, `scan_plan`,
   row counts, and manifest grid fallback values) but do not retain normalized row dictionaries in
   Python memory.

2. Add a row-projection path around the existing JSONL normalization:
   - Continue streaming raw files from `_jsonl_paths`
     (`scripts/build_wq_quarkonly_comparison.py:580-585`) line by line.
   - Preserve existing normalization exactly: raw `seed -> draw_seed`, `params.M_KK / 1000 ->
     mkk_tev`, and top-level `quark_fit_r -> r` with `params.quark_fit_r` only as fallback
     (`scripts/build_wq_quarkonly_comparison.py:167-180`).
   - For each raw row, extract only:
     - the pairing key `(r, mkk_tev, draw_seed)`,
     - compact paired-draw fields needed by `PAIRED_DRAWS_COLUMNS`,
     - survival flags needed by `_accumulate_survival`
       (`scripts/build_wq_quarkonly_comparison.py:364-373`),
     - compact per-draw veto records currently produced by `_draw_veto_rows`,
     - aggregate counter contributions currently handled by `_accumulate_constraint_vetoes`
       (`scripts/build_wq_quarkonly_comparison.py:376-397`).
   - Discard the raw row and full `constraints` mapping immediately after these projections are
     produced. No full per-row constraint payload should survive beyond the current line.

3. Use a temporary SQLite database, placed under the output directory or a sibling temp directory, as
   the pair-validation and deterministic-join index:
   - Create one compact table per run keyed by typed columns, e.g. `r REAL NOT NULL`,
     `mkk_tev REAL NOT NULL`, `draw_seed INTEGER NOT NULL`, with
     `PRIMARY KEY (r, mkk_tev, draw_seed)`. Do not use a stringified/canonical float key for pairing,
     uniqueness, joins, or ordering; insert the normalized Python float values directly into typed
     SQLite columns so the `(r, mkk_tev, draw_seed)` key is lossless relative to the committed
     in-memory tuple key.
   - Add a monotonic ingest/retained-row ordinal column per run. Assign it only to rows that are
     retained after the exact duplicate pairing-key check, matching the insertion order of the
     current `rows` dict built from sorted `_jsonl_paths` and line order.
   - Store compact JSON payloads for paired-draw fields and compact veto rows, not raw constraints.
   - Treat composite primary-key insertion conflicts as duplicate pairing keys, matching the current
     duplicate-key validation intent at `scripts/build_wq_quarkonly_comparison.py:132-137`: count the
     conflict, discard the later exact duplicate row, and raise
     `ComparisonValidationError(f"{root}: duplicate pairing keys: {duplicate_count}")` after the
     stream completes.
   - Track row counts, row-level `config_hash`/`git_sha`, and summary-level `config_hash`/`git_sha`
     without scanning a retained row dictionary (`scripts/build_wq_quarkonly_comparison.py:141-156`).
   - Track manifest grid fallback candidates from the retained minimal rows by typed `r` and
     `mkk_tev` values so `_manifest_payload` can preserve the current fallback behavior at
     `scripts/build_wq_quarkonly_comparison.py:482-484` after `RunRows.rows` is removed.

4. Validate pairing before final artifact publication:
   - Keep the existing scan-plan equivalence check
     (`scripts/build_wq_quarkonly_comparison.py:64-66` and
     `scripts/build_wq_quarkonly_comparison.py:610-611`).
   - Compute `minimal_rows`, `custodial_rows`, `paired_rows`, `minimal_only`, and `custodial_only`
     with SQLite joins/anti-joins on the typed composite key:
     `minimal.r = custodial.r AND minimal.mkk_tev = custodial.mkk_tev AND
     minimal.draw_seed = custodial.draw_seed`.
   - Compute `draw_seed_group_mismatches` by reproducing the current seed-group dict-comprehension
     collapse exactly (`scripts/build_wq_quarkonly_comparison.py:211-223`): for each run, reduce
     retained rows to one `(r, mkk_tev)` per `draw_seed`, and if the same `draw_seed` appears under
     multiple `(r, mkk_tev)` values in one run, the row with the greatest retained-row ordinal wins.
     Then compare only the intersection of minimal/custodial seed sets and count seeds whose winning
     `(r, mkk_tev)` tuples differ. This preserves the current "last row wins" byte-identical
     mismatch count, including the case where repeated seeds are not duplicate pairing keys.
   - If any mismatch exists, raise `ComparisonValidationError` before replacing output artifacts.

5. Stream artifact generation from the validated SQLite join:
   - Iterate `minimal JOIN custodial ON minimal.r = custodial.r AND
     minimal.mkk_tev = custodial.mkk_tev AND minimal.draw_seed = custodial.draw_seed`, ordered by
     the typed minimal `(r, mkk_tev, draw_seed)` columns. This preserves the current deterministic
     sorted-key output order from `scripts/build_wq_quarkonly_comparison.py:255`.
   - For each joined pair, build one `paired_draws` row and write it to a bounded Parquet chunk buffer.
   - For vetoes, emit the minimal row's compact veto records followed by the custodial row's compact
     veto records, preserving the current model order from
     `scripts/build_wq_quarkonly_comparison.py:290-295`. Preserve the existing per-row constraint
     sort inside veto projection (`scripts/build_wq_quarkonly_comparison.py:308`) and the existing
     `hard_not_evaluated` list order (`scripts/build_wq_quarkonly_comparison.py:319`).
   - Accumulate `survival` and `constraint_veto` using the same counters and finalizers, so output
     counts/fractions remain identical to a non-OOM run.

6. Replace bulk Parquet writes with chunked writers:
   - Refactor `_write_paired_draws` and `_write_paired_vetoes` into streaming helpers that accept an
     iterable of row mappings and a fixed chunk size, e.g. `PARQUET_CHUNK_ROWS = 50_000`.
   - Use `pyarrow.parquet.ParquetWriter` with the existing schemas and write `pa.Table.from_pylist`
     only for the current chunk.
   - Ensure an empty `paired_vetoes.parquet` is still created with the correct schema if a valid run
     has zero veto rows.
   - Keep atomic output behavior through `_atomic_write`
     (`scripts/build_wq_quarkonly_comparison.py:710-718`) so failures leave no partial final files.

7. Keep bounded in-memory aggregates:
   - `survival` remains keyed only by `(r, mkk_tev)`.
   - `constraint_veto` remains keyed only by `(ew_model, r, mkk_tev, constraint_id, tag, severity)`.
   - These structures scale with grid size and constraint catalog size, not with total draw count.

8. Preserve manifest, README, schema, and run index behavior:
   - Keep artifact names at `scripts/build_wq_quarkonly_comparison.py:74-82`.
   - Keep manifest run metadata semantics from `scripts/build_wq_quarkonly_comparison.py:467-480`:
     `git_sha` still comes from `_single_or_list` over sorted row/summary SHAs, and
     `config_hashes` remains the sorted union of row-level and summary-level hashes.
   - Preserve manifest grid behavior from `scripts/build_wq_quarkonly_comparison.py:481-486` without
     `minimal.rows.values()`: use `minimal.scan_plan["r_grid"]` and `minimal.scan_plan["mkk_tev"]`
     when present; otherwise fall back to `_sorted_unique` over the retained minimal rows' typed
     SQLite `r` and `mkk_tev` columns, respectively. This is the direct replacement for the current
     `minimal.rows.values()` fallback and keeps the manifest working when `scan_plan.json` omits one
     or both grid fields.
   - Keep schema descriptions for all current artifacts, especially `paired_vetoes.parquet`
     (`scripts/build_wq_quarkonly_comparison.py:537-540`).
   - Keep README pairing/normalization statements at
     `scripts/build_wq_quarkonly_comparison.py:564-567`.

## Determinism and Numerical Equivalence

- Output row order remains deterministic by sorting the SQLite join on `(r, mkk_tev, draw_seed)`,
  matching the current sorted iteration over `minimal.rows`.
- Pairing remains exact on identical normalized keys `(r, mkk_tev, draw_seed)`.
- Survival and veto aggregate numbers remain identical because the same row flags, constraints
  predicates, `_safe_div`, and finalizer logic are reused.
- Per-draw `paired_vetoes` content remains identical except for Parquet row-group boundaries, which
  are not part of the logical schema.

## Tests

- Keep the existing synthetic comparison test green:
  `tests/test_wq_quarkonly_comparison.py::test_wq_quarkonly_comparison_pairs_seed_normalizes_survival_and_veto_enums`.
- Add a scale/memory guard test in `tests/test_wq_quarkonly_comparison.py`:
  - Generate a few thousand paired synthetic rows with repeated rigorous, proxy, and
    `hard_not_evaluated` vetoes.
  - Force a small Parquet chunk size by monkeypatching the module constant, e.g. 64 or 128 rows.
  - Assert `paired_draws.parquet`, `paired_vetoes.parquet`, `survival_by_r_mkk.csv`,
    `constraint_veto_by_r_mkk.csv`, `manifest.json`, `schema.json`, README, and `run_index.json` are
    produced.
  - Assert the `paired_vetoes.parquet` row count equals the expected synthetic veto count.
  - Assert no Parquet write batch exceeds the configured chunk size, either by instrumenting the
    streaming helper or monkeypatching the writer boundary in the test. This directly guards against
    reintroducing a full-list `pa.Table.from_pylist` path.
- Add or extend a validation test for duplicate/missing keys if the SQLite validation code changes
  the shape of `ComparisonValidationError` messages.
- Run targeted verification:
  `python -m pytest -q tests/test_wq_quarkonly_comparison.py`.
- If the implementation touches shared scan or schema helpers unexpectedly, also run the previously
  approved W9 targeted set before review.

## Reviewer Notes

- The fix trades resident memory for temporary disk use. The SQLite index and final
  `paired_vetoes.parquet` can grow with row/veto count, but Python memory should remain bounded by
  SQLite page cache, the aggregate dictionaries, and the fixed Parquet chunk buffers.
- The load-bearing data-engineering choice is explicit: preserve per-draw `paired_vetoes` and stream
  it. Reviewers who prefer a smaller aggregate-only artifact should treat that as a schema/product
  decision, not as a memory requirement.
- No production code has been changed by this plan-authoring run.

PLAN-READY
