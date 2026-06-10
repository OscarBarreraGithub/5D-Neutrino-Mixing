# W9b — IMPLEMENT the DUAL-APPROVED comparison-builder memory fix

The memory-fix plan `.orchestration/runs/W9-CUSTODIAL-SCAN/memfix_plan.md` is DUAL-APPROVED. Implement
it EXACTLY in `scripts/build_wq_quarkonly_comparison.py` (+ tests). If infeasible, STOP and write
`.orchestration/runs/W9-CUSTODIAL-SCAN/memfix_impl_blocker.md`.

## Non-negotiables (from the approved plan)
- Peak memory BOUNDED, independent of total row count: stream raw JSONL; do NOT retain full per-row
  `constraints` payloads for both runs; do NOT materialize a tens-of-millions-row Python list.
- Option (b): per-draw `paired_vetoes.parquet`, SCHEMA-IDENTICAL (`PAIRED_VETOES_COLUMNS` and the
  `_write_paired_vetoes` schema UNCHANGED), streamed via `pyarrow.parquet.ParquetWriter` with a
  module-level `PARQUET_CHUNK_ROWS = 50_000`. Empty-with-schema case preserved. `_atomic_write` kept.
- Disk-backed SQLite join index with a TYPED composite primary key `(r REAL, mkk_tev REAL,
  draw_seed INTEGER)` — lossless, no stringified floats.
- Preserve EXACT dedup-by-`draw_seed` last-row-wins collapse so `draw_seed_group_mismatches` is
  byte-identical; the duplicate `ComparisonValidationError` message unchanged.
- Manifest grid fallback re-sourced from SQLite distinct (or scan_plan) — no `minimal.rows.values()`.
  config_hashes/git_shas semantics preserved.
- Survival + constraint-veto aggregates IDENTICAL (reuse the existing counters/predicates/finalizers,
  same constraint sort + hard_not_evaluated order). Determinism unchanged. NO physics change.

## Tests
- Keep the existing synthetic test green.
- Add the scale/memory guard: monkeypatch a small `PARQUET_CHUNK_ROWS`, assert all artifacts produced,
  paired_vetoes row count == expected, and NO Parquet write batch exceeds the configured chunk size.

## After implementing
```
source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
python -m pytest -q tests/test_wq_quarkonly_comparison.py 2>&1 | tail -15
```
Write `.orchestration/runs/W9-CUSTODIAL-SCAN/memfix_impl_summary.md` (what changed, tests, pass count)
ending `IMPL-READY` (green) or `IMPL-BLOCKED`. Do NOT commit and do NOT run the full 2M build.
