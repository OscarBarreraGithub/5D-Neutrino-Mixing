# W9b — comparison builder MEMORY FIX plan (codex author)

`scripts/build_wq_quarkonly_comparison.py` (committed `9b7bf6a`, dual-approved on a SYNTHETIC tiny
pair) OOMs at ~100 GB on the REAL data (1M minimal + 1M custodial). It wrote `paired_draws.parquet`
(6 MB) then died. Author a fix PLAN (no production code this run). Write it to
`.orchestration/runs/W9-CUSTODIAL-SCAN/memfix_plan.md`, end `PLAN-READY`.

## Diagnosed hotspots (verify against the real file)
1. `_load_run` (line ~115) holds ALL normalized rows for BOTH runs in a dict, each retaining the
   full per-row `constraints` payload — ~1M heavy dicts × 2 runs.
2. `_draw_veto_rows` (line ~305) is called per draw and `veto_rows.extend(...)` builds a Python list
   at per-draw × per-constraint granularity (~tens of millions of dicts) which is then handed to
   `pa.Table.from_pylist` (line ~660) — the dominant blowup.

## Requirements for the fix plan
- Keep the UI deliverable intact: a future interface must still be able to plot minimal-vs-custodial
  survival-vs-M_KK per r, and "which constraint vetoes" per (r, M_KK). Preserve `paired_draws`,
  `survival` (per (r,M_KK) for both runs), the manifest/README/schema.
- Make memory BOUNDED and independent of total row count (stream raw JSONL; do not retain full
  per-row constraint payloads after extracting what's needed; do not materialize a tens-of-millions
  row Python list).
- Decide and STATE explicitly the `paired_vetoes` treatment — this is the load-bearing call:
  EITHER (a) AGGREGATE to per-(r, M_KK, constraint, model) veto counts/fractions (bounded, and the
  natural granularity for the UI "which constraint vetoes" view), OR (b) keep per-draw granularity
  but stream it to disk via a chunked `pyarrow.parquet.ParquetWriter` (schema-identical, bounded
  memory). Justify the choice; note any schema change vs the committed design as an EXPLICIT,
  honestly-flagged delta (this is a data-engineering choice, surface it for the reviewers).
- Pairing must stay correct: identical (r, M_KK, draw_seed) keys; normalize raw `seed`->draw_seed,
  `params.M_KK/1000`->mkk_tev, top-level `quark_fit_r`->r (as already implemented).
- Determinism unchanged; numbers identical to what a non-OOM run would produce.
- Add a memory/scale guard test (e.g. a few-thousand-row synthetic pair asserting bounded behavior
  and that paired_vetoes is produced) and keep the existing synthetic test green.

## Constraints
Dual-gate: a second codex + an Opus review this plan, then the impl. Cite file:line, state the
paired_vetoes decision clearly, list tests. End with `PLAN-READY`.
