# W9b Memory Fix Implementation Summary

Implemented the approved comparison-builder memory fix in
`scripts/build_wq_quarkonly_comparison.py`.

- Replaced retained full-run row dictionaries with a temporary SQLite join index using typed
  `(r REAL, mkk_tev REAL, draw_seed INTEGER)` primary keys.
- Stored only compact JSON projections for paired draw fields, per-draw veto rows, and constraint
  aggregate contributions; raw per-row `constraints` payloads are discarded after each JSONL line.
- Preserved duplicate pairing-key handling and the existing duplicate error message.
- Preserved draw-seed group mismatch semantics with retained-row ordinals so repeated
  `draw_seed` values still collapse last-row-wins.
- Re-sourced manifest grid fallbacks from SQLite distinct retained minimal rows.
- Streamed `paired_draws.parquet` and `paired_vetoes.parquet` through
  `pyarrow.parquet.ParquetWriter` with `PARQUET_CHUNK_ROWS = 50_000`, preserving empty-with-schema
  output and the existing schemas.
- Added a scale/memory guard test that monkeypatches a small chunk size, verifies all artifacts,
  checks `paired_vetoes.parquet` row count, and asserts no Arrow write batch exceeds the chunk size.

Targeted verification run:

```text
source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
python -m pytest -q tests/test_wq_quarkonly_comparison.py 2>&1 | tail -15
```

Result: `2 passed in 10.48s`.

The full 2M build was not run.

IMPL-READY
