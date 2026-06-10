Findings: none.

Verified against the real diff: memory is bounded via temp SQLite typed tables and chunked `ParquetWriter` with `PARQUET_CHUNK_ROWS = 50_000`; no full `RunRows.rows`, `paired_rows`, or `veto_rows` retention remains. `paired_vetoes` columns and PyArrow schema are unchanged. Survival and constraint-veto predicates/finalizers are preserved, including constraint sorting and `hard_not_evaluated` order. Composite key is `(r REAL, mkk_tev REAL, draw_seed INTEGER)`, duplicate key message is unchanged, draw-seed last-row-wins is reproduced with retained-row ordinals, and manifest grid fallback no longer depends on `minimal.rows.values()`.

Requested test passed:

```text
..                                                                       [100%]
2 passed in 3.80s
```

VERDICT: APPROVE