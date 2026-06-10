Reviewed against the real [scripts/build_wq_quarkonly_comparison.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scripts/build_wq_quarkonly_comparison.py:115).

No blocking issues found. The revised plan genuinely incorporates the four fixes:

1. `draw_seed_group_mismatches`: preserves current `rows.values()` dict-comprehension collapse by using retained-row ordinal and greatest ordinal wins.
2. Manifest grid fallback: explicitly replaces `minimal.rows.values()` fallback with retained minimal SQLite `r`/`mkk_tev` values when scan-plan fields are absent.
3. SQLite key: specifies typed `REAL, REAL, INTEGER` composite primary key, avoiding stringified float keys.
4. Manifest metadata: preserves `root`, `ew_model`, `scan_plan`, sorted row+summary `config_hashes`, and sorted row+summary `git_shas` after `RunRows.rows` removal.

Option (b) remains schema-identical, memory-bounded via streaming/chunked Parquet, and aggregate-identical by reusing the existing counter/finalizer semantics.

VERDICT: APPROVE