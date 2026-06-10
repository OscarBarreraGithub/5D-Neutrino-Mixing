**Findings**

No blocking issues in the plan.

I verified it against [scripts/build_wq_quarkonly_comparison.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scripts/build_wq_quarkonly_comparison.py:62). The current code does load both full runs, retains normalized rows including full `constraints`, materializes `paired_rows` and `veto_rows`, then bulk-converts Python lists into Arrow tables. The plan directly addresses all of those with a disk-backed join index plus chunked `ParquetWriter`.

Option (b) preserves the `paired_vetoes.parquet` design: same per-draw long format, same columns/schema, same model order, same per-row constraint ordering, and only Parquet row-group boundaries change. Survival and veto aggregates should remain identical if the implementation reuses the current predicates in `_accumulate_survival` and `_accumulate_constraint_vetoes`.

Pairing correctness is covered: the plan keeps normalized `(r, mkk_tev, draw_seed)`, including raw `seed -> draw_seed`, `params.M_KK / 1000`, and top-level `quark_fit_r` fallback semantics. Determinism is preserved by ordering the validated join on typed `(r, mkk_tev, draw_seed)` and emitting minimal veto rows before custodial veto rows.

Two implementation details should stay explicit during coding: make the SQLite key lossless for floats, or use a typed composite key, and keep manifest grid fallback behavior without relying on `minimal.rows.values()` after `RunRows.rows` is removed. These are implementation guardrails, not plan blockers.

Verification run: `python -m pytest -q tests/test_wq_quarkonly_comparison.py` passes (`1 passed`). No physics change is implied by the plan; it changes artifact construction strategy only.

VERDICT: APPROVE