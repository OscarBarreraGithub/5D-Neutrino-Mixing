# W9 Custodial Scan Implementation Summary

## Files Changed

- `scripts/run_full_catalog_scan.py`
  - Added `--ew-model {minimal_rs,custodial_rs_plr}`.
  - Added `ScanConfig.ew_model` defaulting to `minimal_rs`.
  - Preserved default/minimal config payload byte surface by omitting `ew_model` when it is `minimal_rs`.
  - Threaded `cfg.ew_model` into every scan and universal-c spectrum/builder site with matching `model_label=cfg.ew_model`.
  - Added non-default custodial `ew_model` row, summary, sanity, and provenance metadata.
- `scripts/build_wq_quarkonly_comparison.py`
  - Added paired minimal-vs-custodial comparison builder.
  - Normalizes raw `seed -> draw_seed`, `params.M_KK / 1000 -> mkk_tev`, and top-level `quark_fit_r -> r`.
  - Validates identical scan plans except run identity fields and exact `(r, mkk_tev, draw_seed)` pairing.
  - Writes `README.md`, `manifest.json`, `schema.json`, `paired_draws.parquet`, `paired_vetoes.parquet`, `survival_by_r_mkk.csv`, `constraint_veto_by_r_mkk.csv`, and `run_index.json`.
- `scripts/run_wq_quarkonly_1m_custodial_array.sbatch`
  - Added custodial WQ 1M quark-only array wrapper.
  - Reuses `scripts/wq_quarkonly_1m_plan.py` unchanged for grid and seeds.
  - Uses output root `scan_outputs/wq_quarkonly_1M_custodial_<jobid>/`.
  - Adds `--ew-model custodial_rs_plr --quark-only`; no W8 loop/FCNC flags.
- `tests/test_full_catalog_scan_harness.py`
  - Added pinned WQ quark-only hash guard.
  - Added parser/config payload round-trip coverage for both EW models.
  - Added custodial plumbing coverage for quark-only and full `_evaluate_draw`.
  - Added universal-c sanity propagation coverage.
- `tests/test_rs_ew_custodial_pr1.py`
  - Added injected spectrum/model mismatch guard.
- `tests/test_wq_quarkonly_comparison.py`
  - Added synthetic paired-run comparison test covering seed normalization, skipped rows, survival aggregates, and rigorous/proxy/not-evaluated veto enum output.

## Verification

- Targeted tests: `21 passed in 11.71s`.
- Full suite command:
  `source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH" && python -m pytest -q 2>&1 | tee /tmp/w9_pytest_full.log | tail -30`
- Full suite result: `1752 passed, 1 skipped in 866.48s (0:14:26)`.
- Pinned hashes confirmed unchanged:
  - Default/minimal full hash: `45e21a07585f7489`
  - Canonical WQ quark-only `r=0.05`, shard-00 hash: `c6939cc65d71f86a`
  - `ew_model` absent from both default minimal payloads.

## Custodial Run Path

- Custodial sbatch path: `scripts/run_wq_quarkonly_1m_custodial_array.sbatch`
- Intended run root: `scan_outputs/wq_quarkonly_1M_custodial_<jobid>/`
- No SLURM job was submitted.

IMPL-READY
