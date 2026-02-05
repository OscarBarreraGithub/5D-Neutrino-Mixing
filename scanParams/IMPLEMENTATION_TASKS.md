# Scanner Buildout Task Log

Date: 2026-02-05
Branch: `codex/scanner-theory-priors`

This file tracks atomic implementation tasks for the scanner refactor.

## Planned Tasks

- [x] T1: Establish implementation checklist and task boundaries.
- [x] T2: Refactor `ScanConfig` and scan loops for explicit `k`, `xi_KK`, `MN_mode`, `MN_over_k` and multi-parameter grids.
- [x] T3: Add LFV BR-driven conversion and reproducibility metadata fields in scan output.
- [x] T4: Add anarchic Yukawa sampling/scoring module and integrate optional scoring into scan rows.
- [x] T5: Expand tests for config validation, metadata, KK/LFV behavior, and anarchy determinism.
- [x] T6: Update scanner docs (`README`, theory note references if needed) and run full test suite.
- [x] T7: Replace fixed `% c_L` degeneracy assumption with wavefunction-ratio prior (`max_fL_ratio`) and derived `delta_cL_max` metadata.

## Completion Log

- 2026-02-05: T1 completed.
- 2026-02-05: T2 completed (`scanParams/scan.py` refactor and new scan metadata schema).
- 2026-02-05: T3 completed (LFV BR-driven conversion, `M_KK` override path, metadata fields).
- 2026-02-05: T4 completed (`scanParams/anarchy.py` and optional score integration).
- 2026-02-05: T5 completed (`tests/test_scan.py`, `tests/test_mu_to_e_gamma.py` expanded).
- 2026-02-05: T6 completed (`scanParams/README.md` updated, full tests pass: `42 passed`).
- 2026-02-05: T7 completed (`max_fL_ratio` implemented in `ScanConfig`; per-point `delta_cL_max_{symmetric,one_sided}` derived from `f_IR`; docs and tests updated; full tests pass: `45 passed`).
