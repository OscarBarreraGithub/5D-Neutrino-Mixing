Plan: inspect prior pattern, rewire T015-T017 only, verify v1-zero/live-LFV/absent paths.

Rewired: `zpole_lfv.py` adapter export plus `T015`, `T016`, `T017`; tests updated and shared helper added.

T015/T016/T017 now use `z_delta_g_L/R_e[0,1]`, `[0,2]`, `[1,2]` via `z_lfv_branching_fraction_from_couplings`.

v1 diagonal builder: all three BR_NP = `0.0`, rigorous-tree, non-vetoing.

LFV-live toy: nonzero BRs T015 `2.52e-9`, T016 `4.33e-9`, T017 `8.72e-9`; norm scales by `(3000/6000)^4 = 0.0625`.

Absent old-style point: `evaluated=False`, `missing_extra=rs_ew_couplings`, non-vetoing.

Diagnostics keep `PARTIAL/NEEDS-HUMAN`: loop-induced Z-LFV from `lfv_dipole_spurion` deferred to Phase 7.

Test count: proxy-only Z-LFV tests replaced, targeted collection `38 -> 30`.

Verification: targeted `30 passed`; full `python -m pytest tests/ -q` `1645 passed, 1 skipped`; `git diff --check` clean.

Review: Codex subagent APPROVE; true Opus/Claude tooling was not available here, so I did not claim an Opus approval.

P4B-DONE.