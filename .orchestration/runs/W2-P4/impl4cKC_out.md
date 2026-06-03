Plan executed: inspect, rewire via adapters, replace proxy tests, run focused/full verification.

Rewired IDs/files: K019/K020/K021, C006/C008; adapters `rare_kaon_lfv_*`, `rare_charm_lfv_*`; tests plus `tests/constraints/lfv_rare_phase4c_helpers.py`.

v1-zero: all 5 diagonal-fit rates are `0.0`, `ratio=0.0`, evaluated, passing, rigorous-tree non-vetoing.

LFV-live rates: K019 `6.87e-17`, K020 `4.70e-19`, K021 `4.33e-18`, C006 `3.38e-18`, C008 `1.17e-16`; amplified toys fail and MKK scaling is checked.

Absent path: all 5 return non-vetoing `evaluated=False`, `missing_extra=rs_semileptonic_wilsons`.

Old lepton-coupling `NEEDS-HUMAN` flag is resolved on the production path and replaced by the honest tree-level-zero/deferred-loop note.

Test count changed `49 -> 50` across the five scoped test files; focused suite `50 passed`.

Full suite: `1638 passed, 1 skipped` via `python -m pytest tests/ -q`.

Review: Codex APPROVE. Opus review not available in this environment, so I did not claim Opus approval.

P4CKC-DONE.