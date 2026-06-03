Plan executed: inspect EW003/adapter/tests, implement diagnostic-only CC wiring, verify.

Changed EW003 to keep the scalar inclusive-vs-exclusive pull data-level while optionally adding `rs_charged_current` diagnostics for `epsilon_cb` and `epsilon_ub`.

Added adapter status/note: PARTIAL/data-level, covariance/scheme NEEDS-HUMAN, no naive NP veto.

Pull unchanged: `3.0728851183895105` vs 3.0 budget; scalar fields remain real finite floats.

Universal CC cancellation asserted: synthetic `epsilon_cb=epsilon_ub=0.2` gives `|1+epsilon|=1.2` but ratio/pull multipliers `1.0`.

EW003 tests increased `8 -> 12`; focused run: `12 passed`.

Full run: `python -m pytest tests/ -q` -> `1677 passed, 1 skipped`.

Review: Codex subagent APPROVE; Opus unavailable in this session, so no Opus approval claimed.

P5C-DONE