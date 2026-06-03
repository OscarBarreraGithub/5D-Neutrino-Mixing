Plan: rewire C004/C005/C007 through `c_to_u_ll`, preserve charm LD diagnostics, replace proxy tests.

Rewired IDs/files: C004/C005/C007 in `primary/charm`; adapters `rare_charm_dilepton.py`, `rare_charm_semileptonic.py`; tests plus `tests/rare_charm_phase3d_helpers.py`.

Rigorous example: sample RS-EW C007 has C9+C9′=`27.378+68.050j`, C10+C10′=`-364.459-905.881j` -> BR=`3.2327e-9`.

SM-limit: universal-c point recovers C004=`0`, C005=`0`, C007=`1.0e-12`.

Absent path: old quark-only point gives non-veto `passes=True`, `predicted=None`, `evaluated=False`, `missing_extra="rs_semileptonic_wilsons"`.

C006/C008 and LFV adapters/tests untouched; C007 full-q²/LD/form-factor PARTIAL diagnostics preserved.

Coverage: C004 `10->10`, C005 `12->12`, C007 `11->11`; proxy-only checks replaced by rigorous/SM/absent/no-second-MKK checks.

Reviews: Codex APPROVE; Opus APPROVE.

Tests: targeted `33 passed`; full `python -m pytest tests/ -q` -> `1646 passed, 1 skipped`; `git diff --check` clean.

P3DC-DONE.