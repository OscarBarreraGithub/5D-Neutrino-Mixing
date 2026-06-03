Plan: inspect, rewire adapters/wrappers, replace tests, run full suite + dual review.
Rewired IDs: B022/B023 via `rare_b_nunu`; K004/K005 via `rare_kaon`.
Touched wrappers/adapters/tests plus new `tests/constraints/primary/nunu_phase4d_helpers.py`.
K004 rigorous example: `X_NP_total=8.6008647e-4-1.1326867e-4j`.
K004 BR shifts `8.4725984e-11 -> 4.9286100e-12`, ratio `4.927`, vetoing.
SM-limit universal-c: `X_NP=0`; K004 recovers `8.4725984e-11`.
Absent path: `evaluated=False`, `missing_extra=rs_semileptonic_wilsons`, non-vetoing.
Majorana=Dirac active-ν rates checked equal; rate factor `1.0`.
Old proxy NEEDS-HUMAN resolved for these νν channels; no `_wilson_prefactor`, no second `1/M_KK^2`.
Focused test count changed `40 -> 47` (+7).
Pytest: focused `47 passed`; full `1664 passed, 1 skipped`.
Dual review: Codex APPROVE; second available reviewer APPROVE (Opus not exposed by tool stack).
P4D-DONE.