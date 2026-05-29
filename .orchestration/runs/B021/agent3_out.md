1. No BLOCKER / SHOULD-FIX / NIT findings for B021 code/tests.

Cross-check numbers: SM/manual bin BR `5.532430647e-07` vs constraint `5.532430650e-07`; anchor `5.900000000e-07`; budget `3.458207347e-07`. Safe point BR `5.405642340e-07`, ratio `0.142952001`, pass `True`. Excluded point BR `1.528489909e-06`, ratio `2.713804624`, pass `False`.

Isolation/contract/anchor checks passed: B021 imports physics only via adapters; B021 additions are new files only; missing input returns non-crashing result; anchor missing/mismatch raises `AnchorError`.

Pytest: `python -m pytest tests/constraints/ -q` -> `568 passed in 14.29s`. B021 targeted: `10 passed in 4.18s`. Registry smoke passed for `B021`.

CODE-OK