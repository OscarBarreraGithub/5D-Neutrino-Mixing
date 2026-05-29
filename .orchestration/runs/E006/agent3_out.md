1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation OK: only E006 constraint/test plus new `mercury_edm.py` adapter are added; `deltaf2.py` and `E006.yaml` unchanged. E006 imports no other constraint and reaches physics only via adapter.
5. Contract OK: `E006.py:440,453-461` returns INFO, real float `experimental/ratio/budget`, `predicted=None`, `sm_prediction=None`; empty/missing-coupling path is non-crashing.
6. Anchor OK: `E006.py:247-251` uses `load_anchor`; missing/mismatched probes raise `AnchorError`.
7. Numerical cross-check: YAML limit `7.40000000000000062e-30`, central `2.20000000000000053e-30`, total uncertainty `3.12296333632016914e-30`, ratio `0.29729729729729737`; evaluate ratio matches exactly.
8. Pass/fail check: adapter safe ratio `0.5` passes; over-limit ratio `2.0` fails. E006 itself remains INFO/non-vetoing.
9. Determinism OK: repeated `evaluate()` returns identical result and leaves point extras `{}` unchanged.
10. Tests: focused E006 `8 passed`; full `python -m pytest tests/constraints/ -q` -> `810 passed in 16.37s`.
11. Registry smoke: `79` constraints, E006 registered as `PRIMARY INFO |d_Hg|`, `import_failures {}`.

CODE-OK