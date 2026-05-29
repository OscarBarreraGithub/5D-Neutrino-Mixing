1. BLOCKER: none found.
2. SHOULD-FIX: none found.
3. NIT: none found.
4. Isolation OK: [C004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charm/C004.py:49) reaches physics only through the rare-charm adapter; test core import is only in [test_C004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/charm/test_C004.py:17). No C004 edits to pre-existing tracked functions observed.
5. Contract/anchor OK: HARD severity, real float result fields, complex diagnostics only, missing extra returns non-crashing pass. Anchor probes raised `AnchorError` for missing candidate, missing value, and unit mismatch.
6. Numerical cross-check: `SM_SD constraint=0`, `manual=0`; YAML budget `2.1000000000000002e-09`; LD VMD context `9.45e-13`.
7. Safe point: `left_uc=1e-2`, `predicted=direct=1.3570772436699916e-13`, `ratio=6.4622725889047213e-05`, passes `True`.
8. Excluded point: `left_uc=5`, `predicted=direct=3.3926931091749786e-08`, `ratio=16.155681472261801`, passes `False`.
9. Determinism OK: repeated `evaluate()` results equal; input coupling matrices unchanged.
10. Tests: `python -m pytest tests/constraints/primary/charm/test_C004.py -q` -> `10 passed`; `python -m pytest tests/constraints/ -q` -> `244 passed`.
11. Registry smoke OK: C004 registered as `charm/HARD`; `import_failures=0`.

CODE-OK