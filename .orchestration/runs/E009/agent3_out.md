1. BLOCKER: None.
2. SHOULD-FIX: None.
3. NIT: None.
4. Isolation: E009/adapter import no other constraints; `deltaf2.py` has empty diff; E009 reaches physics only via `physics_adapters.weinberg_operator`.
5. Contract: empty/missing-coupling path returns `E009 INFO passes=True`, `predicted=None`, `sm_prediction=None`, `ratio=None`, `experimental=budget=1.2e-11`.
6. Anchor probes: mismatched `load_anchor` block raises `AnchorError`; missing value id raises `AnchorError`.
7. Independent numbers from YAML: PR `9.12e-13 / 0.022 = 4.1454545454545456e-11`, YAML rounded `4.1e-11`; HH `9.12e-13 / 0.074 = 1.2324324324324325e-11`, YAML rounded `1.2e-11`.
8. Safe/excluded helper check: safe ratio `0.5` passes; excluded ratio `2.0` fails. Constraint remains non-vetoing by design.
9. Determinism: repeated `evaluate()` results equal; `point.extras` not mutated.
10. Registry smoke: `fcc.get("E009")` registered as `edm_neutrino` / `INFO`; no E009 import failures.
11. Pytest: focused E009 `9 passed`; requested `python -m pytest tests/constraints/ -q` -> `959 passed in 20.36s`.

CODE-OK