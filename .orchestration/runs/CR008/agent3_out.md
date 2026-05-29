1. BLOCKER: None.
2. SHOULD-FIX: None.
3. NIT: None.

4. Isolation: PASS. `CR008.py:52` imports only the collider adapter; AST import check found no other constraint/core imports. `CR001-CR007`, adapter, and core tracked diffs are empty.
5. Contract: PASS. Numeric fields are real floats; missing extra returns non-crashing `passes=True`, invalid mass returns non-crashing `passes=False`.
6. Anchor: PASS. `CR008.py:253` routes through `load_anchor`; missing value-id and unit-mismatch probes both raised `AnchorError`.
7. Independent cross-check: PASS. Core `quarkConstraints.collider_resonance` direct math matched constraint output.
8. Numbers: `M_KK=1600 GeV -> pred=1.6 TeV, exp=1.36 TeV, ratio=0.85, pass=True`; `1500 -> ratio=0.906666666667, pass=True`; `1200 -> ratio=1.13333333333, pass=False`; `1360 -> ratio=1.0, pass=True`.
9. Determinism: PASS. Repeated `evaluate()` results identical; point extras unchanged.
10. Tests: `python -m pytest tests/constraints/ -q` -> `862 passed in 17.75s`; CR008 narrow test -> `8 passed in 4.49s`.
11. Registry import smoke: PASS. `fcc.get("CR008")` resolves `flavor_catalog_constraints.primary.collider_rs.CR008`, `HARD`, `collider_rs`.

CODE-OK