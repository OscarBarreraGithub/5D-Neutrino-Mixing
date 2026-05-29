1. No BLOCKER / SHOULD-FIX / NIT findings for C003 code/tests.
2. Isolation: C003 adds only `flavor_catalog_constraints/primary/charm/C003.py`, `flavor_catalog_constraints/physics_adapters/charm_direct_cp.py`, and `tests/constraints/primary/charm/test_C003.py`; AST check found no other-constraint imports and physics access only via the adapter.
3. Contract/anchor probe: result fields are real floats or `None`; empty/missing-couplings point does not crash; `load_anchor` calls hit LHCb/HFLAV anchors; forced missing and mismatched anchors raise `AnchorError`.
4. Cross-check numbers: YAML `-15.4 * 1e-4 = -0.0015400000000000001`, uncertainty `2.9e-4 = 0.00029`, budget `0.0015400000000000001`, ratio `1.0`; constraint output matches exactly.
5. Pass/fail scalar probe: room `1.54e-3` gives ratio `1.0`, pass `True`; room `0.77e-3` gives ratio `2.0`, pass `False`. Point-level exclusion is intentionally absent for this INFO non-vetoing stub.
6. Determinism: repeated `evaluate()` results equal, and `ParameterPoint.extras` unchanged.
7. Tests: `python -m pytest tests/constraints/primary/charm/test_C003.py -q` -> `8 passed`; registry smoke -> `69` constraints, C003 `family=charm`, `severity=INFO`.
8. Full requested run: `python -m pytest tests/constraints/ -q` -> `724 passed in 15.73s`.

CODE-OK