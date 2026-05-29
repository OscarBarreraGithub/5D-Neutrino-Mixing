1. BLOCKER: none found.
2. SHOULD-FIX: none found.
3. NIT: none found.
4. Isolation OK: `CR011.py` imports no sibling constraints and reaches physics only via `physics_adapters.vbs_longitudinal`; `deltaf2.py` is unchanged.
5. Contract OK: outputs are real floats/None; empty/missing raw path returns INFO, `passes=True`, no crash.
6. Anchor OK: active limit loaded from YAML via `load_anchor`; probes for missing active value, bad units, and nonnumeric value all raised `AnchorError`.
7. Cross-check: YAML limit `0.45 fb`; `0.225/0.45 = 0.5`, result ratio `0.5`, pass `True`.
8. Safe/excluded: `0.10/0.45 = 0.2222222222`, pass `True`; `0.90/0.45 = 2.0`, pass `False` with INFO severity.
9. Determinism OK: repeated `evaluate()` results equal; raw input unchanged.
10. Pytest: `tests/constraints/primary/collider_rs/test_CR011.py -q` -> `8 passed`; `tests/constraints/ -q` -> `810 passed`.
11. Registry smoke OK: CR011 registered as `collider_rs INFO`; registry count `79`.

CODE-OK