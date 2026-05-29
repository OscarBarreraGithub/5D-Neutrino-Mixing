1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation OK: `T020.py:46` reaches physics only via `physics_adapters.higgs_lfv`; AST check found 0 removed/changed/added funcs in `physics_adapters/higgs_lfv.py`, `quarkConstraints/higgs_lfv.py`, `physics_adapters/deltaf2.py`.
5. Contract OK: real float result fields at `T020.py:435-450`; complex Yukawas stay in diagnostics at `T020.py:416-417`; missing input returns non-crashing unevaluated result at `T020.py:382-385`.
6. Anchor OK: YAML CMS limit is routed through `load_anchor` at `T020.py:230`; missing/mismatched probes raised `AnchorError` loudly.
7. Independent numerical cross-check: limit `4.4000000000000006e-05`; for `Y_e_mu=1.0e-4+0.4e-4j`, `Y_mu_e=0.5e-4j`, constraint/core/manual BR all `1.7264840809205439e-05`, ratio `0.39238274566375991`, pass `True`.
8. Threshold checks: `Y=1.0e-4` gives BR `1.2244567949791087e-05`, ratio `0.2782856352225247`, pass `True`; `Y=2.0e-4` gives BR `4.897827179916435e-05`, ratio `1.1131425408900988`, pass `False`.
9. Determinism OK: repeated `evaluate()` returned equal results and did not mutate `ParameterPoint.extras`.
10. Registry smoke OK: `fcc.get("T020") -> T020 top_higgs_ew HARD`; registry count `64`.
11. Tests: `python -m pytest tests/constraints/ -q` -> `678 passed in 16.33s`.
CODE-OK