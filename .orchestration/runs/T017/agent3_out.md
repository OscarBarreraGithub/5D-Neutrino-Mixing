1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation OK: scoped status shows only new T017 file, new `zpole_lfv_mutau` adapter, new T017 test; AST found no imports of other constraints, and T017 reaches physics via adapter at `T017.py:46`.
5. Contract OK: valid result numeric fields are floats; complex couplings stay in diagnostics; missing input returns `passes=True, predicted=None, ratio=None`.
6. Anchor OK: YAML/load_anchor limit `6.5e-06`; missing value_id probe raises `AnchorError`.
7. Independent cross-check: `sm_total=3.649563333333334`, `scale=0.0009239087104177778`, `deltaL=(0.00023097717760444444+9.239087104177778e-05j)`, `deltaR=4.619543552088889e-05j`.
8. BR check: manual `3.508394829491628e-08`, constraint `3.508394829491628e-08`, abs diff `0.0`, ratio `0.005397530506910197`.
9. Pass/fail check: overlap `0.10` -> BR `4.677859914890927e-09`, ratio `0.0007196707561370657`, pass; overlap `20.0` -> BR `0.0001870793922231769`, ratio `28.781444957411836`, fail.
10. Determinism OK: repeated `evaluate()` returned equal results and did not mutate the point.
11. Registry smoke OK: `fcc.get("T017") -> T017 top_higgs_ew HARD`.
12. Pytest: targeted T017 `13 passed`; full `python -m pytest tests/constraints/ -q` -> `463 passed in 14.32s`.

CODE-OK