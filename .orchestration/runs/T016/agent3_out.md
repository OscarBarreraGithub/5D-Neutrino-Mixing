1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation OK: T016 imports no other constraint; physics path is via `zpole_lfv_etau`; diff for pre-existing T015/zpole_lfv/deltaf2/quarkConstraints/zpole_lfv is empty.
5. Contract/anchor OK: HARD; result numeric fields are floats; complex deltas only diagnostics; missing extra returns `predicted=None`, `ratio=None`, `passes=True`. Anchor spy used `pdg_or_equivalent.values[0]`, value `5e-06`; mismatch raises `AnchorError`.
6. Independent cross-check: `sm_total=3.649563333333334`; manual BR `3.508394829491628e-08`; constraint BR `3.508394829491628e-08`; abs diff `0`.
7. Deltas: `delta_L=0.0002309771776044444+9.239087104177778e-05j`; `delta_R=0+4.619543552088889e-05j`; norm `6.402054788893809e-08`; ratio `0.007016789658983256`.
8. Pass/fail probes: safe BR `4.677859914890927e-09`, ratio `0.0009355719829781852`, pass `True`; excluded BR `0.0001870793922231769`, ratio `37.41587844463538`, pass `False`.
9. Tests: targeted T016 `13 passed`; full `python -m pytest tests/constraints/ -q` `407 passed in 5.66s`; registry smoke `T016 top_higgs_ew HARD`, registry size `39`.
CODE-OK