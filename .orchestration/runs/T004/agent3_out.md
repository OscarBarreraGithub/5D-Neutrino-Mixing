1. BLOCKER: none. Isolation OK: `T004.py:51-60` imports only scaffold/base/registry and `physics_adapters.top_fcnc`; no staged/unstaged diffs in T001/T002/T003/top_fcnc/deltaf2.
2. SHOULD-FIX: none. Contract OK: missing couplings returns non-crashing `passes=True, predicted=None, ratio=None, budget=8.5e-06`; populated numeric fields are floats, complex dipoles stay in diagnostics.
3. NIT: none. Anchor OK: T004 limits route through `load_anchor` at `T004.py:304`; missing/mismatched probes raised `AnchorError`; active YAML limit is `ATLAS2023:T004:tugamma_left = 8.5e-06`.
4. CROSS-CHECK: direct `quarkConstraints.top_fcnc.photon_dipole_branching_fraction` recompute for `left=1+0.25j, right=0.3j` gives `core=6.0684210348288108e-06`, constraint `predicted=6.0684210348288108e-06`, `ratio=0.71393188645044836`, delta `0`.
5. PASS/FAIL: safe `left=1.0` -> `BR=5.265441244970768e-06`, `ratio=0.6194636758789139`, pass; excluded `left=2.0` -> `BR=2.1061764979883072e-05`, `ratio=2.4778547035156557`, fail.
6. DETERMINISM/REGISTRY: repeated `evaluate()` equality true, input matrices unchanged; registry smoke `all_constraints=39`, `has_T004=True`, `import_failures={}`.
7. TESTS: `python -m pytest tests/constraints/ -q` -> `407 passed in 5.55s`; `test_T004.py -q` -> `9 passed in 4.28s`.

CODE-OK