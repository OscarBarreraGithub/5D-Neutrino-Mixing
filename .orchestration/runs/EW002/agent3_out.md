1. BLOCKER `flavor_catalog_constraints/primary/top_higgs_ew/EW002.py:44`: EW002 bypasses scaffold `load_anchor`, using `load_full_yaml` + `build_anchor` custom lookup at `:130-174`. Fix by routing list entries through `load_anchor` like EW003/K005 and add a spy test.
2. BLOCKER `tests/constraints/primary/top_higgs_ew/test_EW002.py:150`: safe/excluded test calls `evaluate_first_row_sum`, the same adapter used by EW002 at `EW002.py:227`, so it does not exercise `Constraint.evaluate`. Fix by evaluating the constraint with controlled/mocked anchors and independent expected arithmetic.
3. Isolation/contract OK: EW002 imports no other constraint; physics path is only `physics_adapters.ckm_unitarity`; deltaf2 AST check found `changed_existing_functions=[]`; result numeric fields are real floats; empty point is non-crashing; deterministic `r1 == r2`.
4. Anchor probes: missing first-row entry, nonnumeric first-row value, and missing budget all raise `AnchorError`; however this is via the custom loader, not `load_anchor`.
5. Numerical cross-check: YAML first-row sum `0.99829999999999997`, ΔCKM `-0.0017000000000000348`, budget `0.00069999999999999999`, pull `2.4285714285714781`; constraint returns the same values, `passes=False`, `SOFT`.
6. Component arithmetic: `Vud^2=0.94803326890000006`, `Vus^2=0.050314976100000006`, `Vud^2+Vus^2=0.99834824500000008`, quoted-sum minus that `-4.8245000000113336e-05`.
7. Adapter safe/fail probe: `sum=1.0` gives pull `0` and pass; observed sum gives pull `2.4285714285714781` and fail.
8. Tests run: `test_EW002.py -q` = `7 passed`; `python -m pytest tests/constraints/ -q` = `125 passed in 2.97s`; registry smoke `EW002 top_higgs_ew SOFT`.

CODE-NEEDS-FIXES