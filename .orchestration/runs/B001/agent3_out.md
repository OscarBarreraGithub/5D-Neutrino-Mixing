1. SHOULD-FIX `flavor_catalog_constraints/primary/beauty/B001.py:61` / `tests/constraints/primary/beauty/test_B001.py:26`: SM theory sigma `0.062` is hardcoded and duplicated in tests, so the full veto budget is not YAML/anchor-driven. Fix: store/load the HPQCD SM uncertainty from `B001.yaml` and test missing/mismatched failure.
2. NIT `flavor_catalog_constraints/physics_adapters/deltaf2.py:77`: adapter diff also contains non-B001 D0/Bs/K helpers. Existing function bodies are unchanged and `git diff -- quarkConstraints/` is empty; split unrelated helpers if strict B001 isolation is required.

Isolation/contract: B001 imports no other constraint; missing extra returns non-crashing result; numeric fields are floats, complex values stay in diagnostics; deterministic/no mutation test passes.
Anchor probes: missing experimental candidate -> `AnchorError`; mismatched `DELTA_M_BD_EXP_GeV` -> `AnchorError`.
Numerical safe: predicted/direct `2.17789186924378677e-14`, budget `3.36989977945936853e-14`, ratio `0.6462779345898485`, pass=True, agree.
Numerical excluded: predicted/direct `5.44472967310946758e-13`, budget `3.36989977945936853e-14`, ratio `16.156948364746213`, pass=False, agree.
Tests: B001 `7 passed`; `python -m pytest tests/constraints/ -q` -> `51 passed`; registry smoke `B001 beauty HARD Delta m_d`.
CODE-NEEDS-FIXES