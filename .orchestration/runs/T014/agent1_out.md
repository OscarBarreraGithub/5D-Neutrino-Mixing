Plan executed: added adapter helper, implemented T014, added focused tests, ran full suite.

Physics: pure-NP HARD bound on `Z -> bs, bd, sd` FCNC hadronic-width branching proxy, using T014 yaml 95% CL rows from ECFA 2025 / Abu-Ajamieh 2026.

Reused/built: reused `quarkConstraints/zpole.py` via `physics_adapters/zpole.py`; added append-only down-sector off-diagonal Z helper and documented overlap proxy.

Budget: all three channels load `B < 2.9e-3` from `flavor_catalog/processes/secondary/top_higgs_ew/T014.yaml`.

Validation: zero FCNC coupling gives `BR = 0.0`; independent SM hadronic width-weight recomputation is `2.5225098333333333`.

Gap: flagged `NEEDS-HUMAN-PHYSICS` in docstring and diagnostics for missing rigorous RS FCNC-Z matching on `ParameterPoint`.

Files changed: `flavor_catalog_constraints/physics_adapters/zpole.py`, `flavor_catalog_constraints/secondary/top_higgs_ew/T014.py`, `tests/constraints/secondary/top_higgs_ew/test_T014.py`.

Tests: `tests/constraints/secondary/top_higgs_ew/test_T014.py` 11 passed; full `python -m pytest tests/constraints/ -q` 991 passed.