Plan: inspect scaffold/T005/T006 YAML -> implement append-only T006 -> add tests -> run full constraints suite.

Physics: `BR(t -> u g)` via shared top-FCNC chromomagnetic dipole convention, pure-NP vs collider bound; sources loaded from `T006.yaml`.

Reused vs built: reused `physics_adapters.top_fcnc`; no T005 functions or physics modules modified.

Budget: `2.0e-5` from `CMS2017:T006:t_ug`; ATLAS `6.1e-5` retained as context.

Validation: SM anchor `3.6e-14`; independent proxy-width recomputation gives `BR_NP=1.118111915877719e-4` for the test point.

Gap: RS chromo-dipole matching remains `NEEDS-HUMAN-PHYSICS` in docstring and diagnostics.

Files changed: [T006.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T006.py), [test_T006.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/top_higgs_ew/test_T006.py).

Tests: `test_T006.py` 9 passed; `python -m pytest tests/constraints/ -q` 664 passed.