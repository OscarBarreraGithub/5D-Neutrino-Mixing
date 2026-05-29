Plan: studied scaffold/gold examples/T005.yaml, appended adapter proxy, added T005 + tests, ran suite.

Physics: `BR(t -> c g)` via top-FCNC chromomagnetic dipole, pure-NP ratio `BR_NP / BR_limit`.

Reused/built: reused `quarkConstraints.top_fcnc` gluon width; built append-only adapter proxy `t_to_q_gluon_from_couplings`.

Budget: `3.7e-4` from `PDG2025:T005:ATLAS2022:t_cg` in `T005.yaml`.

SM/validation: SM anchor `4.6e-12`; manual proxy check gives `BR_NP=1.118111915877719e-4`, ratio `0.3021924096966808`.

NEEDS-HUMAN-PHYSICS: RS chromodipole matching/prod-recast remains a documented proxy gap in docstring and diagnostics.

Files changed: `physics_adapters/top_fcnc.py`, `primary/top_higgs_ew/T005.py`, `test_T005.py`.

Tests: focused T005 `9 passed`; full `python -m pytest tests/constraints/ -q` `429 passed`.