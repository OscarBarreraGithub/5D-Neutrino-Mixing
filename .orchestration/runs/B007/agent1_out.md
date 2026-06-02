Plan: study scaffold/examples, reuse rare-B dilepton core via adapter, add secondary B007 tests, run full constraints suite.

Physics: `B_q -> e+e-` uses Buras `b -> q ll` C10-dominant leptonic formula with charged-lepton mass set to `m_e`; sources are `B007.yaml` PDG 2026 limits and Bobeth et al. SM rows.

Reused vs built: reused `quarkConstraints/rare_b_dilepton.py`; B007 uses the append-only `physics_adapters/rare_b_electronic.py` electron wrapper. No muon functions modified.

Budget: pure-NP HARD upper-limit budgets from `B007.yaml`: `BR(B_s -> ee) < 9.4e-9`, `BR(B0 -> ee) < 2.5e-9` at 90% CL.

SM validation: independent formula test gives `B_s -> ee = 8.540154879838497e-14`, `B_d -> ee = 2.450169556007537e-15`.

Gap: RS contribution remains `NEEDS-HUMAN-PHYSICS`; current implementation uses documented C9/C10 proxy without full EW KK/Z/Z'/electron/scalar matching.

Files: `flavor_catalog_constraints/secondary/beauty/B007.py`, `flavor_catalog_constraints/physics_adapters/rare_b_electronic.py`, `tests/constraints/secondary/beauty/test_B007.py`.

Pytest: targeted B007 `11 passed`; full `python -m pytest tests/constraints/ -q` -> `959 passed in 19.25s`.