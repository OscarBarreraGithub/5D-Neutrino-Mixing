Plan: reused Z-pole machinery, added charm adapter, implemented T012, added tests, ran constraints suite.

Physics: `R_c = Gamma_c/Gamma_had`, `A_c = (|g_L|^2-|g_R|^2)/(|g_L|^2+|g_R|^2)` from effective Zcc couplings; sources are `T012.yaml` PDG/LEP-SLC anchors.

Reused vs built: reused `quarkConstraints/zpole.py` via adapters; added `zpole_charm.py` for charm radiator calibration and documented Zcc proxy. T010 bottom functions untouched.

Budget: YAML one-sigma anchors: `R_c^0` budget `0.003`, `A_c` budget `0.027`, from `flavor_catalog/processes/top_higgs_ew/T012.yaml`.

SM validation: `R_c = 0.172100`, `A_c = 0.667577` (`A_FB^c = 0.073698` context).

Gap: RS Zcc shift remains `NEEDS-HUMAN-PHYSICS`; proxy uses charm-vs-up overlap non-universality because full EW KK/Z/Z' matching inputs are absent.

Files changed: `flavor_catalog_constraints/physics_adapters/zpole_charm.py`, `flavor_catalog_constraints/primary/top_higgs_ew/T012.py`, `tests/constraints/primary/top_higgs_ew/test_T012.py`.

Pytest: `test_T012.py` 10 passed; `python -m pytest tests/constraints/ -q` 275 passed.