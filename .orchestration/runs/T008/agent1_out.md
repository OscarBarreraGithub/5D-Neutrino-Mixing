Plan executed: studied scaffold, K001/K004, T007, and `T008.yaml`; implemented and tested T008.

Physics: `BR(t -> H u)` uses the shared top-FCNC Higgs scalar-Yukawa width with `L = -H ubar (y_L P_L + y_R P_R)t + h.c.`; SM treated as negligible pure-NP bound.

Reused vs built: reused `physics_adapters.top_fcnc` Higgs-FCNC machinery; built only new append-only `T008.py` and its test. No T007 edits.

Budget: `1.9e-4` from `PDG2026:T008:tHu_headline` in `flavor_catalog/processes/top_higgs_ew/T008.yaml`.

Validation: zero-yukawa SM proxy is `0.0`; independent core recomputation fixture gives `BR = 3.4762570502439738e-6`.

NEEDS-HUMAN-PHYSICS: RS top-Higgs FCNC matching remains a documented scalar-Yukawa proxy, flagged in docstring and diagnostics.

Files changed: `flavor_catalog_constraints/primary/top_higgs_ew/T008.py`, `tests/constraints/primary/top_higgs_ew/test_T008.py`.

Tests: focused T008 `9 passed`; full `python -m pytest tests/constraints/ -q` `745 passed`.