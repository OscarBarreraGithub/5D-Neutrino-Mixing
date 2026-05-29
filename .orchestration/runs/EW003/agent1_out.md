Plan: read scaffold/predecessors, implement EW003 as a SOFT data constraint, add focused tests, run full constraint suite.

Physics: direction-aware Gaussian pull, `|Vxb_inc - Vxb_exc| / sqrt(sigma_inc^2 + sigma_exc^2)`, with YAML uncertainty components. Sources are `EW003.yaml` PDG/HFLAV/FLAG anchors loaded via virtual `load_anchor`.

Built: new `semileptonic_ckm` statistical adapter; no direct physics-module import from constraint.  
Budget: `3.0 sigma` from `PDG 2024 |V_cb| inclusive-exclusive marginal consistency`.

Validation: PDG `|V_cb|` pull `3.072885`, PDG `|V_ub|` pull `1.427005`; max ratio `1.024295`, so EW003 reports `passes=False` with `Severity.SOFT`.

NEEDS-HUMAN-PHYSICS: flagged in diagnostics; no RS NP proxy applied because charged-current EW/WET/covariance inputs are not on `ParameterPoint`.

Files changed: `flavor_catalog_constraints/physics_adapters/semileptonic_ckm.py`, `flavor_catalog_constraints/primary/top_higgs_ew/EW003.py`, `tests/constraints/primary/top_higgs_ew/test_EW003.py`.

Pytest: `python -m pytest tests/constraints/ -q` -> `118 passed in 2.90s`.