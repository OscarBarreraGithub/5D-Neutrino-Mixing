Plan executed: added L002 core, adapter, constraint, and tests.

Physics: `BR(mu->3e) = dipole off-shell factor + coherent Z/box vector-contact proxy`; SM set to `0.0`.

Reused: L001 dipole adapter. Built: `quarkConstraints/lfv_three_body.py` plus `physics_adapters/lfv_three_body.py`.

Budget: `1.0e-12`, loaded from `L002.yaml:primary_current_limit` PDG/SINDRUM.

Validation: independent test recomputation matches; dipole conversion factor `0.00612697641261286`.

Gap: Z/box RS matching is explicitly flagged `NEEDS-HUMAN-PHYSICS`.

Files changed: `quarkConstraints/lfv_three_body.py`, `flavor_catalog_constraints/physics_adapters/lfv_three_body.py`, `flavor_catalog_constraints/primary/charged_lepton/L002.py`, `tests/constraints/primary/charged_lepton/test_L002.py`.

Tests: `python -m pytest tests/constraints/ -q` -> `509 passed`.