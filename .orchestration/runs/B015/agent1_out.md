Plan executed: implemented B015 low-q2 inclusive path, adapter exposure, constraint, and tests.

Physics: active bin is `1 < q^2 < 6 GeV^2` from `B015.yaml`; LO inclusive C7/C9/C10 partonic shape normalized to YAML SM `BR[1,6]_mumu`.

Reused/built: reused `rare_b_dilepton` C9/C10 proxy and B011 C7 proxy; added append-only inclusive machinery plus `rare_b_meson` adapter wrappers.

Budget: HARD budget `9.067886552931955e-07` from B015 low-q2 exp + SM uncertainties plus 30% B016-style proxy-theory envelope.

SM validation: inclusive SM low-q2 prediction validates to `1.62e-06` from `B015.yaml`.

Gap: RS NP remains `NEEDS-HUMAN-PHYSICS` for full EW KK/Z/Z', lepton, NNLO/QED, charm-veto, dipole-loop, and covariance matching.

Files changed: `quarkConstraints/rare_b_dilepton.py`, `flavor_catalog_constraints/physics_adapters/rare_b_meson.py`, `flavor_catalog_constraints/primary/beauty/B015.py`, `tests/constraints/primary/beauty/test_B015.py`.

Tests: `python -m pytest tests/constraints/primary/beauty/test_B015.py -q` → 11 passed; `python -m pytest tests/constraints/ -q` → 296 passed.