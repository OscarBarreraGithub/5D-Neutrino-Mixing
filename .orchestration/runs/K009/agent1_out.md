Plan: reuse K008 y7V/y7A machinery, add K009-only adapter wrapper, add K009 constraint/tests, run full constraint suite.

Physics: K009 uses ISU/BMS `K_L -> pi0 l+l-` direct/indirect/interference/CPC decomposition with K009 muon coefficients from `K009.yaml`.

Reused vs built: reused rare_kaon_dilepton K008 y7V/y7A Wilson proxy; built append-only `rare_kaon_dilepton_muon.py` adapter for muon phase-space coefficients.

Budget: `3.8e-10`, PDG/KTeV 90% CL upper limit from `flavor_catalog/processes/kaon/K009.yaml`.

SM validation: direct CP `2.0151243040372812e-12`; constructive total `1.5268661725244892e-11`, within YAML `(1.5 +/- 0.3)e-11`.

NEEDS-HUMAN-PHYSICS: RS y7V/y7A proxy is reused for muons; full EW KK/Z/Z' and muon-sector matching are absent on `ParameterPoint`.

Files changed: `flavor_catalog_constraints/physics_adapters/rare_kaon_dilepton_muon.py`, `flavor_catalog_constraints/primary/kaon/K009.py`, `tests/constraints/primary/kaon/test_K009.py`.

Pytest: `python -m pytest tests/constraints/ -q` -> `352 passed`.