Plan: inspect E004/E008 + E009.yaml, mirror INFO stub, verify tests.

Physics/parametrization: E009 records neutron-EDM-derived Weinberg bounds only, with Pospelov-Ritz and Haisch-Hala conventions from `E009.yaml`.

Reused vs built: reused scaffold `load_anchor`/registry pattern; added/verified Weinberg stub adapter boundary only, no direct physics imports in constraint.

Budget: default non-vetoing bookkeeping budget is `|C_6| < 1.2e-11 GeV^-2`, source `flavor_catalog/processes/edm_neutrino/E009.yaml` value `PDG2026-HaischHala2019:E009:c6_bound`.

Validation: independent YAML recompute checks `w` bound `4.1e-11 GeV^-2`, `C_6` bound `1.2e-11 GeV^-2`, neutron anchor `1.8e-26 e cm`.

NEEDS-HUMAN-PHYSICS: flagged both hadronic CP-odd gluonic matrix elements and RS CP-odd gluonic matching on `ParameterPoint`.

Files: `flavor_catalog_constraints/primary/edm_neutrino/E009.py`, `flavor_catalog_constraints/physics_adapters/weinberg_operator.py`, `tests/constraints/primary/edm_neutrino/test_E009.py`.

Tests: `test_E009.py` 9 passed; full `python -m pytest tests/constraints/ -q` 937 passed.