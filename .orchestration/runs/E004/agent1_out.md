Plan: read scaffold/C003/K004/E004.yaml, add non-vetoing adapter + constraint, add focused tests, run constraints suite.

Implemented E004 as an INFO-only neutron EDM stub. It loads `canonical_limit` and cross-checks Abel 2020 nEDM@PSI `primary_measurement`; no hadronic or RS quark-dipole calculation is attempted.

Budget: `|d_n| < 1.8e-26 e cm` at 90% CL from `E004.yaml` PDG Live 2026 / Abel et al. PSI anchor.

Validation: measurement central `0.0 e cm`; ratio to limit `0.0`; total uncertainty recomputed as `1.1180339887498949e-26 e cm`.

NEEDS-HUMAN-PHYSICS flags are explicit for both neutron hadronic matrix elements and RS CP-odd quark EDM/CEDM/Weinberg matching.

Files changed: [neutron_edm.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/neutron_edm.py), [E004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/edm_neutrino/E004.py), [test_E004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/edm_neutrino/test_E004.py).

Tests: `python -m pytest tests/constraints/primary/edm_neutrino/test_E004.py -q` -> 8 passed; `python -m pytest tests/constraints/ -q` -> 769 passed.