Plan completed: studied scaffold/E004/K004/E008.yaml, added E008 stub + adapter, added tests, ran constraints suite.

Physics: E008 records YAML qCEDM translation anchors only: neutron and Hg-derived bounds from Pospelov/Ritz 2001 and Olive/Pospelov/Ritz/Santoso 2005, backed by PDG neutron EDM and Graner Hg anchors.

Machinery: reused scaffold anchor helpers; built new small adapter [quark_cedm.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/quark_cedm.py).

Budget: non-vetoing INFO budget is `1.1e-27 cm` from `E008.yaml qcedm_translations.mercury_isovector_combination`.

Validation: independent YAML recompute gives neutron `1.636e-26 cm -> 1.6e-26`, Hg `1.057e-27 cm -> 1.1e-27`.

NEEDS-HUMAN-PHYSICS: flagged both missing hadronic/nuclear matrix elements and missing RS CP-odd qCEDM matching; no ParameterPoint inputs used.

Files changed: [E008.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/edm_neutrino/E008.py), [quark_cedm.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/quark_cedm.py), [test_E008.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/edm_neutrino/test_E008.py).

Tests: targeted E008 `7 passed`; full `python -m pytest tests/constraints/ -q` -> `810 passed`.