Plan executed: studied scaffold/CR007/CR010/K001/K004, appended adapter helpers, added CR012, added tests, ran suite.

Physics: CR012 uses `m_spin1_diboson = kk_ew_mass_gev or M_KK` and `ratio = m_limit / m_proxy`.

Budget: `4.4 TeV` from `PDG2025:CR012:HVTB_Wprime_WZ_mass_lower` in `CR012.yaml` at `95% CL`.

Reused vs built: reused `collider_resonance` mass-limit core and scaffold `load_anchor`; built only CR012-specific adapter helpers.

Validation: independent core recomputation at `m_proxy = 5.0 TeV` gives `ratio = 4.4 / 5.0 = 0.88`, passing.

NEEDS-HUMAN-PHYSICS: documented in docstring and diagnostics for `sigma*BR`, branching surface, width, and acceptance recast.

Files changed: [collider_resonance.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/collider_resonance.py), [CR012.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/collider_rs/CR012.py), [test_CR012.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/collider_rs/test_CR012.py).

Pytest: focused `8 passed`; full `python -m pytest tests/constraints/ -q` `980 passed`.