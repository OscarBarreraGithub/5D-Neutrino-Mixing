Plan: inspected scaffold/K001/K004/L001/L002, implemented core + adapter + L003, added tests, ran full constraints suite.

Physics: KKO-style `CR = omega_conv / Gamma_capture`, with dipole + scalar/vector proton/neutron overlap terms.

Sources: L003 YAML anchors plus KKO nuclear inputs, https://arxiv.org/abs/hep-ph/0203110.

Reused vs built: reused L001 dipole adapter; built `quarkConstraints.mu_e_conversion` and `physics_adapters.mu_e_conversion`.

Budget: `6.7e-17` from Mu2e Al expected 90% C.L. projection in L003.yaml; SINDRUM II Au `7.0e-13` is loaded as diagnostic, not treated as an Al limit.

Validation: independent test recomputation gives `CR = 1.9778911811395536e-18`, ratio `2.9521e-02`; SM set to `0.0`.

NEEDS-HUMAN-PHYSICS: scalar/vector RS matching and unknown dipole-contact phase/chirality remain explicit documented proxies.

Files: `quarkConstraints/mu_e_conversion.py`, `flavor_catalog_constraints/physics_adapters/mu_e_conversion.py`, `flavor_catalog_constraints/primary/charged_lepton/L003.py`, `tests/constraints/primary/charged_lepton/test_L003.py`.

Pytest: `python -m pytest tests/constraints/ -q` -> `626 passed in 15.37s`.