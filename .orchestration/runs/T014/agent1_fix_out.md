- Core helper now lives in [quarkConstraints/zpole.py:493](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/zpole.py:493); BR helper starts at [quarkConstraints/zpole.py:577](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/zpole.py:577).
- Adapter is now thin wrappers only at [zpole.py:176](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/zpole.py:176).
- T014 formula/diagnostics updated at [T014.py:11](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:11) and [T014.py:331](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:331).
- New formula: `B = Gamma_FCNC / (Gamma_Z,total^SM + Gamma_FCNC)`.
- Corrected coupling for `B < 2.9e-3`: `sqrt(|delta_g_L|^2+|delta_g_R|^2) = 0.0420605`.
- Confirmed old tightening: core `Gamma_had/Gamma_Z = 2.5225098/3.6495633 = 0.691`, the same effect as the physical `~0.699` ratio.
- Tests updated to recompute from `quarkConstraints.zpole` core at [test_T014.py:180](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/secondary/top_higgs_ew/test_T014.py:180) and [test_T014.py:300](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/secondary/top_higgs_ew/test_T014.py:300).
- `python -m pytest tests/constraints/secondary/top_higgs_ew/test_T014.py -q`: 11 passed.
- `python -m pytest tests/constraints/ -q`: 1002 passed.
FIX-DONE