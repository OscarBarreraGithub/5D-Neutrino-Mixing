Implemented C005.

- Plan executed: studied scaffold/C004/C005 YAML, added C005, added focused tests, ran suite.
- Physics: `D0 -> e+e-` uses shared `c -> u l l` Hamiltonian with `C10 - C10'` and electron mass, SM SD set by core to `0.0`; LD context is not subtracted.
- Reused: `flavor_catalog_constraints.physics_adapters.rare_charm_dilepton.d0_ee_from_couplings`; no C004 muon functions or core physics changed.
- Budget: `7.9e-8` from `C005.yaml` `canonical_current_limit` PDG Live/API S032.39, HARD upper-limit ratio.
- Validation: independent manual proxy recomputation gives `BR = 1.333296784847654e-17` for the test coupling; helicity suppression vs muon is checked.
- Gap: diagnostics and docstring flag `NEEDS-HUMAN-PHYSICS` for the reused RS Z/KK-penguin proxy.
- Files changed: [C005.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charm/C005.py:1), [test_C005.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/charm/test_C005.py:1).
- Tests: `python -m pytest tests/constraints/ -q` -> `265 passed`.