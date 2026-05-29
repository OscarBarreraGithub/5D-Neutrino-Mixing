Plan implemented: adapter, L001 constraint, focused tests.

Physics: `BR_NP = prefactor_br * |(Y_N_bar Y_N_bar^dagger)_{e mu}|^2 * (3 TeV / M_KK)^4`; SM set to `0.0`.

Reused: `flavorConstraints.muToEGamma.check_mu_to_e_gamma` / raw checker through `physics_adapters/lepton.py`.

Built: thin lepton adapter with explicit proxy input, flagged `NEEDS-HUMAN-PHYSICS`.

Budget: HARD `1.5e-13` from `flavor_catalog/processes/charged_lepton/L001.yaml:primary_current_limit` (MEG II 2025).

Validation: YAML-derived `C = 0.0019364916731037084`; independent proxy recomputation gives excluded `BR = 6.372660764462143e-12`, ratio `42.48440509641429`.

Files changed: `flavor_catalog_constraints/physics_adapters/lepton.py`, `flavor_catalog_constraints/primary/charged_lepton/L001.py`, `tests/constraints/primary/charged_lepton/test_L001.py`.

Pytest: `tests/constraints/primary/charged_lepton/test_L001.py` 9 passed; `tests/constraints/` 125 passed.