Plan executed: traced L001, added a tau-mu dipole adapter, added `L007`, then tested.

Physics: pure-NP `BR(tau -> mu gamma) = prefactor_br * |(Ybar Ybar^dagger)_{mu tau}|^2 * (3 TeV / M_KK)^4`; SM set to `0.0` when evaluated.

Reused vs built: reused `flavorConstraints.muToEGamma` plus L001 proxy validation; built only `physics_adapters/lepton_tau_mu.py` to map `(mu,tau)` via row permutation.

Budget: `4.2e-8` from `L007.yaml` `primary_current_limit` / Belle 2021 current limit. Belle II `5.0e-9` projection is diagnostic only, not the HARD veto.

Validation: manual recomputation for `y=(0.10,0.20,0.30)` gives `lhs=0.02103677462019741`, `BR=1.7701835456839274e-11`, ratio `4.2147e-4`.

NEEDS-HUMAN-PHYSICS: still flagged because full tau-mu lepton-sector RS loop matching/couplings are not on `ParameterPoint`.

Files added: `flavor_catalog_constraints/physics_adapters/lepton_tau_mu.py`, `flavor_catalog_constraints/primary/charged_lepton/L007.py`, `tests/constraints/primary/charged_lepton/test_L007.py`.

Tests: `test_L007.py` `11 passed`; full `python -m pytest tests/constraints/ -q` `841 passed`.