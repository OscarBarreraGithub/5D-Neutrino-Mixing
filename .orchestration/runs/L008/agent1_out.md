Plan executed: reused L001/L007 dipole pattern, added a tau->e row-permutation adapter, wired L008, tested.

Physics: `BR_NP = 4.0e-8 * |(Y_N_bar Y_N_bar^dagger)_{e tau}|^2 * (3000 GeV / M_KK)^4`, `SM=0`.

Reused vs built: reused `flavorConstraints.muToEGamma.check_mu_to_e_gamma_raw`; built only append-only `lepton_tau_e.py`.

Budget: HARD `BR(tau -> e gamma) < 3.3e-8` from `L008.yaml` PDG/BaBar current limit; Belle II `1.2e-8` is carried as projection diagnostics, not the veto.

Validation: tau->e rotation `(0, 2, 1)` gives `|A_e_tau|=0.03365883939`, `BR=4.53166987695e-11`, ratio `1.373e-3` for the test point.

NEEDS-HUMAN-PHYSICS: tau-e dipole proxy is flagged; missing/invalid lepton inputs return `evaluated=False`.

Files changed: `flavor_catalog_constraints/physics_adapters/lepton_tau_e.py`, `flavor_catalog_constraints/primary/charged_lepton/L008.py`, `tests/constraints/primary/charged_lepton/test_L008.py`.

Tests: `test_L008.py` 11 passed; charged-lepton 105 passed; full `tests/constraints/` 898 passed. Unrelated dirty/untracked worktree entries were left untouched.