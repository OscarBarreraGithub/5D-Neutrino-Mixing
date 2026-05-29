1. NIT: Delta-F=2 checks are N/A for K005: this is ΔS=1 semileptonic `K_L -> pi0 nu nubar`, not an `M_12` observable. Code correctly uses only `Im(lambda_t X_t + X_NP)` at [rare_kaon_snd.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rare_kaon_snd.py:263), with no real/magnitude contamination.

2. NIT: QCD-running checklist is N/A: no four-quark `*_with_running` Delta-F=2 evaluator should be used here. K005 calls the rare-kaon adapter at [K005.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/kaon/K005.py:400), and diagnostics explicitly say no Delta-F=2 matrix element is evaluated at [K005.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/kaon/K005.py:452).

3. NIT: Budget is defensible: total BR is compared to the KOTO 90% CL upper limit `2.2e-9`, not a central residual. SM anchor is `2.94e-11`, so limit-minus-SM is `2.1706e-9`; using the full published upper limit is appropriate for a direct upper-bound veto. See [K005.yaml](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/processes/kaon/K005.yaml:107) and [K005.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/kaon/K005.py:302).

4. NIT: Anchor numbers match YAML/snapshots: PDG/KOTO `<2.2e-9`, BGS `2.59(29)e-11`, BV `(2.94±0.15)e-11`; audit confirms these at [factcheck_kaon.md](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/audits/factcheck_kaon.md:106).

5. NIT: SM-limit validation passes: code gives `BR_SM=2.95375989343059e-11`, `kappa_L=2.266439872701711e-10`, reproducing the requested `~3e-11`; targeted tests pass `9 passed`.

PHYSICS-OK