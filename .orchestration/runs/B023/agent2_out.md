1. NIT: Correct-amplitude Δm/CP check is N/A; B023 is a ΔB=1 branching ratio and correctly uses the complex Wilson response via `epsilon,eta` and `r_kstar=(1+1.31 eta) epsilon^2`: `quarkConstraints/rare_b_nunu.py:296`, `quarkConstraints/rare_b_nunu.py:305`; SM gives `epsilon=1`, `eta=0`, `r_kstar=1`.

2. NIT: QCD-running-to-2-GeV check is N/A, not a missing `*_with_running` call; this is semileptonic `b->s nu nubar`, not ΔF=2. The only short-distance multiplier is `eta_x=0.994`, giving `X_t=1.481` vs `X0=1.490` (-0.6%): `quarkConstraints/rare_b_nunu.py:85`, `quarkConstraints/rare_b_nunu.py:243`.

3. NIT: Budget is defensible for the YAML’s combined-limit construction: code uses total-BR 90% CL limit `2.7e-5`, not bare residual `2.7e-5 - 9.2e-6 = 1.78e-5`; SM ratio is `9.2e-6/2.7e-5=0.3407`: `flavor_catalog_constraints/primary/beauty/B023.py:339`, `flavor_catalog_constraints/primary/beauty/B023.py:356`.

4. NIT: Anchors match snapshots/YAML: Belle combined `<2.7e-5`, Buras SM `(9.2±1.0)e-6`, PDG/HFLAV neutral `<1.8e-5`, charged `<4.0e-5`: `flavor_catalog/processes/beauty/B023.yaml:133`, `flavor_catalog/processes/beauty/B023.yaml:145`, `flavor_catalog/references/B023/pdg2025_b0_kst0_nunu_api.txt:31`, `flavor_catalog/references/B023/pdg2025_bp_kstp_nunu_api.txt:31`.

5. NIT: Severity/units/notes are appropriate: HARD veto on dimensionless branching fraction, `NEEDS-HUMAN-PHYSICS` RS proxy is explicit, and Belle 2017 provenance is correctly used because YAML has no Belle II vector result: `flavor_catalog_constraints/primary/beauty/B023.py:388`, `flavor_catalog_constraints/primary/beauty/B023.py:466`, `flavor_catalog/processes/beauty/B023.yaml:204`.

PHYSICS-OK