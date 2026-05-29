1. NIT: Amplitude check is N/A: EW003 is a semileptonic CKM data-tension diagnostic, not ΔF=2 `M_12`; verdict uses `inclusive_exclusive_pull`, so no wrong real/imaginary `M_12` part enters. `flavor_catalog_constraints/primary/top_higgs_ew/EW003.py:387-402`, `flavor_catalog_constraints/physics_adapters/semileptonic_ckm.py:154-195`.

2. NIT: QCD running is N/A and correctly not used: no Wilson coefficients or `mu_had=2 GeV` matrix elements exist for this observable; verdict does not call any non-running ΔF=2 path. Running effect: N/A. `flavor_catalog_constraints/primary/top_higgs_ew/EW003.py:435-439`.

3. NIT: Budget is defensible and uncertainty-aware: `|V_cb|` pull = `2.4/sqrt(0.5^2+0.6^2)=3.073`; `|V_ub|` pull = `0.43/sqrt(0.258^2+0.156^2)=1.427`; compared to PDG `3.0 sigma`, this is a SOFT fail with ratio `1.024`. `EW003.py:387-402`, `EW003.yaml:83-140`, `EW003.yaml:201-210`.

4. NIT: Anchor numbers match snapshots/audit: PDG `42.2±0.5`, `39.8±0.6`, `4.13`, `3.70`; HFLAV `41.97±0.48`, `39.77±0.46`; FLAG `39.23±0.65`, `3.61±0.16`. `flavor_catalog/references/EW003/pdg_2024_vcb_vub.txt:9`, `hflav_2023_semileptonic_report.txt:10`, `flag_2024_b_semileptonic_lattice.txt:9`.

5. NIT: Severity, units, and diagnostics are physically consistent: `Severity.SOFT`, `NEEDS-HUMAN-PHYSICS`, no `ParameterPoint` RS proxy, units are dimensionless `10^-3` and `sigma` rather than GeV. `EW003.py:15-26`, `EW003.py:381-456`.

PHYSICS-OK